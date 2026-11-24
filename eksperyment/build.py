"""Build a private Rev.7 wrapper and core from the verified upstream archive."""
import hashlib
from pathlib import Path
import shutil
import subprocess
import zipfile

ROOT = Path(__file__).resolve().parents[1]
PREFIX = ROOT / '.local-tools/nepal20'
SHA = '8301b4804f45e228633aa823caba4a3387a47f8bd143cf5569fdc5a9f95ad3a7'
URL = 'https://www.landslidemodels.org/r.avaflow/software/r.avaflow.40G.zip'


def replace(text, before, after):
    if text.count(before) != 1:
        raise ValueError('Nieznany wariant źródła Rev.7')
    return text.replace(before, after)


def build():
    PREFIX.mkdir(parents=True, exist_ok=True)
    archive = PREFIX/'rev7.zip'
    if not archive.exists():
        local = ROOT/'vendor/r.avaflow.40G-rev7.zip'
        if local.exists():
            shutil.copyfile(local, archive)
        else:
            tmp = archive.with_suffix('.download')
            subprocess.run(['curl','-fL','--retry','2',URL,'-o',tmp],check=True)
            tmp.replace(archive)
    if hashlib.sha256(archive.read_bytes()).hexdigest() != SHA:
        raise ValueError('Niepoprawna suma SHA-256 źródeł Rev.7')
    with zipfile.ZipFile(archive) as z:
        base = 'r.avaflow.40G/'
        code = z.read(base+'r.avaflow.main/main.c').decode().replace('\r\n','\n')
        wrapper = z.read(base+'r.avaflow/r.avaflow.40G.py').decode()
        (PREFIX/'LICENSE').write_bytes(z.read(base+'LICENSE'))
        for module in ('r.avaflow.background','r.avaflow.paraview','r.avaflow.mult'):
            path = PREFIX/module
            path.write_bytes(z.read(base+module+'/'+module+'.py'))
            path.chmod(0o755)
    upstream = code
    anchor = '    while ( tsum <= tmax && ccontinue == 1 )'
    code = replace(code,anchor,'''    if (!getenv("NEPAL_OUTLET_X") || !getenv("NEPAL_OUTLET_Y")) return 71;
    const double outlet_x = atof(getenv("NEPAL_OUTLET_X")), outlet_y = atof(getenv("NEPAL_OUTLET_Y"));
'''+anchor)
    code = replace(code, '|| ccontinue == 0 ) { // if defined interval',
                   '|| ccontinue == 0 || tsum >= tmax ) { // if defined interval')
    anchor = '        if ( (int)(round( 1000 * tint )) >= (int)( round( 1000 * tout )) || ccontinue == 0 || tsum >= tmax )'
    code = replace(code, anchor, Path(__file__).with_name('dynamics.c.inc').read_text()+'\n'+anchor)
    anchor = '// -- START -- Writing hydrograph infos to files'
    assert code.count(anchor) == 2
    pos = code.rindex(anchor)
    code = code[:pos]+Path(__file__).with_name('monitor.c.inc').read_text()+'\n'+code[pos:]
    # Sparse frames replace the per-step heavyweight visualization exports.
    anchor = 'if ( nout == 1 ) jmin = 0; else jmin = 1;'
    a = code.index('// -- START -- Writing csv files for virtual reality')
    b = code.index(anchor,a)
    code = code[:a]+replace(code[a:b],'if ( sico.MULT == 0 )','if (0)')+code[b:]
    a = code.index('// -- START -- Preparing and writing output raster maps and velocity fields')
    b = code.index('// -- STOP --- Preparing and writing output raster maps and velocity fields',a)
    code = code[:a]+'if (!ccontinue || tsum >= tmax) {\n'+code[a:b]+'}\n'+code[b:]
    main = PREFIX/'main.c'
    main.write_text(code)
    grass = Path(subprocess.check_output(['grass','--config','path'],text=True).strip())
    reference = PREFIX/'reference'
    reference.mkdir(exist_ok=True)
    original = reference/'main.c'
    original.write_text(upstream)
    for dest, main in ((PREFIX, main), (reference, original)):
        binary = dest/'r.avaflow.main'
        subprocess.run(['cc','-O2','-fno-fast-math','-DPACKAGE="grassmods"',
                        '-I'+str(grass/'include'),'-I/usr/include/gdal',main,'-L'+str(grass/'lib'),
                        '-Wl,-rpath,'+str(grass/'lib'),'-lgrass_gis','-lgrass_segment','-lgrass_raster',
                        '-lm','-o',binary],check=True)
        path = dest/'r.avaflow.40G'
        path.write_text(replace(wrapper,'grass.run_command("r.avaflow.main")',
                                 'grass.run_command('+repr(str(binary))+')'))
        path.chmod(0o755)


if __name__ == '__main__':
    build()
