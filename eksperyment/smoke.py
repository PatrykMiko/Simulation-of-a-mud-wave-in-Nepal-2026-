"""Small synthetic site test of both release modes and the private solver."""
import shlex
import tempfile
from .build import PREFIX
from pathlib import Path
import numpy as np
from ravaflow_tools import synthetic, nepal as N, runner
from nepal_gradual.simulation import model_config as gradual
from .design import BASE
from .run import config,simulate
from .render import render


def verify():
    with tempfile.TemporaryDirectory(prefix='nepal20-test-') as tmp:
        folder=Path(tmp)
        cfg=config()
        cfg.update(data_dir=str(folder/'data'),study_dir=str(folder/'work'),volume_m3=200000,
                   output_interval_s=5,control_points=[])
        data=Path(cfg['data_dir'])/'copernicus'
        synthetic.generate(data,volume=cfg['volume_m3'],nrows=96,ncols=96,xll=339990,yll=3110010)
        grid=N.prepare(cfg,1)
        for mode in ('instant','gradual'):
            p={**BASE,'mode':mode,'release_duration_s':4}
            result=simulate(cfg,p,folder/mode,grid,cap=20,outlet=(340750,3111350))
            assert result['time_s']>=20 and result['end_reason']=='safety_cap'
            if mode=='gradual':
                render(folder/mode,data)
                # Compare the actual native grid to the unmodified installed solver.
                ref=folder/'reference'
                ref.mkdir()
                config_ref=gradual(cfg,p,ref,grid,N.ensure_mapset(ref),20)
                config_ref['pre_commands']=['export PATH='+shlex.quote(str(PREFIX/'reference'))+':'+shlex.quote(str(PREFIX))+':"$PATH"']
                standard=runner.run(config_ref)
                from ravaflow_tools import metrics as M
                actual,_=M.read_raster(folder/mode/'hflow_max.tif')
                expected,_=runner.read_output_ascii(Path(standard['ascii_dir'])/'nepal_hflow_max.asc')
                if not np.allclose(actual,expected,atol=.001,rtol=.001):
                    raise ValueError(f'Różnica względem oryginalnego solvera: {np.max(abs(actual-expected))}')
    print('Test syntetyczny zakończony.')


if __name__=='__main__':
    verify()
