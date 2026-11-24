"""Regression for the requested full-domain ranking, batch barriers and retention."""
import csv
import json
from pathlib import Path
from zipfile import ZipFile

from eksperyment import study as S, cli
from ravaflow_tools.checkpoint import atomic, read


def test_campaign(tmp_path,monkeypatch):
    root=tmp_path/'wyniki'
    cfg={'data_dir':str(tmp_path),'dem':'copernicus'}
    grid=tmp_path/'grid'
    atomic(grid/'manifest.json',{})
    monkeypatch.setattr(S,'config',lambda:cfg)
    monkeypatch.setattr(S.N,'prepare',lambda *a:grid)
    monkeypatch.setattr(cli,'signature',lambda:{})
    monkeypatch.setattr(S,'render',lambda folder,data:(folder/'fala.mp4').write_bytes(b'movie'))
    history=[]
    proposer=S.batch
    def batch(first,reports):
        history.append((first,[r['number'] for r in reports]))
        return proposer(first,reports)
    monkeypatch.setattr(S,'batch',batch)
    def worker(number):
        params=read(root/'plans'/f'{number:02}.json')
        iou=0.8 if number in (12,19) else .25
        folder=root/'work'/f't{number:02}'
        target=number!=12
        atomic(folder/'report.json',dict(number=number,params=params,iou_all=iou,
              metrics={'all':{'jaccard':iou},'train':{'jaccard':float(target)}},
              reached=target,time_s=120,arrival_s=120 if target else -1,end_reason='safety_cap'))
        (folder/'raster.tif').write_bytes(str(number).encode())
    monkeypatch.setattr(S,'launch',worker)
    S.study(root)
    assert history==[(n,list(range(1,n))) for n in (1,6,11,16)]
    assert read(root/'finished.json')['winner']==12
    assert sorted(p.name for p in (root/'artifacts').iterdir())==['t12']
    with ZipFile(root/'paczka.zip') as z:
        assert z.testzip() is None
        assert z.read('artifacts/t12/raster.tif')==b'12'
        assert len([p for p in z.namelist() if p.startswith('records/')])==20
        rows=list(csv.DictReader(z.read('ranking.csv').decode().splitlines()))
        assert [row['number'] for row in rows[:2]]==['12','19']
        p=json.loads(z.read('records/01.json'))['params']
        assert p==read(root/'plans/01.json')
        assert p['solid_fraction']==0.6373316047825536
        instant = {k:v for k,v in p.items() if k!='release_duration_s'}
        assert read(root/'plans/04.json')=={**instant,'mode':'instant','solid_fraction':.5}
        assert read(root/'plans/09.json')=={**read(root/'plans/04.json'),'manning':.025}
        assert all(read(plan).get('release_duration_s',0)<=300 for plan in (root/'plans').glob('*.json'))
        assert all(read(plan)['solid_fraction']<=p['solid_fraction'] for plan in (root/'plans').glob('*.json'))


def test_monitor_stop_rules(tmp_path):
    """Exercise the actual monitor with controlled depth histories and geometry."""
    import os
    import subprocess

    src=Path(__file__).resolve().parents[1]/'eksperyment/monitor.c.inc'
    program=tmp_path/'monitor.c'
    program.write_text('''
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
struct { int IMAX, M, N, MODEL; double BDWEST, BDNORTH, CSZ; } sico={8,1,8,7,0,30,30};
int main(int argc,char **argv) {
    int cdomain2[8]={0}, px[8]={0}, py[8]={0,1,2,3,4,5,6,7};
    float awt[8][33]={{0}};
    double outlet_x=atof(argv[1]),outlet_y=15,tmax=43200;
    int ccontinue=1,csuccess=1;
    for (int nout=1;nout<400 && ccontinue;nout++) {
        double tsum=nout*120;
        for(int i=0;i<8;i++) awt[i][16]=awt[i][32]=1;
#include "monitor.c.inc"
    }
    return 0;
}
''')
    binary=tmp_path/'monitor'
    subprocess.run(['cc','-I',str(src.parent),str(program),'-lm','-o',str(binary)],check=True)
    mask=tmp_path/'valid'
    mask.write_bytes(bytes([1])*8)
    for name,x,last_time in [('outlet_stable','150',720),('stalled','4000',10800)]:
        d=tmp_path/name
        d.mkdir()
        subprocess.run([binary,x],check=True,env={**os.environ,'NEPAL_VALID':str(mask),'NEPAL_FRAME_DIR':str(d)})
        last=(d/'timeline.csv').read_text().splitlines()[-1].split(',')
        assert last[-1]==name
        assert float(last[1])==last_time
