"""Regression for constrained IoU ranking, batch barriers and retention."""
import csv
import json
from pathlib import Path
from zipfile import ZipFile

from eksperyment import study as S, cli
from ravaflow_tools.checkpoint import atomic, read
from eksperyment import objective as O


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
        summaries=[dict(id=O.load()['points'][0]['id'],arrival_s=350 if target else 5000,
                        speed_m_s=22.,speed_samples=13,window_complete=True)]
        atomic(folder/'report.json',dict(number=number,params=params,iou_all=iou,
              **O.score(iou,summaries,O.load()),
              metrics={'all':{'jaccard':iou},'train':{'jaccard':float(target)}},
              reached=target,time_s=120,arrival_s=120 if target else -1,end_reason='safety_cap'))
        (folder/'raster.tif').write_bytes(str(number).encode())
    monkeypatch.setattr(S,'launch',worker)
    S.study(root)
    assert history==[(n,list(range(1,n))) for n in (1,6,11,16)]
    assert read(root/'finished.json')['winner']==19
    assert read(root/'finished.json')['dynamics_pass']
    assert sorted(p.name for p in (root/'artifacts').iterdir())==['t19']
    with ZipFile(root/'paczka.zip') as z:
        assert z.testzip() is None
        assert z.read('artifacts/t19/raster.tif')==b'19'
        assert len([p for p in z.namelist() if p.startswith('records/')])==20
        rows=list(csv.DictReader(z.read('ranking.csv').decode().splitlines()))
        assert rows[0]['number']=='19'
        assert rows[-1]['number']=='12'
        p=json.loads(z.read('records/01.json'))['params']
        assert p==read(root/'plans/01.json')
        assert p['solid_fraction']==0.6373316047825536
        instant = {k:v for k,v in p.items() if k!='release_duration_s'}
        assert read(root/'plans/04.json')=={**instant,'mode':'instant','solid_fraction':.5}
        assert read(root/'plans/09.json')=={**read(root/'plans/04.json'),'delta':6.,'manning':.022,'entrainment_log10':-9.}
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


def test_rolling_save_survives_failed_replacement_and_restores_files(tmp_path,monkeypatch):
    from types import SimpleNamespace
    from eksperyment import autosave as A
    import pytest
    root=tmp_path/'wyniki'
    signature={'root':'test'}
    atomic(root/'checkpoint/state.json',dict(state='saved',saved_at=1,generation=1,
                                             files=[],signature=signature))
    atomic(root/'work/t01/state.json',{'t':120})
    monkeypatch.setattr(A.os,'sync',lambda:None)
    monkeypatch.setattr(A.shutil,'disk_usage',lambda p:SimpleNamespace(free=10**12))
    old=A.publish(root)
    atomic(root/'work/t01/state.json',{'t':240})
    with monkeypatch.context() as patch:
        def disk_full(*args):
            raise OSError('disk full')
        patch.setattr(A,'copy',disk_full)
        with pytest.raises(OSError,match='disk full'):
            A.publish(root)
    assert A.latest(root)==old
    A.restore(root,signature)
    assert read(root/'work/t01/state.json')=={'t':120}
    atomic(root/'work/t01/state.json',{'t':360})
    assert A.publish(root)['slot']!=old['slot']
    assert not (A.store(root)/f"copy-{old['slot']}").exists()
    with pytest.raises(RuntimeError,match='Środowisko'):
        A.restore(root,{'root':'moved'})
    assert read(root/'work/t01/state.json')=={'t':360}


def test_vm_shutdown_requires_confirmed_save(tmp_path, monkeypatch):
    import os
    import subprocess
    executable=tmp_path/'multipass'
    executable.write_text('#!/bin/sh\necho "$*" >> "$VM_TEST_LOG"\n'
                          'if [ "$1" = exec ]; then exit "$VM_TEST_RESULT"; fi\n')
    executable.chmod(0o755)
    log=tmp_path/'log'
    monkeypatch.setenv('PATH',str(tmp_path)+os.pathsep+os.environ['PATH'])
    monkeypatch.setenv('VM_TEST_LOG',str(log))
    wrapper=Path(__file__).resolve().parents[1]/'scripts/eksperyment-vm.sh'
    monkeypatch.setenv('VM_TEST_RESULT','1')
    assert subprocess.run(['bash',wrapper,'save-stop']).returncode!=0
    assert 'stop nepal30' not in log.read_text().splitlines()
    monkeypatch.setenv('VM_TEST_RESULT','0')
    subprocess.run(['bash',wrapper,'save-stop'],check=True)
    assert log.read_text().splitlines()[-1]=='stop nepal30'


def test_supervisor_triggers_automatic_save_by_time_or_batch(tmp_path,monkeypatch):
    from types import SimpleNamespace
    from eksperyment import autosave as A

    monkeypatch.setattr(cli.C, 'alive', lambda _:True)
    monkeypatch.setattr(cli,'SESSION',SimpleNamespace(state=lambda:{'state':'running'}))
    monkeypatch.setattr(cli.disk,'members',lambda root:[])
    monkeypatch.setattr(cli.disk,'estimate',lambda *args:cli.GIB)
    monkeypatch.setattr(cli.shutil,'disk_usage',lambda _:SimpleNamespace(free=100*cli.GIB))
    requests=[]
    def save_once(required):
        requests.append(required)
        return False
    monkeypatch.setattr(cli,'periodic_save',save_once)
    for case in ('hour','batch'):
        root=tmp_path/case
        root.mkdir()
        monkeypatch.setattr(cli,'OUT',root)
        clock=[0]
        def tick(_):
            assert clock[0]==0, 'Nadzór przegapił pierwszy zaplanowany zapis.'
            if case=='batch':
                for n in range(5):
                    atomic(root/'records'/f'{n:02}.json',{})
                clock[0]=1
            else:
                clock[0]=A.PERIOD
        monkeypatch.setattr(cli,'time',SimpleNamespace(monotonic=lambda:clock[0],sleep=tick))
        cli.supervise()
    assert requests==[cli.GIB,cli.GIB]
