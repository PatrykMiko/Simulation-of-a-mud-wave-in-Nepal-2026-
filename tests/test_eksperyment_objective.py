import csv
import math
import os
from pathlib import Path
import subprocess

import pytest

from eksperyment import objective as O, design


def summary(arrival=350, speed=22, **changes):
    return dict(id=O.load()['points'][0]['id'],arrival_s=arrival,speed_m_s=speed,
                speed_samples=13,window_complete=True,**changes)


def test_feasibility_then_iou():
    spec=O.load()
    bad=O.score(.9,[summary(5000)],spec)
    good=O.score(.2,[summary()],spec)
    fast=O.score(.3,[summary(200,40)],spec)
    assert bad['objective_score'] < good['objective_score'] < fast['objective_score']
    assert not bad['dynamics_pass'] and fast['dynamics_pass']
    assert O.score(.3,[summary()],spec)['objective_score']==fast['objective_score']
    for arrival,speed in [(350/5,22/5),(1750,110)]:
        assert O.score(.3,[summary(arrival,speed)],spec)['dynamics_pass']
    for arrival,speed in [(1751,22),(350,111),(None,None)]:
        assert not O.score(.9,[summary(arrival,speed)],spec)['dynamics_pass']


def test_missing_points_and_invalid_scores():
    spec=O.load()
    for iou in (math.nan,math.inf,-.1,1.1):
        with pytest.raises(ValueError): O.score(iou,[summary()],spec)
    with pytest.raises(ValueError): O.score(.5,[],spec)
    with pytest.raises(ValueError): O.ranking_key(dict(number=1,iou_all=.9))
    spec['points'][0]['arrival_tolerance_factor']=1
    with pytest.raises(ValueError): O.validate(spec)


def test_first_wave_only_and_dry_series(tmp_path):
    spec=O.load()
    path=tmp_path/'dynamics.csv'
    with path.open('w') as f:
        writer=csv.writer(f)
        for t in range(5,1001,5):
            wet=3 if t>=350 else 0
            speed=22 if t<=410 else 1000  # Late spike cannot rescue or spoil first-wave score.
            writer.writerow([t,0,wet,wet*2,wet*2*speed])
    result=O.evaluate(path,.4,spec)
    assert result['dynamics_pass']
    assert result['dynamics'][0]['speed_m_s']==22
    assert result['dynamics'][0]['arrival_s']==350
    path.write_text('5,0,0,0,0\n10,0,0,0,0\n')
    assert not O.evaluate(path,.9,spec)['dynamics_pass']
    path.write_text('350,0,3,6,132\n')
    assert not O.evaluate(path,.9,spec)['dynamics_pass']
    path.write_text('5,0,3,6,nan\n')
    with pytest.raises(ValueError): O.evaluate(path,.9,spec)
    path.write_text('5,0,3,6,132\n5,0,3,6,132\n')
    with pytest.raises(ValueError): O.evaluate(path,.9,spec)


def test_optuna_uses_same_score(monkeypatch):
    import optuna
    original=optuna.create_study
    studies=[]
    def capture(*args,**kwargs):
        s=original(*args,**kwargs);studies.append(s);return s
    monkeypatch.setattr(optuna,'create_study',capture)
    history=[]
    for n,p in enumerate(design.batch(1,[])+design.batch(6,[]),1):
        history.append(dict(number=n,params=p,iou_all=.8,
            **O.score(.8,[summary(5000 if n==1 else 350)],O.load())))
    assert len(design.batch(11,history))==5
    assert [t.value for t in studies[0].trials[:10]]==[r['objective_score'] for r in history]
    assert studies[0].trials[0].user_attrs['dynamics_constraints'][0]>0
    assert all(v<=0 for v in studies[0].trials[1].user_attrs['dynamics_constraints'])


def test_no_acceptable_trial_retains_diagnostic_winner(tmp_path,monkeypatch,capsys):
    from eksperyment import study as S
    from ravaflow_tools.checkpoint import atomic,read
    atomic(tmp_path/'protocol.json',{'objective':O.load()})
    monkeypatch.setattr(S,'render',lambda *args:None)
    for n,iou in [(1,.2),(2,.3)]:
        folder=tmp_path/'artifacts'/f't{n:02}'
        atomic(folder/'report.json',dict(number=n,iou_all=iou,metrics={'all':{'jaccard':iou}},
            **O.score(iou,[summary(5000)],O.load())))
        S.accept(tmp_path,n)
    S.finish(tmp_path,{'data_dir':str(tmp_path),'dem':'unused'})
    assert not (tmp_path/'artifacts/t01').exists()
    assert read(tmp_path/'finished.json')['model_status']=='no_acceptable_dynamics'
    assert read(tmp_path/'finished.json')['winner']==2
    assert 'BRAK ZAAKCEPTOWANEJ DYNAMIKI' in capsys.readouterr().out


def test_custom_vm_name_is_used(tmp_path):
    executable=tmp_path/'multipass'
    executable.write_text('#!/bin/sh\nprintf "%s\\n" "$*"\n')
    executable.chmod(0o755)
    script=Path(__file__).resolve().parents[1]/'scripts/eksperyment-vm.sh'
    env={**os.environ,'PATH':str(tmp_path)+':'+os.environ['PATH'],'NEPAL_VM_NAME':'nepal30-dynamics'}
    result=subprocess.run(['bash',str(script),'status'],env=env,text=True,capture_output=True,check=True)
    assert 'exec nepal30-dynamics --' in result.stdout
    env['NEPAL_VM_NAME']='bad/name'
    assert subprocess.run(['bash',str(script),'status'],env=env,capture_output=True).returncode!=0


def test_native_frequent_sampler(tmp_path):
    src=Path(__file__).resolve().parents[1]/'eksperyment'
    spec=tmp_path/'spec'
    spec.write_text('1 5 1\n3 0 1 2\n')
    program=tmp_path/'test.c'
    program.write_text('''
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
struct { int M,N,IMAX,MODEL; } sico={1,3,3,7};
int main(void) {
 int cdomain2[3]={0},px[3]={0},py[3]={0,1,2},ccontinue=1;
 float awt[3][33]={{0}};
 double tmax=30;
 for(int i=0;i<3;i++) {awt[i][0]=1;awt[i][6]=1;awt[i][16]=2;awt[i][13]=20;awt[i][15]=24;}
 for(double tsum=1;tsum<=30;tsum++) {
#include "dynamics.c.inc"
 }
 return 0;
}
''')
    binary=tmp_path/'test'
    subprocess.run(['cc','-I',str(src),str(program),'-lm','-o',str(binary)],check=True)
    out=tmp_path/'dynamics.csv'
    subprocess.run([binary],check=True,env={**os.environ,'NEPAL_DYNAMICS_SPEC':str(spec),'NEPAL_DYNAMICS_CSV':str(out)})
    rows=list(csv.reader(out.open()))
    assert [float(r[0]) for r in rows]==[1,5,10,15,20,25,30]
    assert all(float(r[4])/float(r[3])==22 for r in rows)
