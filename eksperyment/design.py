"""Trial numbering starts at one. Persist each batch before launching workers."""
import yaml
import optuna
from optuna.distributions import FloatDistribution as FD, CategoricalDistribution as CD
from .build import ROOT

BASE = yaml.safe_load((ROOT/'configs/nepal_v2_final30.yaml').read_text())['fixed_params']
CHANGES = [
    {}, {'mode':'instant'}, {'solid_fraction':.5}, {'release_duration_s':600.0},
    {'delta':6.0}, {'manning':.022}, {'entrainment_log10':-8.0},
    {'entrainment_log10':-9.0}, {'mode':'instant','solid_fraction':.5,'manning':.025},
    {'release_duration_s':120.,'delta':7.,'phi':33.,'entrainment_log10':-8.7},
]
SPACE = dict(phi=FD(30,38), delta=FD(5,12), solid_fraction=FD(.5,BASE['solid_fraction']),
             manning=FD(.02,.045,log=True), entrainment_log10=FD(-9,-7.5),
             mode=CD(('instant','gradual')), release_duration_s=FD(60,600,log=True))


def batch(start, history):
    if start <= 10:
        plans = [{**BASE,'mode':'gradual',**c} for c in CHANGES[start-1:start+4]]
        for p in plans:
            if p['mode']=='instant':
                del p['release_duration_s']
        return plans
    sampler=optuna.samplers.TPESampler(seed=4200+start,n_startup_trials=0,constant_liar=True)
    study=optuna.create_study(direction='maximize',sampler=sampler)
    for report in history:
        p=report['params']
        study.add_trial(optuna.trial.create_trial(params=p,distributions={k:SPACE[k] for k in p},
                                                  value=report['iou_all']))
    plans=[]
    for _ in range(5):
        t=study.ask()
        for key,d in SPACE.items():
            if key=='release_duration_s' and t.params['mode']=='instant':
                continue
            if isinstance(d,CD):
                t.suggest_categorical(key,d.choices)
            else:
                t.suggest_float(key,d.low,d.high,log=d.log)
        plans.append(t.params)
    return plans
