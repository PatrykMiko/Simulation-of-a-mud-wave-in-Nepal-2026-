#!/usr/bin/env python3
"""Lightweight perf gate (±2%): compare C++ vs NumPy on 120x120."""
import time, sys, numpy as np
from sciddica_ss3.automaton import SciddicaSS3, ModelParameters
def run(use_cpp, steps=200):
    rng=np.random.default_rng(0)
    alt=(rng.random((120,120))*50).astype(np.float32)
    th=np.zeros((120,120),np.float32); th[55:65,55:65]=5.0
    soil=np.zeros((120,120),np.float32); p=ModelParameters(pc=30.0)
    m=SciddicaSS3(alt.copy(),th.copy(),soil.copy(),params=p,use_cpp=use_cpp)
    for _ in range(5): m.step()
    m=SciddicaSS3(alt.copy(),th.copy(),soil.copy(),params=p,use_cpp=use_cpp)
    t0=time.perf_counter(); [m.step() for _ in range(steps)]; return time.perf_counter()-t0
py=run(False); cpp=run(True); ratio=py/cpp
print(f"py={py:.4f}s cpp={cpp:.4f}s speedup={ratio:.2f}x")
# ±2% gate vs previous baseline (pass baseline as arg)
if len(sys.argv)>1:
    base=float(sys.argv[1]); delta=(cpp-base)/base
    print(f"delta_vs_base={delta*100:+.2f}%"); sys.exit(0 if abs(delta)<=0.02 else 1)
