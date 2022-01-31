# SPDX-FileCopyrightText: 2022 UChicago Argonne, LLC
# SPDX-License-Identifier: MIT

from statistics import stdev
from scipy.stats import qmc
from example4 import *

sampler = qmc.LatinHypercube(d=2)
sample = sampler.random(n=10)
l_bounds = [0.5, 0.5]
u_bounds = [1.0, 1.0]
sample_scaled = qmc.scale(sample, l_bounds, u_bounds)

print (sample_scaled)
results = {}
results["keff"] = []
results["max_Tf"] = []
results["avg_Tf"] = []

for X in sample_scaled:
    res = calc_workflow(X)
    results["keff"].append(res[0])
    results["max_Tf"].append(res[1])
    results["avg_Tf"].append(res[2])

print (results)
print (f"keff: mean = { mean(results["keff"]) } , stdev = {stdev(results["keff"])}")
print (f"max_Tf: mean = { mean(results["max_Tf"]) } , stdev = {stdev(results["max_Tf"])}")
print (f"avg_Tf: mean = { mean(results["avg_Tf"]) } , stdev = {stdev(results["avg_Tf"])}")
