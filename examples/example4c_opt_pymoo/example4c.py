# SPDX-FileCopyrightText: 2022 UChicago Argonne, LLC
# SPDX-License-Identifier: MIT

from example4 import *
import numpy as np
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.optimize import minimize
from pymoo.core.problem import Problem

# Problem definition
class fitness_calc(Problem):
    """simple definition of the optimization problem with initialization and evaluation methods"""
    def __init__(self):
        super().__init__(n_var=2, n_obj=2, xl=np.array([0.5, 0.5]), xu=np.array([1.0, 1.0])) # initialization of the problem

    def _evaluate(self, x, out, *args, **kwargs):
        print(x) # this is the matrix of the size of the population
    #    out["F"] = [] 
        for s in x:
            (keff, max_Tf, avg_Tf) = calc_workflow(s) # the workflow is called here and we are saving the results in out['F']
            print(keff, max_Tf, avg_Tf) 
            out["F"].append([keff, max_Tf/avg_Tf]) 

algorithm = NSGA2(pop_size=2) # multicriteria algorithm applied

res = minimize(fitness_calc(),
               algorithm,
               ('n_gen', 2),
               seed=1,
               verbose=False) # this runs the optimization algorithm
print(res.x)
print(res.F)
