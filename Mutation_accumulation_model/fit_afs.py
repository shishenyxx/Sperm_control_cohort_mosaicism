## This script describes the model fitting for the stepwise exponential regression for mutation accumulation.

import os
import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
#fit a function: y_t = AF_0/(E^t) for each stage t, t = 1,2,...
#stage width x_t= floor(1/y_t)
#y_0 = AF_0, x_0 = 1
def simulation_to_afs(afs, lower_CI, upper_CI):
    #e ranges from 1.001 to 2, end included
    e_choices = np.linspace(1.001,2,1000)
    #af0 ranges from binomial lower_CI to upper_CI of the MAF, step size = 0.001
    af0_choices = np.arange(lower_CI,upper_CI + 0.001, 0.001)
    best_residual = float("inf")
    best_af0 = 0
    best_e = 0
    best_simulations = []
    #loop through all combinations
    for af0 in af0_choices:
        for e in e_choices:
            i = 0
            simulations = [af0]
            while len(simulations) < len(afs):
                i += 1
                yn = af0/(e**i)
                xn = math.floor(1/yn)
                simulations += [yn] * xn
            simulations = np.array(simulations[:len(afs)])
            diff = afs - simulations
            residual = sum(abs(diff)) #residual uses L1 loss
            #record the best result
            if residual <= best_residual:
                best_residual = residual
                best_af0 = af0
                best_e = e
                best_simulations = simulations
    #compute r_square from the pearson correlation coefficient
    pearson_r = pearsonr(afs, best_simulations)
    r_square = pearson_2 ** 2
    return best_simulations, best_af0, best_e, best_residual, r_square
#construct input dict
afs_dict = {}
for f in os.listdir(dirpath):
    data = pd.read_csv(dirpath + f, header = None)
    afs = data[0].values
    lower_CI, upper_CI = list(map(float, f.split("_")[-1][:-4].split("-")))
    afs_dict[f] = [afs, lower_CI, upper_CI]
    
#for each file fit the simulated line
afs_results = {}
for file in afs_dict.keys():
    afs, lower_CI, upper_CI = afs_dict[file]
    simulations, af0, e, residual, r_square = simulation_to_afs(afs, lower_CI, upper_CI)
    afs_results[file] = [simulations, af0, e, residual, r_square]
