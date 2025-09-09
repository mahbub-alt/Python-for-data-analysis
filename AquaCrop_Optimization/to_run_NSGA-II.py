from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.operators.sampling.rnd import FloatRandomSampling
from pymoo.operators.crossover.sbx import SBX
from pymoo.operators.mutation.pm import PM
from pymoo.termination.default import DefaultMultiObjectiveTermination
from pymoo.core.problem import Problem
from pymoo.optimize import minimize
from pymoo.termination import get_termination
import warnings
warnings.filterwarnings('ignore')
import numpy as np
import pandas as pd

import Aq_Optimization_functons
from Aq_Optimization_functons import *

class AquaCropCalibration(Problem):
    def __init__(self):
        super().__init__(n_var=10, n_obj=2, n_constr=0, 
                         xl=[30, 50, 50, 30, 0.85, 30, 1.0, 0.45, 0.5, 1.20], 
                         xu=[60, 65, 65, 60, 0.98, 35, 1.10, 0.55, 10, 2.5])

    
    def _evaluate(self, X, out, *args, **kwargs):
        results = []

        for params in X:
            smt_1, smt_2, smt_3, smt_4, CCx, WP, Kcb, HI0, a_HI, Zmax = params
            smt = [smt_1, smt_2, smt_3, smt_4]

            # Run AquaCrop model (assuming `for_objf()` is correctly defined)
            _, all_yld_df, simul_reported = for_objf(smt=smt, CCx=CCx, WP=WP, Kcb=Kcb, HI0=HI0, a_HI=a_HI, Zmax=Zmax, train = False)

            # Compute R² values
            R2_yield = calculate_r2(simul_reported, 'Yield (tonne/ha)', 'Reported_Yield')
            R2_irrigation = calculate_r2(simul_reported, 'Seasonal irrigation (mm)', 'Reported_Irrigation')

            # Transform R² for minimization: (1 - R²)
            obj_yield = 1 - R2_yield
            obj_irrigation = 1 - R2_irrigation

            results.append([obj_yield, obj_irrigation])

        out["F"] = np.array(results)  # Ensure correct shape



def calculate_r2(df, sim_col, obs_col):
    df = df.dropna(subset=[sim_col, obs_col])
    obs = df[obs_col].values
    sim = df[sim_col].values
    ss_total = np.sum((obs - np.mean(obs)) ** 2)
    ss_residual = np.sum((obs - sim) ** 2)
    return 1 - (ss_residual / ss_total) if ss_total != 0 else np.nan  # Avoid division by zero



problem = AquaCropCalibration()

# Configure the algorithm
algorithm = NSGA2(
    pop_size=5,
    sampling=FloatRandomSampling(),
    crossover=SBX(prob=0.9, eta=15),
    mutation=PM(eta=20),
    eliminate_duplicates=False
)

termination = get_termination("n_gen", 1)

# Run the optimization
res = minimize(problem, algorithm, termination, seed=1, verbose=True)


# Find the best solution (highest R² for Yield)
best_index = np.argmin(res.F[:, 0])  # Since we minimize (1 - R²), lower is better
best_params = res.X[best_index]


print("\nOptimized Parameters (Best Solution):", best_params)

F = res.F  # Shape (n_solutions, 2), columns: [1-R2_Yield, 1-R2_Irrigation]

# Convert back to original R² values
R2_Yield = 1 - F[:, 0]
R2_Irrigation = 1 - F[:, 1]

results_df = pd.DataFrame(
    {'R2 (Yield)': R2_Yield, 'R2 (Irrigation)': R2_Irrigation}
)


# Add SMTs and other parameters
results_df['smt_1'] = res.X[:, 0]
results_df['smt_2'] = res.X[:, 1]
results_df['smt_3'] = res.X[:, 2]
results_df['smt_4'] = res.X[:, 3]
results_df['CCx'] = res.X[:, 4]
results_df['WP'] = res.X[:, 5]
results_df['Kcb'] = res.X[:, 6]
results_df['HI0'] = res.X[:, 7]
results_df['a_HI'] = res.X[:, 8]
results_df['Zmax'] = res.X[:, 9]

# Print the DataFrame with the results
print("Optimization Results:",results_df)
results_df.to_csv('output_data/NSGA_II_onlyNW.csv', index=False)