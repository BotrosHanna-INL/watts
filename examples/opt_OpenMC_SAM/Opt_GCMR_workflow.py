from math import cos, pi
import watts
from scipy.optimize import minimize
from statistics import mean
from openmc_template import build_openmc_model


watts.tee_stdout = False # This doesn't seem to work.
watts.tee_stderr = False

params = watts.Parameters()
# TH params
params['He_inlet_temp'] = 600 + 273.15  # K
params['He_outlet_temp'] = 850 + 273.15 # K
params['He_cp'] = 5189.2 # J/kg-K
params['He_K'] =  0.32802   # W/m-K
params['He_density'] = 3.8815   # kg/m3
params['He_viscosity'] = 4.16e-5 # Pa.s
params['He_Pressure'] = 7e6    # Pa
params['Tot_assembly_power'] = 250000 # W

for i in range(1, 6):
    params[f'Init_P_{i}'] = 1 # Fraction

# Core design params
params['ax_ref'] = 20 # cm
params['num_cool_pins'] = 1*6+2*6+6*2/2
params['num_fuel_pins'] = 6+6+6+3*6+2*6/2+6/3
params['Height_FC'] = 2.0 # m
params['Lattice_pitch'] = 2.0
params['Assembly_pitch'] = 7.5 * 2 * params['Lattice_pitch'] / (cos(pi/6) * 2)
params['lbp_rad'] = 0.25 # cm
params['mod_ext_rad'] = 0.90 # cm
params['shell_thick'] = 0.05   # FeCrAl
params['liner_thick'] = 0.007  # Cr
params['control_pin_rad'] = 0.99 # cm
# Control use of S(a,b) tables
params['use_sab'] = True
params['use_sab_BeO'] = True
params['use_sab_YH2'] = False
# OpenMC params
params['cl'] = params['Height_FC']*100 - 2 * params['ax_ref'] # cm
params['pf'] = 40 # percent
params['num_cpu'] = 60

# initial X values
X = [0.9, 0.6]

def calc_workflow(X):
    """ example of workflow calculation that includes SAM and OpenMC """

    # params to optimize
    params['FuelPin_rad'] = X[0] # cm
    params['cool_hole_rad'] = X[1] # cm
    params['Coolant_channel_diam'] = (params['cool_hole_rad'] * 2)/100 # in m
    params['Graphite_thickness'] = (params['Lattice_pitch'] - params['FuelPin_rad'] - params['cool_hole_rad']) # cm

    print("FuelPin_rad / cool_hole_rad", X[0], X[1])

    # SAM workflow
    sam_plugin = watts.PluginSAM('../initial_SAM/sam_template')
    sam_plugin.sam_exec = "/home/rhu/projects/SAM/sam-opt"
    sam_result = sam_plugin.workflow(params)
    max_Tf = max(sam_result.csv_data[f'max_Tf_{i}'][-1] for i in range(1, 6))
    avg_Tf = mean(sam_result.csv_data[f'avg_Tf_{i}'][-1] for i in range(1, 6))
    print("MaxTfuel / AvgTfuel= ", max_Tf, avg_Tf)

    # get temperature from SAM results
    params['temp'] = mean([sam_result.csv_data[f'avg_Tgraphite_{i}'][-1] for i in range(1, 6)])
    for i in range(1, 6):
        params[f'temp_F{i}'] = sam_result.csv_data[f'avg_Tf_{i}'][-1]

    # Run OpenMC plugin
    openmc_plugin = watts.PluginOpenMC(build_openmc_model)
    openmc_result = openmc_plugin.workflow(params)
    print("KEFF = ", openmc_result.keff)

    fitness = abs(openmc_result.keff.n - 1) + (max_Tf/avg_Tf)
    return fitness


# optimization function - only 10 maximum iterations to make it run quick!
res = minimize(calc_workflow, X, method ='SLSQP', bounds=((0.5, 1.0), (0.5, 0.99)), options={'maxiter': 10, 'iprint': 1, 'disp': False, 'eps': 0.01})
params.show_summary(show_metadata=True, sort_by='time')


X = res.x
print("optimum X(FuelPin_rad, cool_hole_rad) = ", X)
