"""
This example demonstrates how to use WATTS to perform
OpenMC calculation.
"""

import time
time_start = time.time()

import numpy as np
import watts
import openmc
import openmc.deplete
from openmc_template import build_openmc_model

import warnings
warnings.filterwarnings("ignore")



params = watts.Parameters()

# # Core design params

params['common_temperature'] = 600
params['enrichment'] = 0.1975
params["H_Zr_ratio"] = 1.6
params['U_met_wo']  = 0.3


# fuel pin dims (these variables' names need to be reviewed!)
params['fuel_radius'] = 0.28575
params['first_gap_radius'] =  0.3175
params['fuel_meat_radius'] =  1.5113
params['second_gap_radius'] = 1.5367
params['cladding_radius'] = 1.5875
params['drum_Absorber_thickness'] = 1
params["pin_gap_distance"] =  0.1
params['assembly_rings'] = 12
params['drum_radius_to_lattice_radius'] = 0.22784810068
params['drum_height_to_lattice_height'] = 1.24

params['power_MW_th'] = 20
params['thermal_efficiency'] = 0.31
params['extra_reflector'] = 14

# params.show_summary(show_metadata=False, sort_by='time')

# # Create OpenMC plugin
openmc_plugin = watts.PluginOpenMC(build_openmc_model, show_stderr=True) # show only error

# Run OpenMC plugin, instructing it to plot the geometry and run a simulation
def run_func():
    #openmc.plot_geometry()
    openmc.run()
    lattice_geometry = openmc.Geometry.from_xml()
    settings = openmc.Settings.from_xml()
    openmc.config['cross_sections'] = "/home/hannbn/projects/MARVEL_MRP/Github_repos/openmc_data/endfb-viii.0-hdf5/cross_sections.xml"
    operator = openmc.deplete.CoupledOperator( openmc.Model(geometry=lattice_geometry, \
        settings=settings), chain_file='/home/hannbn/projects/MARVEL_MRP/Github_repos/openmc_data/simplified_thermal_chain11.xml')
    
    time_steps = [30] * 6
    burnup_step= np.array([0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 60.0, 80.0, 100.0, 120.0, 140.0]) #MWd/kg

    burnup = np.diff ( burnup_step, prepend =0.0 )

    integrator = openmc.deplete.PredictorIntegrator(operator, burnup, 1000000 * params['power_MW_th'], timestep_units='MWd/kg')
    depletion = integrator.integrate(output= 'False')
    
openmc_result = openmc_plugin(params, function=run_func)


# print("KEFF = ", openmc_result.keff)
# print(openmc_result.inputs)
# print(openmc_result.outputs)
# print(openmc_result.tallies[0].get_pandas_dataframe())

# # Open folder in order to view plot images that were produced
# openmc_result.open_folder()

# params.show_summary(show_metadata=True, sort_by='time')

elapsed_time = (time.time() - time_start)/60
print('Execution time:', np.round(elapsed_time,0), 'minutes')

# print(vars(openmc_result))