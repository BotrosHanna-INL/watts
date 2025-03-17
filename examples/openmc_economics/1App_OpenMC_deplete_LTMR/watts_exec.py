"""
This example demonstrates how to use WATTS to perform
OpenMC calculation.
"""


#####
# Section 0
# Importing libraries and modules

import time
import numpy as np

# importing WATTS and OpenMC
import watts
import openmc
import openmc.deplete

#Using the OpenMc template
from openmc_template import build_openmc_model

from utils import *
from pins_arrangement import rings_1

import warnings
warnings.filterwarnings("ignore")
# Record the current time (in seconds) 
time_start = time.time()

#####
# Section 1
# The User-Defined Parameters
params = watts.Parameters()

# The location of the OpenMC Data
params['cross_sections_xml_location'] = '/home/hannbn/projects/MARVEL_MRP/Github_repos/openmc_data/endfb-viii.0-hdf5/cross_sections.xml'

# System Variables
params['power_MW_th'] = 20
params['thermal_efficiency'] = 0.31

# # Core design params
params['common_temperature'] = 600  # Temperature of the material in Kelvin
params['enrichment'] = 0.1975 # The Uranium enrichment (wt%). Its range is 0 to 1
params["H_Zr_ratio"] = 1.6 # #The proportion of hydrogen atoms relative to zirconium atoms

params['U_met_wo']  = 0.3  # Uranium Weight/ Fuel Weight. Its range is 0 to 1

# The fuel pin and moderator pin
params['fuel_pin_radii'] = [0.28575, 0.3175, 1.5113, 1.5367, 1.5875]
params['fuel_pin_materials'] = ['Zr', None,'TRIGA_fuel', None, 'SS304' ]
params['moderator_pin_radii'] = [params['fuel_pin_radii'][3], params['fuel_pin_radii'][4]]
params['moderator_pin_materials'] = ['ZrH' , 'SS304']
params["pin_gap_distance"] =  0.1
params['assembly_rings'] = 12 # number of rings
params['lattice_radius'] = calculate_lattice_radius(params['fuel_pin_radii'][-1],
                                                    params["pin_gap_distance"], params['assembly_rings'] )
params['lattice_height'] = 2 * params['lattice_radius']
params['extra_reflector'] = 14
#Hexagon area : https://en.wikipedia.org/wiki/Hexagon
params['hex_area'] = 2.598* params['lattice_radius']  * params['lattice_radius']
params['core_radius'] =   params['lattice_radius']  + params['extra_reflector']

# Reflector and coolant
params['reflector'] = 'BeO'
params['coolant'] = 'NaK'


# Control drums
params["deviation angle between drums"] = (np.pi/14 )
params['drum_Absorber_thickness'] = 1
params['control_drum_absorber'] = 'B4C_nat'
params['control_drum_reflector'] = 'Be'
params['drum_radius_to_lattice_radius'] = 0.22784810068
params['Drum_Radius'] = params['drum_radius_to_lattice_radius'] * params['lattice_radius'] # default value is 0.22784810068
params['drum_height_to_lattice_height'] = 1.24
params['drum_height'] = params['drum_height_to_lattice_height'] * params['lattice_height']   # since in MARVEL, the control drum height is 1.24* active height
params['angle_between_drums_pairs'] = 60
params['drum_gap_distance'] = params['Drum_Radius']/90
params['drum_tube_radius'] = params['Drum_Radius'] + params['drum_gap_distance'] 
params['distance between control drums'] = 0.86602540378 * params['lattice_radius']  + params['drum_tube_radius']
params['all_drums_volume'], params['drum_absorp_all_mass'] , params['drum_refl_all_mass'] = \
    calculate_drum_volume(params['Drum_Radius'], params['drum_height'], params['drum_Absorber_thickness'] , params['angle_between_drums_pairs'])
    
    
params['tot_drum_area_all'] =  params['all_drums_volume'] /params['drum_height']

    
params['all_drums_volume', 'drum_absorp_all_mass', 'drum_refl_all_mass'] =\
    calculate_drum_volume(params['Drum_Radius'], params['drum_height'], params['drum_Absorber_thickness'] , params['angle_between_drums_pairs'])



params['reflector_mass'] = calculate_reflector_mass(params['hex_area'], params['core_radius'], params['tot_drum_area_all'], params['drum_height']) 

# assembly

params['rings'] = rings_1
params['heat_flux'] = calculate_heat_flux(params['fuel_pin_radii'][-1],
                                          params['lattice_height'],
                                          params['rings'], params['power_MW_th'])
params['heat_flux_criteria'] = 0.9 

# Count occurrences of the placeholder for fuel_pin
fuel_pin_count = sum(row.count("FUEL") for row in params['rings'] ) 

   
# Run OpenMC plugin, instructing it to plot the geometry and run a simulation
def run_func():
    openmc.run()

openmc_plugin = watts.PluginOpenMC(build_openmc_model, show_stderr=True) # show only error
# params.show_summary(show_metadata=True, sort_by='time')


if params['heat_flux'] <= params['heat_flux_criteria']:

    openmc_result = openmc_plugin(params, function=run_func)

    elapsed_time = (time.time() - time_start)/60
    print('Execution time:', np.round(elapsed_time, 1), 'minutes')

elif params['heat_flux'] > params['heat_flux_criteria']:
    print(f"\033[91mHIGH HEAT FLUX: {params['heat_flux']} MW/m^2.\033[0m")
    