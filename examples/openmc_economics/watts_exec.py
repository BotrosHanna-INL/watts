"""
This example demonstrates how to use WATTS to perform
OpenMC calculation.
"""

# Importing libraries and modules
import time
import numpy as np
import watts
import openmc
import openmc.deplete
from core_design.openmc_template import build_openmc_model
from core_design.utils import calculate_lattice_radius, calculate_reflector_mass, calculate_heat_flux, openmc_depletion, calculate_drum_volume
from reactor_engineering_evaluation.operation import refueling_operation
from core_design.pins_arrangement import rings_1
import warnings

warnings.filterwarnings("ignore")

# Record the current time (in seconds)
time_start = time.time()

# Define user parameters
params = watts.Parameters()

params['plotting'] = "No"

params['cross_sections_xml_location'] =\
    '/home/hannbn/projects/MARVEL_MRP/Github_repos/openmc_data/endfb-viii.0-hdf5/cross_sections.xml'
params['simplified_chain_thermal_xml'] =\
    '/home/hannbn/projects/MARVEL_MRP/Github_repos/openmc_data/simplified_thermal_chain11.xml'

params['power_MW_th'] = 20
params['thermal_efficiency'] = 0.31

params['common_temperature'] = 600 # Kelvins
params['enrichment'] = 0.1975
params["H_Zr_ratio"] = 1.6
params['U_met_wo'] = 0.3

params['fuel_pin_radii'] = [0.28575, 0.3175, 1.5113, 1.5367, 1.5875]
params['fuel_pin_materials'] = ['Zr', None, 'TRIGA_fuel', None, 'SS304']
params['fuel'] = 'TRIGA_fuel'
params['moderator_pin_radii'] = [params['fuel_pin_radii'][3], params['fuel_pin_radii'][4]]
params['moderator_pin_materials'] = ['ZrH', 'SS304']

params["pin_gap_distance"] = 0.1
params['assembly_rings'] = 12
params['lattice_radius'] = \
    calculate_lattice_radius(params['fuel_pin_radii'][-1], params["pin_gap_distance"], params['assembly_rings'])
params['lattice_height'] = 2 * params['lattice_radius']
params['extra_reflector'] = 14
params['hex_area'] = 2.598 * params['lattice_radius'] * params['lattice_radius']
params['core_radius'] = params['lattice_radius'] + params['extra_reflector']
params['reflector'] = 'BeO'
params['coolant'] = 'NaK'

params["deviation angle between drums"] = (np.pi / 14)
params['drum_Absorber_thickness'] = 1
params['control_drum_absorber'] = 'B4C_nat'
params['control_drum_reflector'] = 'Be'
params['drum_radius_to_lattice_radius'] = 0.22784810068
params['Drum_Radius'] = \
    params['drum_radius_to_lattice_radius'] * params['lattice_radius']
params['drum_height_to_lattice_height'] = 1.24
params['drum_height'] = params['drum_height_to_lattice_height'] * params['lattice_height']
params['angle_between_drums_pairs'] = 60
params['drum_gap_distance'] = params['Drum_Radius'] / 90
params['drum_tube_radius'] = params['Drum_Radius'] + params['drum_gap_distance']
params['distance between control drums'] = \
    0.86602540378 * params['lattice_radius'] + params['drum_tube_radius']

params['all_drums_volume'], params['drum_absorp_all_mass'], params['drum_refl_all_mass'] =\
    calculate_drum_volume(params['Drum_Radius'], params['drum_height'],\
        params['drum_Absorber_thickness'], params['angle_between_drums_pairs'])
params['tot_drum_area_all'] = params['all_drums_volume'] / params['drum_height']

params['burnup_steps_MWd_per_Kg'] = [0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 15.0, 20.0,
                                     30.0, 40.0, 50.0, 60.0, 80.0, 100.0, 120.0, 140.0]

params['reflector_mass'] = calculate_reflector_mass(params['hex_area'],\
    params['core_radius'], params['tot_drum_area_all'], params['drum_height'])

params['rings'] = rings_1

params['heat_flux'] = calculate_heat_flux(params['fuel_pin_radii'][-1],\
    params['lattice_height'], params['rings'], params['power_MW_th'])
params['heat_flux_criteria'] = 0.9

params['fuel_pin_count'] = sum(row.count("FUEL") for row in params['rings'])

# operation
params['num_people_required_per_refueling'] = 5
params['levelization_period_years'] = 60 # in years
params['refueling_period_days'] = 15

# Cost parameters
params['FTE_Cost'] = 170000 # $ per FTE

# Function to run OpenMC plugin
def run_openmc(params):
    openmc_plugin = watts.PluginOpenMC(build_openmc_model, show_stderr=True)
    
    def run_func():
        openmc.run()
        lattice_geometry = openmc.Geometry.from_xml()
        settings = openmc.Settings.from_xml()
        depletion_results = openmc_depletion(params, lattice_geometry, settings)
        
        params['fuel_lifetime_days'] = depletion_results[0]
        params['mass_U235'] = depletion_results[1]
        params['mass_U238'] = depletion_results[2]

    if params['heat_flux'] <= params['heat_flux_criteria']:
        openmc_result = openmc_plugin(params, function=run_func)
        
        elapsed_time = (time.time() - time_start) / 60
        print('Execution time:', np.round(elapsed_time, 1), 'minutes')
    
    else:
        print(f"\033[91mHIGH HEAT FLUX: {params['heat_flux']} MW/m^2.\033[0m")


def design_evaluations(params):
    run_openmc(params)
    params['people_by_days_refueling'] = refueling_operation(params)
    params.show_summary(show_metadata=True, sort_by='time')


# Main execution flow
if __name__ == "__main__":
    design_evaluations(params)
    


