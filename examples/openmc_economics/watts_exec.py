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
from core_design.utils import calculate_lattice_radius, calculate_reflector_mass,\
    calculate_heat_flux, openmc_depletion, calculate_drum_volume
from reactor_engineering_evaluation.tools import cylinder_annulus_mass
from reactor_engineering_evaluation.operation import reactor_operation
from reactor_engineering_evaluation.fuel_calcs import fuel_calculations
from reactor_engineering_evaluation.vessels_calcs import vessels_specs
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
params['extra_reflector'] = 14 # reflector thickness
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


#Shielding
params['in_vessel_shield_thickness'] = 10 #cm
params['in_vessel_shielding_inner_radius'] = params['core_radius'] 
params['in_vessel_shielding_outer_radius'] = params['core_radius'] + params['in_vessel_shield_thickness']
params['in_vessel_material'] = 'boron_carbide' 


params['out_of_vessel_shield_thickness'] = 39.37 #cm
params['out_vessel_shield_material'] = 'water_extended_polymer'
# The out of vessel shield is not fully made of the out of vessel material (e.g. WEP) so we use an effective density factor
params['out_vessel_shield_effective_density_factor'] = 0.5

# Vessels parameters
params['vessel_radius'] = params['core_radius'] + params['in_vessel_shield_thickness']   + 20 # cm
params['vessel_thickness'] = 2 # cm
params['vessel_lower_plenum_height'] = 30 # cm
params['vessel_upper_plenum_height'] = 60 # cm
params['vessel_upper_gas_gap'] = 10 # cm
params['vessel_bottom_depth'] = 35  # cm
params['vessel_material'] ='stainless_steel'

params['gap_between_vessel_and_guard_vessel'] = 2 # cm
params['guard_vessel_thickness'] = 0.5
params['guard_vessel_material'] ='stainless_steel'

params['gap_between_guard_vessel_and_cooling_vessel'] = 5 # cm
params['cooling_vessel_thickness'] = 5 # cm
params['cooling_vessel_material'] ='stainless_steel'

params['gap_between_cooling_vessel_and_intake_vessel'] = 0.3 # cm
params['intake_vessel_thickness'] = 0.3 # cm
params['intake_vessel_material'] ='stainless_steel'


# operation
params['num_people_required_per_refueling'] = 5
params['num_people_required_per_startup'] = 4

params['levelization_period_years'] = 60 # in years
params['refueling_period_days'] = 15
params['number_of_unanticipated_shutdowns_per_year']= 0.8

params['duration_to_startup_after_refueling_days'] = 7
params['duration_to_startup_after_shutdown_days'] = 14
params['reactors_monitored_by_one_person'] = 5 
params['FTEs_for_security_staff'] = 5 







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
        
        params['fuel_lifetime_days'] = depletion_results[0] # days
        params['mass_U235'] = depletion_results[1] # grams
        params['mass_U238'] = depletion_results[2] # grams

    if params['heat_flux'] <= params['heat_flux_criteria']:
        openmc_result = openmc_plugin(params, function=run_func)
        
        elapsed_time = (time.time() - time_start) / 60
        print('Execution time:', np.round(elapsed_time, 1), 'minutes')
    
    else:
        print(f"\033[91mHIGH HEAT FLUX: {params['heat_flux']} MW/m^2.\033[0m")


## TEMPORARY  ## DELETE LATER!!!!!!!!!!!!!!!!!!!!
params['fuel_lifetime_days'] = 2078 # days
params['mass_U235'] = 67711.4 # grams
params['mass_U238'] = 278650.8  # grams



def design_evaluations(params):
    # run_openmc(params)
    params['people_by_days_refueling_per_year'], params['people_by_days_startup_per_year'],\
        params['capacity_factor'] = reactor_operation(params)

    
    params['natural_U_mass_consumption_Kg'], params['fuel_tail_waste_mass_Kg'], params['SWU_kg'] =\
        fuel_calculations(params)
        
    params['vessels_total_radius'], params['vessel_height'] , params['vessels_total_height'],\
        params['vessel_mass_kg'], params['guard_vessel_mass_kg'] ,\
            params['cooling_vessel_mass'], params['intake_vessel_mass_kg'] = vessels_specs(params)
            
    
    # in vessel shielding mass (kilograms)
    params['in_vessel_shielding_mass'] = cylinder_annulus_mass(params['in_vessel_shielding_outer_radius'],\
    params['in_vessel_shielding_inner_radius'], params['vessel_height'], params['in_vessel_material'] )  

    params['out_of_vessel_shielding_mass'] = params['out_vessel_shield_effective_density_factor'] * cylinder_annulus_mass(params['out_of_vessel_shield_thickness']+ params['vessels_total_radius'],\
        params['out_of_vessel_shield_thickness'], params['vessels_total_height'], params['out_vessel_shield_material']) 

    params.show_summary(show_metadata=True, sort_by='time')


# Main execution flow
if __name__ == "__main__":
    design_evaluations(params)
    


