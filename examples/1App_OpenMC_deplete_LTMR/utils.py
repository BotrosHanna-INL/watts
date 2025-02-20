import numpy as np
import openmc


def create_cells(regions:dict, materials:list)->dict:
    return {key:openmc.Cell(name=key, fill=mat, region=value) for (key,value), mat in zip(regions.items(), materials)}


def circle_area(r):
    return (np.pi) * r **2


def circle_perimeter(r):
    return 2*(np.pi)*r


def cylinder_radial_shell(r, h):
    return circle_perimeter(r) * h


def calculate_lattice_radius(cladding_radius, pin_gap_distance, number_of_assembly_rings):
    lattice_radius = (cladding_radius * 2 + pin_gap_distance) * number_of_assembly_rings
    return lattice_radius

def calculate_heat_flux(cladding_radius, lattice_height, rings, power_MW_th):
    fuel_number =  sum(r.count("FUEL") for r in rings)
    heat_transfer_surface = cylinder_radial_shell(cladding_radius, lattice_height ) * fuel_number  * 1e-4 # convert from cm2 to m2
    
    return power_MW_th/heat_transfer_surface # MW/m^2


def calculate_drum_volume(DRUM_RADIUS, drum_height, absorber_thickness, angle_between_drums_pairs):
    drum_volume = 3.14*DRUM_RADIUS * DRUM_RADIUS *drum_height
    
    drum_absorp_vol = (3.14*( DRUM_RADIUS * DRUM_RADIUS - (DRUM_RADIUS-absorber_thickness)*(DRUM_RADIUS-absorber_thickness) )*drum_height)/3
    drum_refl_vol = drum_volume - drum_absorp_vol 
    
    number_of_drums = 2* (360/angle_between_drums_pairs)
    all_drums_volume = drum_volume * number_of_drums
    
    drum_absorp_vol_all = drum_absorp_vol  * number_of_drums  
    drum_refl_vol_all = drum_refl_vol  * number_of_drums  
    
    drum_absorp_all_mass = drum_absorp_vol_all * 2.52/1000 # B4C (in Kg)
    drum_refl_all_mass = drum_refl_vol_all  * 3.02/1000 # BeO (in Kg)
    
    return all_drums_volume, drum_absorp_all_mass, drum_refl_all_mass


def calculate_reflector_mass(hex_area, core_radius, area_of_all_drums, drum_height):
    # I assume for now that the drums are always fully inside the reflector
    
    area_reflector = 3.14 * core_radius * core_radius - hex_area  - area_of_all_drums # cm2
    vol_reflector = area_reflector * drum_height # cm^3
    mass_reflector = vol_reflector * 3.02/1000 # mass in Kg
    return mass_reflector