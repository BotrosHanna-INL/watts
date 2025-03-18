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


def create_fuel_pin_regions(params):
    # Read what the used decided for the dimensions of the fuel pin
    ## Fuel (these variables' names need to be reviewed!)
    fuel_radii = {'insert': params['fuel_pin_radii'][0],
                'gap1': params['fuel_pin_radii'][1],
                'fuel_meat': params['fuel_pin_radii'][2],
                'gap2': params['fuel_pin_radii'][3],
                'cladding': params['fuel_pin_radii'][4]
    }

    # # Creating surfaces
    shells = [openmc.ZCylinder(r=r) for r in fuel_radii.values()]
    
    region = {'insert': -shells[0],
            'gap1': +shells[0] & -shells[1],
            'fuel_meat': +shells[1] & -shells[2],
            'gap2': +shells[2] & -shells[3],
            'cladding': +shells[3] & -shells[4],
            'coolant': +shells[4]
    }

    return region


def create_moderator_pin_regions(params):
        ## Reflector
    moderator_radii = {'moderator': params['moderator_pin_radii'][0],
                    'cladding': params['moderator_pin_radii'][1]
                    }
    shells = [openmc.ZCylinder(r=r) for r in moderator_radii.values()]
    
    region = {'moderator': -shells[0],
            'cladding': +shells[0] & -shells[1],
            'coolant': +shells[1]
    }
    return region



def create_drums_universe(absorber_thickness, drum_radius,
                          control_drum_absorber_material,
                          control_drum_reflector_material,
                          angle_between_drums_pairs):

    absorber_arc = np.pi/3
    REFERENCE_ANGLE = 0
    rotation_angle = 180

    cd_inner_shell = openmc.ZCylinder(r= drum_radius - absorber_thickness)
    cd_outer_shell = openmc.ZCylinder(r= drum_radius)

    cutting_plane_1 = openmc.Plane(a=1, b=absorber_arc/2)
    cutting_plane_2 = openmc.Plane(a=1, b=-absorber_arc/2)

    drum_absorber = +cd_inner_shell & -cd_outer_shell & -cutting_plane_1 & -cutting_plane_2
    drum_reflector = -cd_outer_shell & ~drum_absorber
    drum_outside = +cd_outer_shell
    drum_absorber = openmc.Cell(name='drum_absorber', fill= control_drum_absorber_material, region=drum_absorber)
    drum_reflector = openmc.Cell(name='drum_reflector', fill= control_drum_reflector_material, region=drum_reflector)
    drum_exterior = openmc.Cell(name='drum_outside', region=drum_outside)

    drum_reference = openmc.Universe(cells=(drum_reflector, drum_absorber, drum_exterior))
    
    drum_cells = []
    for r in range(0, 360, angle_between_drums_pairs):
        dc = openmc.Cell(name=f'drum{r}', fill=drum_reference)
        dc.rotation = [0, 0, REFERENCE_ANGLE + r + rotation_angle]
        drum_cells.append(dc)

    drums = [openmc.Universe(cells=(dc,)) for dc in drum_cells]    
    return drums



def create_assembly_universe(params, fuel_pin_universe, moderator_pin_universe, pin_pitch, reflector_material, outer_universe):

    assembly = openmc.HexLattice()

    assembly.center = (0., 0.)
    assembly.pitch = (pin_pitch,)
    assembly.outer = outer_universe # coolant universe is the outer universe probably
    rings = params['rings']
 
    for i in range(len(rings)):
        for j in range(len(rings[i])):
            if rings[i][j] == 'FUEL':
                rings[i][j] = fuel_pin_universe
            elif rings[i][j] == 'MODERATOR':
                rings[i][j] = moderator_pin_universe
    
    rings = rings[-params['assembly_rings']:]
    assembly.universes = rings
    
    assembly_boundary = openmc.model.hexagonal_prism(edge_length=\
        pin_pitch*(params['assembly_rings']-1)+pin_pitch*0.6, corner_radius = (params['fuel_pin_radii'])[-1] \
            + params["pin_gap_distance"])

    fuel_assembly_cell = openmc.Cell(fill=assembly, region=assembly_boundary)
    reflector_cell = openmc.Cell(fill = reflector_material, region=~assembly_boundary)

    assembly_universe = openmc.Universe(cells=[fuel_assembly_cell, reflector_cell])
    return assembly_universe



def create_control_drums_positions( params, number_of_drums):
        
    # Placement of drums happen by tracing a line through the core apothems
    # then 2 drums are place after each apothem by deviating from this line
    # by a deviation angle
    sector = (params['angle_between_drums_pairs']/180) * np.pi
    
    deviation = params["deviation angle between drums"]  
    positions = []
    for s in range(number_of_drums):
        positions.append(s*sector-deviation)
        positions.append(s*sector+deviation)
    return positions 





def create_core_geometry(params, drums, drums_positions, assembly_universe ):
    cd_distance = params['distance between control drums']
    drum_tube_radius = params['drum_tube_radius']
    drum_universes = []
    for d in drums:
        drum_universes.append(d)
        drum_universes.append(d)

    drum_shells = []
    drum_cells = []
    for p, du in zip(drums_positions, drum_universes):
        x, y = np.cos(p)*cd_distance, np.sin(p)*cd_distance
        drum_shell = openmc.ZCylinder(x0=x, y0=y, r=drum_tube_radius)
        drum_shells.append(drum_shell)
        drum_cell = openmc.Cell(fill=du, region=-drum_shell)
        drum_cell.translation = (x, y, 0)  # translates the center of the drum universe to match the cylinder position
        drum_cells.append(drum_cell)
    
    drums_outside = +drum_shells[0]
    for d in drum_shells[1:]:
        drums_outside = drums_outside & +d

    outer_surface = openmc.ZCylinder(r=params['core_radius'] , boundary_type='vacuum')

    core_cell = openmc.Cell(fill= assembly_universe, region=-outer_surface & drums_outside)

    core_geometry = openmc.Geometry([core_cell] + drum_cells)  
    return core_geometry 


def create_universe_plot(pin_universe, pin_plot_width, num_pixels, font_size,\
    title, fig_size, output_file_name):
    
    pin_plot = pin_universe.plot(width = ( pin_plot_width, pin_plot_width),
                                 pixels=(num_pixels, num_pixels))
    pin_plot.set_xlabel('x [cm]', fontsize= font_size)
    pin_plot.set_ylabel('y [xm]', fontsize= font_size)
    pin_plot.set_title(title, fontsize= font_size)

    pin_plot.tick_params(axis='x', labelsize= font_size)
    pin_plot.tick_params(axis='y', labelsize= font_size)
    
    # Retrieve the figure from the Axes object
    fig = pin_plot.figure
    fig.set_size_inches(fig_size, fig_size) 
    fig.tight_layout()
    # Save the figure to a file
    fig.savefig(output_file_name) 
    
def openmc_depletion(params, lattice_geometry, settings):
    
    openmc.config['cross_sections'] = params['cross_sections_xml_location'] 
    
    # depletion operator, performing transport simulations, is created using the geometry and settings xml files
    operator = openmc.deplete.CoupledOperator(openmc.Model(geometry=lattice_geometry, 
            settings=settings),
            chain_file= params['simplified_chain_thermal_xml'])
    burnup_steps_list_MWd_per_Kg = params['burnup_steps_MWd_per_Kg']
    
    #MWd/kg (MW-day of energy deposited per kilogram of initial heavy metal)
    burnup_step = np.array(burnup_steps_list_MWd_per_Kg)     #np.array([0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 60.0, 80.0, 100.0, 120.0, 140.0]) 
    burnup = np.diff (burnup_step, prepend =0.0 )
    
    # Deplete using a first-order predictor algorithm.
    integrator = openmc.deplete.PredictorIntegrator(operator, burnup,
                                                    1000000 * params['power_MW_th'], timestep_units='MWd/kg')
    integrator.integrate()
    results = openmc.deplete.Results("./depletion_results.h5")
    time, k = results.get_keff()

    time /= (24 * 60 * 60)  # convert back to days from second
    for j, ki in enumerate(k):
        if ki[0] < 1.0:
            i = j-1
            break
    fuel_lifetime_days = (time[i])   
    orig_material = results.export_to_materials(0)

    mass_U235 = orig_material[0].get_mass('U235')
    mass_U238 = orig_material[0].get_mass('U238')
    return fuel_lifetime_days, mass_U235, mass_U238     
     