# Importing libraries
import numpy as np
import openmc
import openmc.model
from utils import create_cells, circle_area, create_fuel_pin_regions, create_moderator_pin_regions
from openmc_materials_database import collect_materials_data, color_materials

"""
An OpenMC function that accepts an instance of "parameters" 
and generates the necessary XMl files
"""

def build_openmc_model(params):
    
    materials_database = collect_materials_data(params)
    
    fuel = materials_database [params['fuel']]
    fuel_pin_filler_rod = materials_database[params['fuel_pin_filler_rod']]
    cladding = materials_database[params['cladding']]
    coolant = materials_database[params['coolant']]
    moderator = materials_database[params['moderator_pin']] 
    reflector = materials_database[params['reflector']]
    control_drum_absorber = materials_database[params['control_drum_absorber']]
    control_drum_reflector = materials_database[params['control_drum_reflector']]
    

    
    
    # **************************************************************************************************************************
    #                                                Sec. 2 : GEOMETRY: FUEL, MODERATOR, ASSEMBLY
    # **************************************************************************************************************************
    
    region = create_fuel_pin_regions(params)
    
    # ## Fuel (these variables' names need to be reviewed!)
    # fuel_radii = {'insert': params['fuel_radius'],
    #             'gap1': params['first_gap_radius'],
    #             'fuel_meat': params['fuel_meat_radius'],
    #             'gap2': params['second_gap_radius'],
    #             'cladding': params['cladding_radius']
    # }

    # shells = [openmc.ZCylinder(r=r) for r in fuel_radii.values()]
    
    # region = {'insert': -shells[0],
    #         'gap1': +shells[0] & -shells[1],
    #         'fuel_meat': +shells[1] & -shells[2],
    #         'gap2': +shells[2] & -shells[3],
    #         'cladding': +shells[3] & -shells[4],
    #         'coolant': +shells[4]
    # }

    fuel_materials = [fuel_pin_filler_rod, None,\
        fuel, None, cladding , coolant]

    fuel_cells = create_cells(region, fuel_materials)

    fuel_pin = openmc.Universe(cells=fuel_cells.values())

    
    
    
    ## Reflector
    moderator_pin_region = create_moderator_pin_regions(params)

    # moderator_radii = {'moderator': fuel_radii['gap2'],
    #                 'cladding': fuel_radii['cladding']}

    # shells = [openmc.ZCylinder(r=r) for r in moderator_radii.values()]

    # region = {'moderator': -shells[0],
    #         'cladding': +shells[0] & -shells[1],
    #         'coolant': +shells[1]
    # }

    moderator_materials = [moderator, cladding, coolant]

    moderator_cells = create_cells(moderator_pin_region, moderator_materials)

    moderator_pin = openmc.Universe(cells=moderator_cells.values())

    reflector_cell = openmc.Cell(fill=reflector)
    reflector_universe = openmc.Universe(cells=(reflector_cell,))

    coolant_cell = openmc.Cell(fill=coolant)
    coolant_universe = openmc.Universe(cells=(coolant_cell,))


    # Assembly
    pin_pitch = (params['fuel_pin_radii'][-1] )* 2 + params["pin_gap_distance"]
    
    
    # **************************************************************************************************************************
    #                                                Sec. 3 : CONTROL DRUMS
    # **************************************************************************************************************************

    
    

    ABSORBER_THICKNESS = params['drum_Absorber_thickness']

    absorber_arc = np.pi/3
    REFERENCE_ANGLE = 0

    rotation_angle = 180

    cd_inner_shell = openmc.ZCylinder(r= params['Drum_Radius'] - ABSORBER_THICKNESS)
    cd_outer_shell = openmc.ZCylinder(r= params['Drum_Radius'])

    cutting_plane_1 = openmc.Plane(a=1, b=absorber_arc/2)
    cutting_plane_2 = openmc.Plane(a=1, b=-absorber_arc/2)

    drum_absorber = +cd_inner_shell & -cd_outer_shell & -cutting_plane_1 & -cutting_plane_2
    drum_reflector = -cd_outer_shell & ~drum_absorber
    drum_outside = +cd_outer_shell
    drum_absorber = openmc.Cell(name='drum_absorber', fill=control_drum_absorber, region=drum_absorber)
    drum_reflector = openmc.Cell(name='drum_reflector', fill= control_drum_reflector, region=drum_reflector)
    drum_exterior = openmc.Cell(name='drum_outside', region=drum_outside)

    drum_reference = openmc.Universe(cells=(drum_reflector, drum_absorber, drum_exterior))

    drum_cells = []
    for r in range(0, 360, params['angle_between_drums_pairs']):
        dc = openmc.Cell(name=f'drum{r}', fill=drum_reference)
        dc.rotation = [0, 0, REFERENCE_ANGLE + r + rotation_angle]
        drum_cells.append(dc)

    drums = [openmc.Universe(cells=(dc,)) for dc in drum_cells]

    # for d in drums:
    #     d.plot(width=(DRUM_RADIUS*2, DRUM_RADIUS*2), color_by='material', colors=colors)

    # **************************************************************************************************************************
    #                                                Sec. 4 : ASSEMBLY & RINGS
    # **************************************************************************************************************************
    
    
    
    assembly = openmc.HexLattice()

    assembly.center = (0., 0.)
    assembly.pitch = (pin_pitch,)
    assembly.outer = coolant_universe

    rings = params['rings']
 
   
    for i in range(len(rings)):
        for j in range(len(rings[i])):
            if rings[i][j] == 'FUEL':
                rings[i][j] =fuel_pin
            elif rings[i][j] == 'MODERATOR':
                rings[i][j] = moderator_pin
    
    rings = rings[-params['assembly_rings']:]

    
    assembly.universes = rings
   
    # Number of fuel elements and moderator elements
    fuel_number = sum(r.count(fuel_pin) for r in rings)

   
    # **************************************************************************************************************************
    #                                                Sec. 5 : VOLUME INFO for Depletion
    # **************************************************************************************************************************
    
    
    
    # fissile_area = circle_area(fuel_radii['fuel_meat']) - circle_area(fuel_radii['gap1'])

    # fuel.volume = fissile_area *params['lattice_height']  * fuel_number
    

    # heat_transfer_surface = cylinder_radial_shell(fuel_radii['cladding'], params['lattice_height'] ) * fuel_number  * 1e-4 # convert from cm2 to m2

    power_MW_th = params['power_MW_th']
    power_MW_e = power_MW_th * params['thermal_efficiency']

    materials = openmc.Materials([fuel, moderator, coolant, fuel_pin_filler_rod,\
        cladding, control_drum_reflector, reflector, control_drum_absorber])
    # fuel_pin, moderator_pin, coolant, fuel_pin_filler_rod, cladding, control_drum_reflector, reflector , control_drum_absorber])
   
    openmc.Materials.cross_sections = params['cross_sections_xml_location']
    materials.export_to_xml()
    
    
    assembly_boundary = openmc.model.hexagonal_prism(
        edge_length=pin_pitch*(params['assembly_rings']-1)+
        pin_pitch*0.6, corner_radius=params['fuel_pin_radii'][-1] +
        params["pin_gap_distance"])

    fuel_assembly_cell = openmc.Cell(fill=assembly, region=assembly_boundary)
    reflector_cell = openmc.Cell(fill=reflector, region=~assembly_boundary)

    assembly_universe = openmc.Universe(cells=[fuel_assembly_cell, reflector_cell])
    
    
    
    # **************************************************************************************************************************
    #                                                Sec. 6 : CORE DRUM REPLACEMENT
    # **************************************************************************************************************************
    
    
    
    drum_gap_distance =  params['Drum_Radius']/90 # it was 0.1 and I made it as a ratio
    drum_tube_radius = params['Drum_Radius'] + drum_gap_distance

    # Placement of drums happen by tracing a line through the core apothems
    # then 2 drums are place after each apothem by deviating from this line
    # by a deviation angle

    sector = np.pi/3
    deviation = (np.pi/14 )

    # i replaced the formula of Rodrigo with another one
    cd_distance = 0.86602540378 * params['lattice_radius']  + drum_tube_radius 
    positions = []
    for s in range(6):
        positions.append(s*sector-deviation)
        positions.append(s*sector+deviation)

    drum_universes = []
    for d in drums:
        drum_universes.append(d)
        drum_universes.append(d)

    drum_shells = []
    drum_cells = []
    for p, du in zip(positions, drum_universes):
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

    core_cell = openmc.Cell(fill=assembly_universe, region=-outer_surface & drums_outside)

    core_geometry = openmc.Geometry([core_cell] + drum_cells)
    core_geometry.export_to_xml()
    
    
    # **************************************************************************************************************************
    #                                                Sec. 7 : SIMULATION
    # **************************************************************************************************************************
    
    
    
    point = openmc.stats.Point((0, 0, 0))
    source = openmc.Source(space=point)

    settings = openmc.Settings()
    settings.source = source
    settings.batches = 100
    settings.inactive = 50
    settings.particles = 500

    settings.export_to_xml()
    
    
    cell_filter = openmc.CellFilter(fuel_cells['fuel_meat'])

    tally = openmc.Tally(1)
    tally.filters = [cell_filter]

    tally.nuclides = ['U235']
    tally.scores = ['total', 'fission', 'absorption', '(n,gamma)']

    tallies = openmc.Tallies([tally])
    
    
    energies = np.logspace(np.log10(1e-5), np.log10(20.0e6), 501)

    e_filter = openmc.EnergyFilter(energies)

    # Create tally with energy filter
    equal_lethargy_tally = openmc.Tally()
    equal_lethargy_tally.filters = [e_filter]
    equal_lethargy_tally.scores = ['flux']

    # Set model tallies
    tallies.append(equal_lethargy_tally)

    tallies.export_to_xml()