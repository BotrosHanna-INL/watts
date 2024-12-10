# SPDX-FileCopyrightText: 2022-2023 UChicago Argonne, LLC
# SPDX-License-Identifier: MIT

import openmc
import openmc.model
import numpy as np
from utils import create_cells, circle_area, cylinder_radial_shell

def build_openmc_model(params):
    """ OpenMC Model """
    
    
    """
    ***************************************************************************************************************************
                                                    Sec. 1 : MATERIALS
    ***************************************************************************************************************************
    """
    
    
     # Materials properties
    common_temperature = params['common_temperature']

    ## TRIGA Fuel

    # The fuel is declared following a mechanism that is close to the logic of the
    # fabrication of the actual TRIGA fuel. This type of fuel is specified for example
    # as "45/20" fuel, which implies that 45% weight is composed of Uranium (metal)
    # that is 20% enriched in a matrix of ZrH. TRIGA fuel can also contain 3% weight of
    # Erbium as a burnable absorber in the fuel meat.

    # First let's declare the individual components of the fuel

    
    U_met = openmc.Material(name="U_met")
    U_met.set_density("g/cm3", 19.05)

    U_met.add_nuclide("U235", params['enrichment'])
    U_met.add_nuclide("U238", 1 - params['enrichment'])

    ZrH_fuel = openmc.Material(name="ZrH_fuel")
    ZrH_fuel.set_density("g/cm3", 5.63)

    ZrH_fuel.add_element("zirconium", 1.0)
    ZrH_fuel.add_nuclide("H1", params["H_Zr_ratio"])

    # Now we mix the components in the right weight %
    # The NRC seems to license up to TRIGA fuel with up to 45% weight Uranium

    TRIGA_fuel = openmc.Material.mix_materials(
        [U_met, ZrH_fuel], [params['U_met_wo'], 1 - params['U_met_wo']], "wo", name="UZrH"
    )

    TRIGA_fuel.temperature = common_temperature
    TRIGA_fuel.add_s_alpha_beta("c_H_in_ZrH")

    # Let's also make a version with 3% Erbium in the meat

    Er = openmc.Material(name="Er", temperature=common_temperature)
    Er.set_density("g/cm3", 9.2)
    Er.add_element("erbium", 1.0)

    BA_wo = 0.03

    TRIGA_fuel_BA = openmc.Material.mix_materials(
        [U_met, Er, ZrH_fuel], [params['U_met_wo'], BA_wo, 1 - params['U_met_wo'] - BA_wo], "wo", name="UZrH-Er"
    )

    TRIGA_fuel_BA.temperature = common_temperature
    TRIGA_fuel_BA.add_s_alpha_beta("c_H_in_ZrH")


    ## Moderator pins
    ZrH = openmc.Material(name="ZrH", temperature=common_temperature)
    ZrH.set_density("g/cm3", 5.6)

    ZrH.add_nuclide("H1", 1.85)
    ZrH.add_element("zirconium", 1.0)
    ZrH.add_s_alpha_beta("c_H_in_ZrH")

    ## Coolant
    NaK = openmc.Material(name="NaK", temperature=common_temperature)
    NaK.set_density("g/cm3", 0.75)
    NaK.add_nuclide("Na23", 2.20000e-01)
    NaK.add_nuclide("K39", 7.27413e-01)
    NaK.add_nuclide("K41", 5.24956e-02)

    ## Reflectors
    Be = openmc.Material(name="Be")
    Be.add_element("beryllium", 1.0)
    Be.add_s_alpha_beta("c_Be")
    Be.set_density("g/cm3", 1.84)
    Be.temperature = common_temperature
    BeO = openmc.Material(name="BeO", temperature=common_temperature)
    BeO.set_density("g/cm3", 3.01)
    BeO.add_element("beryllium", 1.0)
    BeO.add_element("oxygen", 1.0)
    BeO.add_s_alpha_beta("c_Be_in_BeO")

    ## Structural materials
    Zr = openmc.Material(name="Zr", temperature=common_temperature)
    Zr.set_density("g/cm3", 6.49)
    Zr.add_element("zirconium", 1.0)
    SS304 = openmc.Material(name="SS304", temperature=common_temperature)
    SS304.set_density("g/cm3", 7.98)
    SS304.add_element("carbon", 0.04)
    SS304.add_element("silicon", 0.50)
    SS304.add_element("phosphorus", 0.023)
    SS304.add_element("sulfur", 0.015)
    SS304.add_element("chromium", 19.00)
    SS304.add_element("manganese", 1.00)
    SS304.add_element("iron", 70.173)
    SS304.add_element("nickel", 9.25)

    SS316 = openmc.Material(name="SS316", temperature=common_temperature)
    SS316.set_density("g/cm3", 7.98)
    SS316.add_element("carbon", 0.041)
    SS316.add_element("silicon", 0.507)
    SS316.add_element("phosphorus", 0.023)
    SS316.add_element("sulfur", 0.015)
    SS316.add_element("chromium", 17.00)
    SS316.add_element("manganese", 1.014)
    SS316.add_element("iron", 66.90)
    SS316.add_element("nickel", 12.00)
    SS316.add_element("molybdenum", 2.50)
    
    
    ## Absorbers
    B4C_nat = openmc.Material(name="B4C", temperature=common_temperature)
    B4C_nat.add_element("boron", 4)
    B4C_nat.add_element("carbon", 1)
    B4C_nat.set_density("g/cm3", 2.52)
    B4C_rich = openmc.Material(name="B4C_rich", temperature=common_temperature)
    B4C_rich.add_element("boron", 4, enrichment=50.0, enrichment_target="B10")
    B4C_rich.add_element("carbon", 1)
    B4C_rich.set_density("g/cm3", 2.52)

    # Adding Hafnium in the Zr insert of the fuel (so non-purified zirconium)
    Hf = openmc.Material(name="Hf")
    Hf.set_density("g/cm3", 13.31)
    Hf.add_element("hafnium", 1.0)
    Hf_impurity = 0.03  # About 2-3% Hf in Zirconium
    ZrHf_nat = openmc.Material.mix_materials(
        [Zr, Hf], [1 - Hf_impurity, Hf_impurity], "ao", name="ZrHfnat"
    )

    ZrHf_nat.temperature = common_temperature

    Gd = openmc.Material(name="Gd")
    Gd.set_density("g/cm3", 7.9)
    Gd.add_element("gadolinium", 1.0)

    # Adding Gadolinium to the Zr insert of the fuel
    Gd_addition = 0.07
    ZrGd = openmc.Material.mix_materials(
        [Zr, Gd], [1 - Gd_addition, Gd_addition], "ao", name="ZrGd"
    )
    ZrGd.temperature = common_temperature
    # Zirconium Boride absorber coating

    ZrB2 = openmc.Material(name="ZrB2", temperature=common_temperature)
    ZrB2.set_density("g/cm3", 6.1)

    ZrB2.add_element("zirconium", 1.0)
    ZrB2.add_element("boron", 2.0)


    # Samarium Oxide Burnable Absorber (NAA-SR-9642, pg. 14)

    Sm2O3 = openmc.Material(name="Sm2O3", temperature=common_temperature)
    Sm2O3.set_density("g/cm3", 8.35)

    Sm2O3.add_element("samarium", 2.0)
    Sm2O3.add_element("oxygen", 3.0)


    materials = openmc.Materials(
        [
            TRIGA_fuel,
            TRIGA_fuel_BA,
            ZrH,
            NaK,
            Zr,
            SS304,
            SS316,
            Be,
            BeO,
            B4C_nat,
            B4C_rich,
            ZrHf_nat,
            ZrGd,
            ZrB2,
            Sm2O3,
        ]
    )
    materials.export_to_xml()
    
    # Materials
    fuel = TRIGA_fuel
    moderator = ZrH
    reflector = BeO
    
    
    """
    ***************************************************************************************************************************
                                                    Sec. 2 : GEOMETRY: FUEL, MODERATOR, ASSEMBLY
    ***************************************************************************************************************************
    """
    
    
    ## Fuel (these variables' names need to be reviewed!)
    fuel_radii = {'insert': params['fuel_radius'],
                'gap1': params['first_gap_radius'],
                'fuel_meat': params['fuel_meat_radius'],
                'gap2': params['second_gap_radius'],
                'cladding': params['cladding_radius']
    }

    shells = [openmc.ZCylinder(r=r) for r in fuel_radii.values()]

    region = {'insert': -shells[0],
            'gap1': +shells[0] & -shells[1],
            'fuel_meat': +shells[1] & -shells[2],
            'gap2': +shells[2] & -shells[3],
            'cladding': +shells[3] & -shells[4],
            'coolant': +shells[4]
    }

    fuel_materials = [Zr, None, fuel, None, SS304, NaK]

    fuel_cells = create_cells(region, fuel_materials)

    fuel_pin = openmc.Universe(cells=fuel_cells.values())

    # colors to be used in plotting the fuel pin
    colors = {Zr: 'green',
            SS304: 'pink',
            NaK: 'blue',
            TRIGA_fuel: 'red',
            TRIGA_fuel_BA: 'goldenrod',
            ZrH: 'orange',
            Be: 'moccasin',
            BeO: 'seagreen',
            B4C_nat: 'black'}

    
    ## Reflector
    moderator_radii = {'moderator': fuel_radii['gap2'],
                    'cladding': fuel_radii['cladding']}

    shells = [openmc.ZCylinder(r=r) for r in moderator_radii.values()]

    region = {'moderator': -shells[0],
            'cladding': +shells[0] & -shells[1],
            'coolant': +shells[1]
    }

    moderator_materials = [moderator, SS304, NaK]

    moderator_cells = create_cells(region, moderator_materials)

    moderator_pin = openmc.Universe(cells=moderator_cells.values())

    reflector_cell = openmc.Cell(fill=reflector)
    reflector_universe = openmc.Universe(cells=(reflector_cell,))

    coolant_cell = openmc.Cell(fill=NaK)
    coolant_universe = openmc.Universe(cells=(coolant_cell,))


    # Assembly
    pin_pitch = fuel_radii['cladding'] * 2 + params["pin_gap_distance"]

    lattice_radius = pin_pitch * params['assembly_rings']

    lattice_height = 2*lattice_radius
    
    
    """
    ***************************************************************************************************************************
                                                    Sec. 3 : CONTROL DRUMS
    ***************************************************************************************************************************
    """
    
    
    DRUM_RADIUS = params['drum_radius_to_lattice_radius'] * lattice_radius  # default value is 0.22784810068

    ABSORBER_THICKNESS = params['drum_Absorber_thickness']

    absorber_arc = np.pi/3
    REFERENCE_ANGLE = 0

    rotation_angle = 180

    cd_inner_shell = openmc.ZCylinder(r=DRUM_RADIUS - ABSORBER_THICKNESS)
    cd_outer_shell = openmc.ZCylinder(r=DRUM_RADIUS)

    cutting_plane_1 = openmc.Plane(a=1, b=absorber_arc/2)
    cutting_plane_2 = openmc.Plane(a=1, b=-absorber_arc/2)

    drum_absorber = +cd_inner_shell & -cd_outer_shell & -cutting_plane_1 & -cutting_plane_2
    drum_reflector = -cd_outer_shell & ~drum_absorber
    drum_outside = +cd_outer_shell
    drum_absorber = openmc.Cell(name='drum_absorber', fill=B4C_nat, region=drum_absorber)
    drum_reflector = openmc.Cell(name='drum_reflector', fill=Be, region=drum_reflector)
    drum_exterior = openmc.Cell(name='drum_outside', region=drum_outside)

    drum_reference = openmc.Universe(cells=(drum_reflector, drum_absorber, drum_exterior))

    drum_cells = []
    for r in range(0, 360, 60):
        dc = openmc.Cell(name=f'drum{r}', fill=drum_reference)
        dc.rotation = [0, 0, REFERENCE_ANGLE + r + rotation_angle]
        drum_cells.append(dc)

    drums = [openmc.Universe(cells=(dc,)) for dc in drum_cells]

    # for d in drums:
    #     d.plot(width=(DRUM_RADIUS*2, DRUM_RADIUS*2), color_by='material', colors=colors)

    drum_height  = params['drum_height_to_lattice_height'] * lattice_height   # since in MARVEL, the control drum height is 1.24* active height
    tot_drum_vol = 3.14*DRUM_RADIUS * DRUM_RADIUS *drum_height 
    drum_absorp_vol = (3.14*( DRUM_RADIUS * DRUM_RADIUS - (DRUM_RADIUS-1)*(DRUM_RADIUS-1) )*drum_height)/3
    drum_refl_vol = tot_drum_vol - drum_absorp_vol 

    tot_drum_vol_all =  tot_drum_vol * 12 
    tot_drum_area_all =  tot_drum_vol_all /drum_height

    drum_absorp_vol_all = drum_absorp_vol  * 12 
    drum_refl_vol_all = drum_refl_vol  * 12 

    drum_absorp_all_mass = drum_absorp_vol_all * 2.52/1000 # B4C (in Kg)
    drum_refl_all_mass = drum_refl_vol_all  * 3.02/1000 # BeO (in Kg)
    

    """
    ***************************************************************************************************************************
                                                    Sec. 4 : ASSEMBLY & RINGS
    ***************************************************************************************************************************
    """ 
    
    
    assembly = openmc.HexLattice()

    assembly.center = (0., 0.)
    assembly.pitch = (pin_pitch,)
    assembly.outer = coolant_universe
    
    rings = [[moderator_pin, fuel_pin, fuel_pin, fuel_pin, moderator_pin, fuel_pin, fuel_pin, moderator_pin, fuel_pin, fuel_pin, fuel_pin]*6,
         [fuel_pin, fuel_pin, moderator_pin, fuel_pin, fuel_pin, moderator_pin, fuel_pin, fuel_pin, moderator_pin, fuel_pin]*6,
         [moderator_pin, fuel_pin, fuel_pin, fuel_pin, fuel_pin, fuel_pin, fuel_pin, fuel_pin, fuel_pin]*6,
         [fuel_pin, fuel_pin, moderator_pin, fuel_pin, fuel_pin, fuel_pin, moderator_pin, fuel_pin]*6,
         [moderator_pin, fuel_pin, fuel_pin, fuel_pin, fuel_pin, fuel_pin, fuel_pin]*6,
         [fuel_pin, fuel_pin, moderator_pin, fuel_pin, moderator_pin, fuel_pin]*6,
         [moderator_pin, fuel_pin, fuel_pin, fuel_pin, fuel_pin]*6,
         [fuel_pin, fuel_pin, moderator_pin, fuel_pin]*6,
         [moderator_pin, fuel_pin, fuel_pin]*6,
         [fuel_pin, moderator_pin]*6,
         [fuel_pin]*6,
         [moderator_pin]
         ]

    rings = rings[-params['assembly_rings']:]
    assembly.universes = rings
    
    # Number of fuel elements and moderator elements
    fuel_number = sum(r.count(fuel_pin) for r in rings)
    moderator_number = sum(r.count(moderator_pin) for r in rings)
    
    
    """
    ***************************************************************************************************************************
                                                    Sec. 5 : VOLUME INFO for Depletion
    ***************************************************************************************************************************
    """ 
    
    
    fissile_area = circle_area(fuel_radii['fuel_meat']) - circle_area(fuel_radii['gap1'])

    fuel.volume = fissile_area * lattice_height * fuel_number

    heat_transfer_surface = cylinder_radial_shell(fuel_radii['cladding'], lattice_height) * fuel_number  * 1e-4 # convert from cm2 to m2

    power_MW_th = params['power_MW_th']
    power_MW_e = power_MW_th * params['thermal_efficiency']

    print(f'Average heat flux = {np.round(power_MW_th/heat_transfer_surface, 2)} MW/m2')

    if (power_MW_th/heat_transfer_surface) > 0.9:
        print(" \n WARNING : HIGHT HEAT FLUX \n")

    materials = openmc.Materials([fuel, ZrH, NaK, Zr, SS304, Be, BeO, B4C_nat])

    openmc.Materials.cross_sections = "/home/hannbn/projects/MARVEL_MRP/Github_repos/openmc_data/endfb-viii.0-hdf5/cross_sections.xml"
    materials.export_to_xml()
    
    
    assembly_boundary = openmc.model.hexagonal_prism(edge_length=pin_pitch*(params['assembly_rings']-1)+pin_pitch*0.6, corner_radius=fuel_radii['cladding'] + params["pin_gap_distance"])

    fuel_assembly_cell = openmc.Cell(fill=assembly, region=assembly_boundary)
    reflector_cell = openmc.Cell(fill=reflector, region=~assembly_boundary)

    assembly_universe = openmc.Universe(cells=[fuel_assembly_cell, reflector_cell])
    
    
    """
    ***************************************************************************************************************************
                                                    Sec. 6 : CORE DRUM REPLACEMENT
    ***************************************************************************************************************************
    """ 
    
    
    drum_gap_distance =  DRUM_RADIUS/90 # it was 0.1 and I made it as a ratio
    drum_tube_radius = DRUM_RADIUS + drum_gap_distance

    # Placement of drums happen by tracing a line through the core apothems
    # then 2 drums are place after each apothem by deviating from this line
    # by a deviation angle

    sector = np.pi/3
    deviation = (np.pi/14 )

    # i replaced the formula of Rodrigo with another one
    cd_distance = 0.86602540378 * lattice_radius  + drum_tube_radius 
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

    core_radius = lattice_radius + params['extra_reflector']

    outer_surface = openmc.ZCylinder(r=core_radius, boundary_type='vacuum')

    core_cell = openmc.Cell(fill=assembly_universe, region=-outer_surface & drums_outside)

    core_geometry = openmc.Geometry([core_cell] + drum_cells)
    core_geometry.export_to_xml()

    # Hexagon area : https://en.wikipedia.org/wiki/Hexagon
    hex_area = 2.598*lattice_radius*lattice_radius
    # I assume for now that the drums are always fully inside the reflector
    area_reflector = (3.14*core_radius*core_radius) - hex_area - tot_drum_area_all # cm2
    vol_reflector = area_reflector * drum_height # cm^3
    mass_reflector = vol_reflector * 3.02/1000
    
    
    """
    ***************************************************************************************************************************
                                                    Sec. 7 : SIMULATION
    ***************************************************************************************************************************
    """ 
    
    
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
        
        
        


        
        
        
        
        
        
        
        
        
        





        
    
 