# Importing libraries
import openmc


def collect_materials_data(params):
    
    # **************************************************************************************************************************
    #                                               Sec. 1 : MATERIALS
    # **************************************************************************************************************************
    
    
    # """""""""""""""""""""
    # Sec. 1.1 : TRIGA Fuel
    # """""""""""""""""""""

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
    ZrH_fuel.add_nuclide("H1", params["H_Zr_ratio"]) # The proportion of hydrogen atoms relative to zirconium atoms


    # Now we mix the components in the right weight %
    # The NRC seems to license up to TRIGA fuel with up to 45% weight Uranium

    TRIGA_fuel = openmc.Material.mix_materials(
        [U_met, ZrH_fuel], [params['U_met_wo'], 1 - params['U_met_wo']], "wo", name="UZrH"
    )

    TRIGA_fuel.temperature = params['common_temperature']
    TRIGA_fuel.add_s_alpha_beta("c_H_in_ZrH")

    # Let's also make a version with 3% Erbium in the meat
    Er = openmc.Material(name="Er", temperature= params['common_temperature'])
    Er.set_density("g/cm3", 9.2)
    Er.add_element("erbium", 1.0)



    # """""""""""""""""""""
    # Sec. 1.2 : Zirconium Hydride 
    # """""""""""""""""""""
       
    ZrH = openmc.Material(name="ZrH", temperature= params['common_temperature'])
    ZrH.set_density("g/cm3", 5.6)

    ZrH.add_nuclide("H1", 1.85)
    ZrH.add_element("zirconium", 1.0)
    ZrH.add_s_alpha_beta("c_H_in_ZrH")

    # """""""""""""""""""""
    # Sec. 1.3 : NaK (Coolant)
    # """""""""""""""""""""
    
    NaK = openmc.Material(name="NaK", temperature= params['common_temperature'])
    NaK.set_density("g/cm3", 0.75)
    NaK.add_nuclide("Na23", 2.20000e-01)
    NaK.add_nuclide("K39", 7.27413e-01)
    NaK.add_nuclide("K41", 5.24956e-02)

    # """""""""""""""""""""
    # Sec. 1.4 : Beryllium 
    # """""""""""""""""""""
    Be = openmc.Material(name="Be")
    Be.add_element("beryllium", 1.0)
    Be.add_s_alpha_beta("c_Be")
    Be.set_density("g/cm3", 1.84)
    Be.temperature =  params['common_temperature']
    
    
    # """""""""""""""""""""
    # Sec. 1.5 : Beryllium Oxide
    # """""""""""""""""""""
    
    BeO = openmc.Material(name="BeO", temperature= params['common_temperature'])
    BeO.set_density("g/cm3", 3.01)
    BeO.add_element("beryllium", 1.0)
    BeO.add_element("oxygen", 1.0)
    BeO.add_s_alpha_beta("c_Be_in_BeO")


    # """""""""""""""""""""
    # Sec. 1.6 : Zirconium
    # """""""""""""""""""""
    
    Zr = openmc.Material(name="Zr", temperature= params['common_temperature'])
    Zr.set_density("g/cm3", 6.49)
    Zr.add_element("zirconium", 1.0)
    
    # """""""""""""""""""""
    # Sec. 1.7 : SS304
    # """""""""""""""""""""
    
    SS304 = openmc.Material(name="SS304", temperature= params['common_temperature'])
    SS304.set_density("g/cm3", 7.98)
    SS304.add_element("carbon", 0.04)
    SS304.add_element("silicon", 0.50)
    SS304.add_element("phosphorus", 0.023)
    SS304.add_element("sulfur", 0.015)
    SS304.add_element("chromium", 19.00)
    SS304.add_element("manganese", 1.00)
    SS304.add_element("iron", 70.173)
    SS304.add_element("nickel", 9.25)

  
    # """""""""""""""""""""
    # Sec. 1.8 : Boron Carbide
    # """""""""""""""""""""
    B4C_nat = openmc.Material(name="B4C", temperature= params['common_temperature'])
    B4C_nat.add_element("boron", 4)
    B4C_nat.add_element("carbon", 1)
    B4C_nat.set_density("g/cm3", 2.52)

    materials = openmc.Materials(
        [
            TRIGA_fuel,
            ZrH,
            NaK,
            Zr,
            SS304,
            Be,
            BeO,
            B4C_nat,
        ]
    )
    
    # I am not sure if this is necessary but keeping it for now
    materials.export_to_xml()
    
    return {"Zr": Zr, "TRIGA_fuel": TRIGA_fuel,\
        "SS304": SS304, "NaK": NaK, "ZrH": ZrH,\
            "Be": Be, "BeO": BeO, "B4C_nat": B4C_nat}
    
    
    
def color_materials():
    colors = {'Zr': 'green',
            'SS304': 'pink',
            'NaK': 'blue',
            'TRIGA_fuel': 'red',
            'ZrH': 'orange',
            'Be': 'moccasin',
            'BeO': 'seagreen',
            'B4C_nat': 'black'}
    return colors

def extract_list_lof_materials_properties(params, materials_list):
    material_properties = collect_materials_data(params)
    materials_properties_list = []
    # the input is a list like ("Zr", "SS")
    for item in materials_list:
        if item == None:
            materials_properties_list.append(None)
        else:    
            
            material_properties = (collect_materials_data(params))[item]
            materials_properties_list.append(material_properties)
        
    return materials_properties_list 
    