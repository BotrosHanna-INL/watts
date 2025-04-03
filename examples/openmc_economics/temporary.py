

def extra_design_specs(enrichment, core_radius, Fuel_Mass_kg, U_mass, fuel_lifetime_days, refueling_period,\
    levelization_period_years, people_per_refueling, FTE_cost, duration_to_startup)   : 
    
    # how many times you add the fuel
    add_fuel_num = int(np.floor(365*levelization_period_years / (refueling_period + fuel_lifetime_days)))

    num_of_refuel_days_per_year = refueling_period *add_fuel_num/levelization_period_years
    
    People_by_days_refueling = people_per_refueling * num_of_refuel_days_per_year# per year
    refueling_cost_per_year = (People_by_days_refueling /365)*FTE_cost 


    Capacity_factor  = 1 - ((num_of_refuel_days_per_year + duration_to_startup )/365)
    people_by_days_startup = people_per_startup *duration_to_startup
    startup_cost_per_year_after_emergency_shutdown = FTE_cost  * people_by_days_startup /365

    days_per_year_to_refuel = 365*refueling_period / (refueling_period + fuel_lifetime_days)
    people_by_days_to_startup_after_refuel = people_per_startup_after_refuel * days_per_year_to_refuel 
    startup_cost_per_year_after_refuel = FTE_cost  * people_by_days_to_startup_after_refuel  /365

    

    monitor_cost = 5 * FTE_cost /reactors_monitored_by_one_person 
    sec_cost = 5 *FTE_cost /2
    
    
    startup_cost_per_year =  startup_cost_per_year_after_refuel  + startup_cost_per_year_after_emergency_shutdown 
    
    nat_u_consum = U_mass*(enrichment -0.0025)/(0.0071-0.0025) # Kg
    tail_waste = nat_u_consum - U_mass # Kg

    # value functions
    f_val_fun = (1-2*enrichment)*np.log((1-enrichment)/enrichment)
    tail_waste_val_fun = 5.96
    nat_u_waste_val_fun = 4.87


    kg_SWU = (U_mass*f_val_fun+tail_waste*tail_waste_val_fun- nat_u_consum *nat_u_waste_val_fun)# Check the units here! everything should be Kg
    
    if  0 <enrichment < 0.1:
        premium = 1
        print(f"HALEU premium is {premium}")
        
    if 0.1<=enrichment <=   0.2:
        premium = 1.15 
 
        
    

    vessel_tot_height = (vessel_calcs (core_radius , lattice_height, boron_carbide_shield_thickness, extra_reflector ))[0] # cm
    guard_vessel_mass = (vessel_calcs (core_radius , lattice_height, boron_carbide_shield_thickness, extra_reflector ))[1]
    vessel_mass = (vessel_calcs (core_radius , lattice_height, boron_carbide_shield_thickness, extra_reflector ))[2]
    vessel_height = (vessel_calcs (core_radius , lattice_height, boron_carbide_shield_thickness, extra_reflector ))[3]
    vessels_tot_radius = (vessel_calcs (core_radius , lattice_height, boron_carbide_shield_thickness, extra_reflector ))[4]
    cooling_vessel_mass = (vessel_calcs (core_radius , lattice_height, boron_carbide_shield_thickness, extra_reflector ))[5]
    intake_vessel_mass = (vessel_calcs (core_radius , lattice_height, boron_carbide_shield_thickness, extra_reflector ))[6]



    B4C_shield_vol = 3.14 *(( core_radius + boron_carbide_shield_thickness )**2  - core_radius**2 ) * vessel_height
    B4C_shield_mass =   B4C_shield_vol  * 2.52 / 1000
    
    WEP_shield_vol = 3.14 * (((WEP_shield_thickness + vessels_tot_radius)**2) -  vessels_tot_radius**2 ) * vessel_tot_height
    
    WEP_shield_mass = (WEP_shield_vol  * 1.1 /2)/1000 # divide by 2 because it is not fully WEP (Kg)  
    
    return refueling_cost_per_year , Capacity_factor , startup_cost_per_year_after_emergency_shutdown, startup_cost_per_year_after_refuel,\
        monitor_cost, sec_cost, startup_cost_per_year, nat_u_consum , tail_waste , kg_SWU, premium, \
            guard_vessel_mass, vessel_mass, cooling_vessel_mass, intake_vessel_mass, vessel_tot_height, B4C_shield_mass, WEP_shield_mass, add_fuel_num









def ellipsoid_shell(a, b, c):

    return 4*np.pi*np.power(((a*b)**1.6 + (a*c)**1.6 + (b*c)**1.6)/3, 1/1.6)

# Vessel Calcs
def vessel_calcs (core_radius_including_reflector , lattice_height, boron_carbide_shield_thickness, extra_reflector ):



    vessel_radius = core_radius_including_reflector + boron_carbide_shield_thickness+ 20# I added 10 cm (4 inches) for the boron carbide in addition to the 20 cm that were already there
    vessel_thickness = 2

    lower_plenum_height = 30
    upper_plenum_height = 60
    upper_gas_gap = 10
    bottom_depth = 35
 
    vessel_height = lattice_height + 2.5*extra_reflector + lower_plenum_height + upper_plenum_height + upper_gas_gap # This is the first vessel
    vessel_volume = (ellipsoid_shell(vessel_radius, vessel_radius, bottom_depth)/2)*vessel_thickness + (circle_area(vessel_radius + vessel_thickness) - circle_area(vessel_radius))*vessel_height
    vessel_mass_kg = vessel_volume  * 8/1000

    gap_vessel = 2 # cm
    guard_vessel_thickness = 0.5
    guard_vessel_radius = vessel_radius + vessel_thickness + gap_vessel
    guard_bottom_depth = bottom_depth + vessel_thickness + gap_vessel
    guard_vessel_volume = (ellipsoid_shell(guard_vessel_radius, guard_vessel_radius, guard_bottom_depth)/2)*guard_vessel_thickness + (circle_area(guard_vessel_radius + guard_vessel_thickness) - circle_area(guard_vessel_radius))*vessel_height
    guard_vessel_mass_kg = guard_vessel_volume * 8/1000

    # cooling vessel
    gap_cooling = 5 # cm
    cooling_vessel_thickness = 0.3
    cooling_vessel_radius = guard_vessel_radius + gap_cooling # cm
    cooling_bottom_depth = guard_bottom_depth + guard_vessel_thickness + gap_cooling
    cooling_vessel_volume = (ellipsoid_shell(cooling_vessel_radius, cooling_vessel_radius, cooling_bottom_depth)/2)*cooling_vessel_thickness + (circle_area(cooling_vessel_radius + cooling_vessel_thickness) - circle_area(cooling_vessel_radius))*vessel_height
    cooling_vessel_mass = cooling_vessel_volume *8/1000
    
    gap_intake = 3

    # This is the intake vessel
    intake_vessel_thickness = 0.3
    intake_vessel_radius = cooling_vessel_radius + gap_intake
    intake_bottom_depth = cooling_bottom_depth + cooling_vessel_thickness + gap_intake
    intake_vessel_volume = (ellipsoid_shell(intake_vessel_radius, intake_vessel_radius, intake_bottom_depth)/2)*intake_vessel_thickness + (circle_area(intake_vessel_radius + intake_vessel_thickness) - circle_area(intake_vessel_radius))*vessel_height
    intake_vessel_mass = intake_vessel_volume *8/1000
    
    total_vessel_height = intake_bottom_depth + vessel_height
    full_radius = intake_vessel_radius + intake_vessel_thickness
    
    return total_vessel_height, guard_vessel_mass_kg, vessel_mass_kg, vessel_height, full_radius,  cooling_vessel_mass, intake_vessel_mass












    