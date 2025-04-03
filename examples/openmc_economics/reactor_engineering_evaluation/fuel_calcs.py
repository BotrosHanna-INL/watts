import numpy as np
def fuel_calculations(params):
    
    U_mass = (params['mass_U235'] + params['mass_U238']) / 1000 # The mass of Uranium only (in Kg)
    nat_u_consum = U_mass*(params['enrichment'] -0.0025)/(0.0071-0.0025) # Kg
    tail_waste = nat_u_consum - U_mass # Kg

    # value functions
    f_val_fun = (1-2*params['enrichment'])*np.log((1-params['enrichment'])/params['enrichment'])
    tail_waste_val_fun = 5.96
    nat_u_waste_val_fun = 4.87


    kg_SWU = (U_mass*f_val_fun+tail_waste*tail_waste_val_fun- nat_u_consum *nat_u_waste_val_fun)
    
    if  0 <params['enrichment'] < 0.1:
        premium = 1
        print(f"HALEU premium is {premium}")
        
    if 0.1<= params['enrichment'] <=   0.2:
        premium = 1.15 
        return nat_u_consum , tail_waste , kg_SWU
 