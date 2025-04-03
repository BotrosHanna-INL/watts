
import numpy as np 

def refueling_operation(params):
    #this function returns how many people are required for the refueling process multiplied by how many days of refueling
    
    # how many times you add the fuel over the entire reactor lifetime
    add_fuel_num = int(np.floor(365*params['levelization_period_years']/ 
                                (params['refueling_period_days'] + params['fuel_lifetime_days'])))

    num_of_refuel_days_per_year = params['refueling_period_days'] *\
        add_fuel_num/params['levelization_period_years']
    
    people_by_days_refueling = params['num_people_required_per_refueling'] * num_of_refuel_days_per_year
    return people_by_days_refueling      
