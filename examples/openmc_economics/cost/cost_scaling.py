import pandas as pd
# Base Costs




# Unit Costs



#Cost Scaling methods
def cost_scale_general(db, account_number, old_var, new_var, scaling_exponent, old_cost, constant_cost):
    if old_var == None and new_var == None and old_cost == None:
        new_cost = constant_cost
    else:    
        new_cost == old_cost * ((new_var/old_var)**scaling_exponent) + constant_cost
       #add new cost to the table
    db.loc[db.Account == account_number, 'Estimated Cost (2023 $)'] =  new_cost      
    
def scale_unit_costs(db, account_number, unit_cost, new_var):
    new_cost =   unit_cost * new_var
    db.loc[db.Account == account_number, 'Estimated Cost (2023 $)'] =  new_cost  

def add_customized_cost_estimate(db, account_number, customized_cost_estimate):
    db.loc[db.Account == account_number, 'Estimated Cost (2023 $)'] =  customized_cost_estimate 
    # print ( f"The {db.loc[db.Account == account_number, 'Account Title'].values[0] } account has been updated assuming {scaling_method} scaling \n")
    
    
    

def cost_estimate(cost_ref) :
      #Create new database
    db = pd.DataFrame()
    cost_new_reactor = cost_ref[['Account', 'Account Title', 'Estimated Cost (2023 $)']].copy()  


