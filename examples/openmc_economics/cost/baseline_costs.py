
from cost.cost_utils import adjust_for_inflation

# **************************************************************************************************************************
#                                                Sec. 1 : Baseline  Fixed Costs (dollars)
# **************************************************************************************************************************

### Sec. 1.1 Costs From MARVEL


### Sec. 1.2 Costs From other sources



# **************************************************************************************************************************
#                                                Sec. 2 : Baseline  Unit Costs (dollars per unit)
# **************************************************************************************************************************


### Sec. 2.1 Unit Costs From MARVEL



### Sec. 2.2 Unit Costs From other sources

# from the EBD report
#https://inldigitallibrary.inl.gov/sites/sti/sti/Sort_46104.pdf
land_unit_cost = ("Land Unit Cost", 3800, 2021) # $/acre  in 2021


    


# **************************************************************************************************************************
#                                                Sec. 3 : Add all the costs to the cost dictionary
# **************************************************************************************************************************

# add all the costs to the cost dictionary
costs_list = [land_unit_cost]

adjusted_cost_list = []
for item in costs_list:
    item_adjusted_cost = adjust_for_inflation(item)
adjusted_cost_list.append(item_adjusted_cost)
cost_dictionary_data = dict(adjusted_cost_list)


