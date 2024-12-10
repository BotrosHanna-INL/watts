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