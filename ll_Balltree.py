import math
import numpy as np
from sklearn import preprocessing
from sklearn.neighbors import BallTree

#earth radius in km 
EARTH_RADIUS_KM = 6371


"""
params:
    idx_mask - mask indexed by (lat, lon)
    lats - list of latitude values
    lons - list of longitude values
    
returns:
    list of lat and lon for the masked values
"""
def get_cells(idx_mask, lats, lons):
    if idx_mask is None:
        return np.array([lats, lons])
    
    cell_lats = lats[idx_mask[0]]
    cell_lons = lons[idx_mask[1]]
    return np.array([cell_lats, cell_lons], dtype = np.float64)

"""
params:
    cells - list of lat and lon values
returns:
    BallTree object
"""
def build_ll_BallTree(cells):
    cells_rad = (np.transpose(cells) * np.pi) / 180

    return BallTree(cells_rad, leaf_size = 2, metric = 'haversine') 

def query_ll_BallTree(tree, cells):
    cells_rad  = (np.transpose(cells) * np.pi) / 180 

    dist, ind = tree.query(cells_rad, k = 1)
    
    return dist.T[0]*EARTH_RADIUS_KM, ind.T[0]