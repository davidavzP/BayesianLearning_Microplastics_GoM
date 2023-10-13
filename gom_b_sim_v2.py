'''
@David Pojunas

Runs backwards simulations
'''

import xarray as xr
import numpy as np
import pandas as pd
import math 
import time
import sys
import os

import parcels.rng as ParcelsRandom
from parcels import (FieldSet, Field, ParticleSet, JITParticle, ErrorCode, Variable ,GeographicPolar,Geographic, VectorField)


from glob import glob
from datetime import timedelta 

import warnings
warnings.filterwarnings("ignore")


######################
### SET FIELDSETS  ###
######################

data_dir = "data/input/"
cs = {'time': ('time', 1), 'lat': ('latitude', 256), 'lon': ('longitude', 256)}

args = sys.argv[1:]
outfile = args[0]
release_year = int(args[1])
end_month = float(args[2])
start_month = float(args[3])
advection_time = float(args[4])


def get_hycom_fieldset(indices = {}):
    hycom_files = sorted(glob(data_dir + 'hycom/' + '*.nc'))    
    filenames = {'U': hycom_files, 'V': hycom_files}
    variables = {'U': 'water_u', 'V': 'water_v'}
    dimensions = {'lat': 'lat', 'lon': 'lon', 'time': 'time'}
    
    return FieldSet.from_netcdf(filenames, variables, dimensions, allow_time_extrapolation = True)

def get_stokes(indices = {}):
    stokes_files = sorted(glob(data_dir + 'stokes/' + '*_GOM_STOKES_raw.nc'))
    filenames = {'Ust': stokes_files, 'Vst': stokes_files}
    variables = {'Ust': 'VSDX', 'Vst': 'VSDY'}
    dimensions = {'lat': 'latitude', 'lon': 'longitude', 'time': 'time'}
    interp_method = {'Ust' : 'linear','Vst' : 'linear'}
    return FieldSet.from_netcdf(filenames, variables, dimensions, allow_time_extrapolation = True)
    
def set_stokes_sumfieldset(fieldset):
    stokes_fieldset = get_stokes()
    stokes_fieldset.Ust.units = GeographicPolar()
    stokes_fieldset.Vst.units = Geographic()
    fieldset = FieldSet(U=fieldset.U + stokes_fieldset.Ust,
                        V=fieldset.V + stokes_fieldset.Vst)
    return fieldset

def set_stokes_fieldset(fieldset):
    stokes_fieldset = get_stokes()
    stokes_fieldset.Ust.units = GeographicPolar()
    stokes_fieldset.Vst.units = Geographic()
    
    fieldset.add_field(stokes_fieldset.Ust)
    fieldset.add_field(stokes_fieldset.Vst)
    
    vectorfield_stokes = VectorField('UVst', fieldset.Ust, fieldset.Vst)
    fieldset.add_vector_field(vectorfield_stokes)
    return fieldset

def set_displacement_field(fieldset):
    # ADD BOUNDRY PROPERTIES
    #fieldset.interp_method = {'U': 'freeslip', 'V': 'freeslip'}
    
    # ADD DISPLACEMENT FIELD PROPERTIES
    gom_masks = xr.open_dataset(data_dir + 'gom_masks_w_inputs.nc')
    u_displacement = gom_masks.disp_vx.values
    v_displacement = gom_masks.disp_vy.values
    d2s = gom_masks.d2s.values

    fieldset.add_field(Field('dispU', data=u_displacement,
                            lon=fieldset.U.grid.lon, lat=fieldset.U.grid.lat,
                            mesh='spherical', allow_time_extrapolation = True))
    fieldset.add_field(Field('dispV', data=v_displacement,
                            lon=fieldset.U.grid.lon, lat=fieldset.U.grid.lat,
                            mesh='spherical', allow_time_extrapolation = True))
    fieldset.dispU.units = GeographicPolar()
    fieldset.dispV.units = Geographic()
    
    fieldset.add_field(Field('distance2shore', d2s,
                            lon=fieldset.U.grid.lon, lat=fieldset.U.grid.lat, 
                            mesh='spherical', allow_time_extrapolation = True))
    return fieldset

def set_smagdiff_fieldset(fieldset, diff = 0.1):
    fieldset.add_field(Field(name='cell_areas', data=fieldset.U.cell_areas(), lon=fieldset.U.grid.lon, lat=fieldset.U.grid.lat))
    fieldset.add_constant('Cs', diff)
    return fieldset 

#########################
### LOAD PARTICLE SET ###
#########################

def monte_carlo_multi_pr(vals, n, size = CELL_SIZE, seed = 1001):
    #np.random.seed(seed)
    assert(len(vals) == len(n)) # sanity check
    release_locs = np.repeat(vals, n, axis = 0)
    lls = release_locs[:, :2]
    r = size * np.sqrt(np.random.rand(np.shape(lls)[0], np.shape(lls)[1]))
    theta = 2 * math.pi * (np.random.rand(np.shape(lls)[0]))
    theta = np.array([np.sin(theta), np.cos(theta)]).T
    lls = (lls + theta * r).T
    # RETURNS LATS, LONS, TIME
    return lls[0], lls[1], release_locs.T[2]

class Nurdle(JITParticle):
    dU = Variable('dU')
    dV = Variable('dV')
    d2s = Variable('d2s', initial=1e3)
    age = Variable('age', initial=0.0)
    # beached : 0 sea, 1 beached, 2 after non-beach dyn, 3 after beach dyn, 4 please unbeach
    beached = Variable('beached', initial = 0.0)
    unbeachCount = Variable('unbeachCount', initial=0.0)

def get_particle_set(fieldset):
    # load nurdle dataset
    df_all_nurdles = pd.read_csv(data_dir + 'nurdle_release.csv')
    df_all_nurdles['patrol_date'] = pd.to_datetime(df_all_nurdles['patrol_date'])
    df_all_nurdles.iloc[df_all_nurdles[df_all_nurdles.nurdle_count > 10000].index, 4] = 10000 

    #find the year of release
    df_all_nurdles = df_all_nurdles[(df_all_nurdles.patrol_date.dt.year == release_year) & (df_all_nurdles.patrol_date.dt.month >= start_month) & (df_all_nurdles.patrol_date.dt.month <= end_month)]
    print(df_all_nurdles.groupby([df_all_nurdles.patrol_date.dt.year, df_all_nurdles.patrol_date.dt.month]).nurdle_count.sum())
    print(df_all_nurdles.nurdle_count.sum())
    
    vals = df_all_nurdles[['lat', 'lon', 'time']].values
    n = df_all_nurdles.nurdle_count.values
    lats, lons, time = monte_carlo_multi_pr(vals, n)
    return ParticleSet.from_list(fieldset = fieldset, pclass = Nurdle, lon = lons, lat = lats, time = time,)

###############
### KERNELS ###
###############
def AdvectionRK4(particle, fieldset, time):
    if particle.beached == 0:
        (u1, v1) = fieldset.UV[time, particle.depth, particle.lat, particle.lon]
        lon1, lat1 = (particle.lon + u1*.5*particle.dt, particle.lat + v1*.5*particle.dt)

        (u2, v2) = fieldset.UV[time + .5 * particle.dt, particle.depth, lat1, lon1]
        lon2, lat2 = (particle.lon + u2*.5*particle.dt, particle.lat + v2*.5*particle.dt)

        (u3, v3) = fieldset.UV[time + .5 * particle.dt, particle.depth, lat2, lon2]
        lon3, lat3 = (particle.lon + u3*particle.dt, particle.lat + v3*particle.dt)

        (u4, v4) = fieldset.UV[time + particle.dt, particle.depth, lat3, lon3]
        particle.lon += (u1 + 2*u2 + 2*u3 + u4) / 6. * particle.dt
        particle.lat += (v1 + 2*v2 + 2*v3 + v4) / 6. * particle.dt
        particle.beached = 2
        
def StokesUV(particle, fieldset, time):
    if particle.beached == 0:
        (u_uss, v_uss) = fieldset.UVst[time, particle.depth, particle.lat, particle.lon]
        particle.lon += u_uss * particle.dt
        particle.lat += v_uss * particle.dt
        particle.beached = 3

def SmagDiffBeached(particle, fieldset, time):
    if particle.beached == 0:
        dx = 0.01
        # gradients are computed by using a local central difference.
        updx, vpdx = fieldset.UV[time, particle.depth, particle.lat, particle.lon+dx]
        umdx, vmdx = fieldset.UV[time, particle.depth, particle.lat, particle.lon-dx]
        updy, vpdy = fieldset.UV[time, particle.depth, particle.lat+dx, particle.lon]
        umdy, vmdy = fieldset.UV[time, particle.depth, particle.lat-dx, particle.lon]

        dudx = (updx - umdx) / (2*dx)
        dudy = (updy - umdy) / (2*dx)
        
        dvdx = (vpdx - vmdx) / (2*dx)
        dvdy = (vpdy - vmdy) / (2*dx)

        A = fieldset.cell_areas[time, 0, particle.lat, particle.lon]
        sq_deg_to_sq_m = (1852*60)**2*math.cos(particle.lat*math.pi/180)
        A = A / sq_deg_to_sq_m
        Kh = fieldset.Cs * A * math.sqrt(dudx**2 + 0.5*(dudy + dvdx)**2 + dvdy**2)
        
        dlat = ParcelsRandom.normalvariate(0., 1.) * math.sqrt(2*math.fabs(particle.dt)* Kh) 
        dlon = ParcelsRandom.normalvariate(0., 1.) * math.sqrt(2*math.fabs(particle.dt)* Kh) 

        particle.lat += dlat
        particle.lon += dlon
        
        particle.beached = 3

def BeachTesting(particle, fieldset, time):
    if particle.beached == 2 or particle.beached == 3:
        (u, v) = fieldset.UV[time, particle.depth, particle.lat, particle.lon]
        if u == 0 and v == 0:
            if particle.beached == 2:
                particle.beached = 4
            else:
                dispUab, dispVab = fieldset.dispU[time, particle.depth, particle.lat,particle.lon], fieldset.dispV[time, particle.depth, particle.lat,particle.lon]
                dtt = -1*particle.dt
                particle.lon += dispUab*dtt
                particle.lat += dispVab*dtt
                particle.beached = 1    
        else:
            particle.beached = 0

def Unbeach(particle, fieldset, time):    
    if particle.beached == 4:
        dispUab, dispVab = fieldset.dispU[time, particle.depth, particle.lat,particle.lon], fieldset.dispV[time, particle.depth, particle.lat,particle.lon]
        dtt = -1*particle.dt
        particle.lon += dispUab*dtt
        particle.lat += dispVab*dtt
        particle.beached = 0

def Ageing(particle, fieldset, time):
    if particle.age < -15552000.0:
        particle.delete()
    particle.age += particle.dt

def DeleteParticle(particle, fieldset, time):
    particle.delete()

def OutOfBounds(particle, fieldset, time): 
    particle.lon = 0
    particle.lat = 0
    #if particle.beached = 5, it is not advected anymore by any kernel and is thus frozen
    particle.beached = 5
    
######################
### RUN SIMULATION ###
######################
def run_gom_mp_backwards(outdir, release_year = None, end_month = None,  start_month = None, advection_time = None, sim_dt_hours = 2, output_dt_hours = 24, fw = -1):
 
    outfile = os.path.join(outdir, "{}".format('GOM' + '_rt_' + str(release_year)+ '_' + str(end_month) + '-' + str(start_month) + '_at_' + str(advection_time)))

    # FIELD SET
    fieldset = get_hycom_fieldset()
    fieldset = set_stokes_fieldset(fieldset)
    fieldset = set_displacement_field(fieldset)
    fieldset = set_smagdiff_fieldset(fieldset)
    fieldset.field_chunksize = cs
    
    # PARTICLE SET
    pset = get_particle_set(fieldset)

            
    # KERNELS
    kernels = (pset.Kernel(AdvectionRK4) + pset.Kernel(BeachTesting) + pset.Kernel(Unbeach) + 
               pset.Kernel(StokesUV) + pset.Kernel(BeachTesting) + 
               pset.Kernel(SmagDiffBeached) + 
               pset.Kernel(Ageing) + pset.Kernel(BeachTesting))
    
    # Release Time
    release_time = (30 * ((end_month-start_month)+ 1))

    pfile = pset.ParticleFile(name=outfile, outputdt=timedelta(hours=output_dt_hours))
    pset.execute(kernels, 
                 runtime=timedelta(days=release_time + advection_time), 
                 dt=fw*timedelta(hours=sim_dt_hours), 
                 output_file=pfile, 
                 recovery={ErrorCode.ErrorOutOfBounds: OutOfBounds})
    pfile.close()
    
    return fieldset, pset

run_gom_mp_backwards(outfile, release_year = release_year, end_month = end_month,  start_month = start_month, advection_time = advection_time)