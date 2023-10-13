
import xarray as xr
import pandas as pd
import numpy as np
import cftime
import cartopy

from pylab import *
from numpy import *
from glob import glob
from os import path
import time
import dask
import dask.array as da
from dask.diagnostics import ProgressBar

import sys
args = sys.argv[1:]

#specify the directory in which the MPI code wrote its output
inputDir=str(args[0])

#chunksize and where the output zarr file should go; also set chunksize of output file
chunksize={'trajectory':5*int(1e4),'obs':10};
outputDir=str(args[1])

#do not overwrite existing data sets
if path.exists(outputDir):
    print('the ouput path',outputDir,'exists')
    print('please delete if you want to replace it')
    assert False,'stopping execution'
    
varType = {
         'lat':dtype('float32'),
         'lon':dtype('float32'),
          'time':dtype('float64'), #to avoid bug in xarray
         'z':dtype('float32'),
          }


#dataIn = xr.open_zarr(inputDir, decode_times=False)

print('opening data from multiple process files')
tic=time.time()
files = glob(path.join(inputDir, "*.zarr"));
print(files)
dataIn = xr.concat([xr.open_zarr(f,decode_times=False) for f in files], dim='trajectory',
                     compat='no_conflicts', coords='minimal')
print('   done opening in %5.2f'%(time.time()-tic))

for v in varType.keys():
    dataIn[v]=dataIn[v].astype(varType[v])
    
print('re-chunking')
tic=time.time()
for v in dataIn.variables:
    if 'chunks' in dataIn[v].encoding:
        del dataIn[v].encoding['chunks']
dataIn=dataIn.chunk(chunksize)
print('   done in',time.time()-tic)

delayedObj=dataIn.to_zarr(outputDir,compute=False)
with ProgressBar():
        results=delayedObj.compute()
