#import sys
#sys.path.append("/Users/cpatrizio/repos/")
from netCDF4 import Dataset
import glob
import numpy as np
import AMO.misc_fns
from AMO.misc_fns import spatial_ave
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

fkernels = '/Volumes/GoogleDrive/My Drive/data/cam5-kernels/kernels/'
fdata = '/Volumes/GoogleDrive/My Drive/data/cam5-kernels/MERRA2/
fout = '/Volumes/GoogleDrive/My Drive/PhD/figures/CAM5kernels/'

#CAM5 KERNELS FOR MERRA2
falbkern = Dataset(fkernels + 'alb.kernel.nc')
fqkern = Dataset(fkernels + 'q.kernel.MERRA2plev.nc')
ftkern = Dataset(fkernels + 't.kernel.MERRA2plev.nc')
#ftkerngw = Dataset.open(fin + 't.kernel.nc')
ftskern = Dataset(fkernels + 'ts.kernel.nc')

pdiff = Dataset(fdata + 'dp_MERRA2plev.nc')['dp'][:]/100.

#NOTE: plev kernels have units W m^-2 hPa^-1

#FLNSt = ftkern('FLNS')
#gw = ftkerngw('gw')

lats = falbkern['lat'][:]
lons = falbkern['lon'][:]
t = falbkern['time'][:]
nlat = lats.size
nlon = lons.size
nt = t.size
p_Pa = ftkern['plev'][:]
nplevs = p_Pa.size
p = np.zeros((nplevs, nlat, nlon)).T
p[:,:,:] = p_Pa/100.
p = p.T
p = np.repeat(p[np.newaxis,...], nt, axis=0)

#grid = Dataset.createGenericGrid(lats,lons)

#MERRA-2 DATA
fnamet = glob.glob(fdata + 'MERRA2_t_monthly1970to2017.nc')[0]
fnameqv = glob.glob(fdata + 'MERRA2_qv_monthly1970to2017.nc')[0]
fnametsps = glob.glob(fdata + 'MERRA2_tsps_monthly1970to2017.nc')[0]
fnamerad = glob.glob(fdata + 'MERRA2_rad_monthly1970to2017.nc')[0]

tfile = Dataset(fnamet)
#qvfile = Dataset(fnameqv)
#tspsfile = Dataset(fnametsps)
#radfile = Dataset(fnamerad)

tinit=12*1

t = tfile['T'][-tinit:-1,:]
#qv = qvfile['QV'][:]
#ts = tspsfile['TS'][:]