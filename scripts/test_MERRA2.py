#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  7 16:30:47 2018

@author: cpatrizio
"""

import sys
sys.path.append("/Users/cpatrizio/repos/")
import cdms2 as cdms2
import glob
import numpy as np
from AMO.misc_fns import spatial_ave
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, shiftgrid, maskoceans
import cdutil
import matplotlib

matplotlib.rcParams.update({'font.size': 28})
matplotlib.rcParams.update({'figure.figsize': (8,6)})
matplotlib.rcParams.update({'lines.linewidth': 3})
matplotlib.rcParams.update({'legend.fontsize': 22})
matplotlib.rcParams.update({'mathtext.fontset': 'cm'})

fdata = '/Users/cpatrizio/data/MERRA2/'
fout = '/Volumes/GoogleDrive/My Drive/PhD/figures/MERRA2/'
fkernels = '/Users/cpatrizio/repos/cam5-kernels/kernels/'

#CAM5 KERNELS FOR MERRA2
falbkern = cdms2.open(fkernels + 'alb.kernel.nc')
fqkern = cdms2.open(fkernels + 'q.kernel.MERRA2plev.nc')
ftkern = cdms2.open(fkernels + 't.kernel.MERRA2plev.nc')
#ftkerngw = cdms2.open(fin + 't.kernel.nc')
ftskern = cdms2.open(fkernels + 'ts.kernel.nc')

pdiff = cdms2.open(fdata + 'dp_MERRA2plev.nc')['dp'][:]/100.

#NOTE: plev kernels have units W m^-2 hPa^-1

FLNSt = ftkern('FLNS')
#gw = ftkerngw('gw')

lats = FLNSt.getLatitude()[:]
lons = FLNSt.getLongitude()[:]
t = FLNSt.getTime()[:]
nlat = lats.size
nlon = lons.size
nt = t.size
p_Pa = ftkern('plev')
nplevs = p_Pa.size
p = np.zeros((nplevs, nlat, nlon)).T
p[:,:,:] = p_Pa/100.
p = p.T
p = np.repeat(p[np.newaxis,...], nt, axis=0)
x=np.cos(np.deg2rad(lats))
p_tropopause_zonalmean=300-200*x

p_tropopause = np.zeros((nlat, nlon)).T
p_tropopause[:,:] = p_tropopause_zonalmean
p_tropopause = p_tropopause.T
p_tropopause = np.repeat(p_tropopause[np.newaxis,...],nplevs,axis=0)
p_tropopause = np.repeat(p_tropopause[np.newaxis,...],nt,axis=0)

grid = cdms2.createGenericGrid(lats,lons)

#MERRA-2 DATA
fnamet = glob.glob(fdata + 'MERRA2_t_monthly1970to2017.nc')[0]
fnameqv = glob.glob(fdata + 'MERRA2_qv_monthly1970to2017.nc')[0]
fnametsps = glob.glob(fdata + 'MERRA2_tsps_monthly1970to2017.nc')[0]
fnamerad = glob.glob(fdata + 'MERRA2_rad_monthly1970to2017.nc')[0]

tfile = cdms2.open(fnamet)
qvfile = cdms2.open(fnameqv)
tspsfile = cdms2.open(fnametsps)
radfile = cdms2.open(fnamerad)

tinit = 12*5


#cdutil.setTimeBoundsMonthly(swdwns)
#cdutil.setTimeBoundsMonthly(swnets_cs)
#cdutil.setTimeBoundsMonthly(swdwns_cs)

#NEED t, ts, qv, AND albedo TO DO CLEAR-SKY ENERGY CLOSURE OF SEASONAL CYCLE
#ADDITIONALLY: NEED all-sky and clear-sky LW, SW radiation TO DO ALL-SKY ENERGY CLOSURE OF SEASONAL CYCLE
ta = tfile('T')[-tinit:-1,:,:,:]
cdutil.setTimeBoundsMonthly(ta)
ta= cdutil.ANNUALCYCLE.climatology(ta)

ts = tspsfile('TS')[-tinit:-1,:]
cdutil.setTimeBoundsMonthly(ts)
ts =  cdutil.ANNUALCYCLE.climatology(ts)

qv = qvfile('QV')[-tinit:-1,:]
cdutil.setTimeBoundsMonthly(qv)
qv = cdutil.ANNUALCYCLE.climatology(qv)


LW_net_surf = radfile('LWGNT')
cdutil.setTimeBoundsMonthly(LW_net_surf)
LW_net_surf = cdutil.ANNUALCYCLE.climatology(LW_net_surf)
#lwnets_cs = radfile('LWGNTCLR')
LW_net_TOA = radfile('LWTUP')
cdutil.setTimeBoundsMonthly(LW_net_TOA)
LW_net_TOA = cdutil.ANNUALCYCLE.climatology(LW_net_TOA)
#lwnettoa_cs = radfile('LWTUPCLR')
SW_net_surf = radfile('SWGNT')
#swnets_cs = radfile('SWGNTCLR')
#swnettoa = radfile('SWTNT')
#swnettoa_cs = radfile('SWTNTCLR')
SW_dwn_surf = radfile('SWGDN')
#swdwns_cs = radfile('SWGDNCLR')
#psvar = tspsfile['PS']

alb = 1 - SW_net_surf/SW_dwn_surf
cdutil.setTimeBoundsMonthly(alb)
alb = cdutil.ANNUALCYCLE.departures(alb)


merralats = ts.getLatitude()[:]
merralons = ts.getLongitude()[:]

#TODO: GET ALL SKY AND CLEAR SKY RADIATIVE VARIABLES

#TODO: MASK LAND FROM VARIABLES?

#HORIZONTALLY INTERPOLATE TO KERNEL GRID
ta = ta.regrid(grid, regridTool="esmf", regridMethod = "linear")
ts = ts.regrid(grid, regridTool="esmf", regridMethod = "linear")
qv = qv.regrid(grid, regridTool="esmf", regridMethod = "linear")
alb = alb.regrid(grid, regridTool="esmf", regridMethod = "linear")

tabar = np.ma.average(ta, axis=0)
tabar_grid = np.zeros(ta.shape)
tabar_grid[:] = tabar
tsbar = np.ma.average(ts, axis=0)
tsbar_grid = np.zeros(ts.shape)
qvbar = np.ma.average(qv, axis=0)
qvbar_grid = np.zeros(qv.shape)
albbar = np.ma.average(alb,axis=0)
albbar_grid = np.zeros(alb.shape)

ta_anom = ta - tabar_grid
ts_anom = ts - tsbar_grid
qv_anom = qv - qvbar_grid
alb_anom = alb - albbar_grid


#TODO: IMPLEMENT calcdq1k.ncl (change in moisture for +1 K perturbation at constant RH)

#TODO: DECOMPOSE ALL-SKY AND CLEAR-SKY RADIATION INTO DIFFERENT COMPONENTS (t, qv, albedo, and cloud)

#TEMPERATURE COMPONENT
dts=ts_anom
#dts_globalmean = np.ma.average(spatial_ave(dts,lats))

ts_kernel_TOA = ftskern('FLNT')
ts_kernel_surf = ftskern('FLNS')
ts_kernel_TOA_cs = ftskern('FLNTC')
ts_kernel_surf_cs = ftskern('FLNSC')

dLW_ts_TOA=ts_kernel_TOA*dts  
dLW_ts_surf = ts_kernel_surf*dts
dLW_ts_TOA_cs = ts_kernel_TOA_cs*dts
dLW_ts_surf_cs = ts_kernel_surf_cs*dts

dta = ta_anom
ta_kernel_TOA = ftkern('FLNT')
ta_kernel_surf = ftkern('FLNS')
ta_kernel_TOA_cs = ftkern['FLNTC']
ta_kernel_surf_cs = ftkern['FLNSC']

mask_temp = dta.mask
#mask stratospheric values
dta = np.ma.MaskedArray(dta, mask=np.bitwise_and(mask_temp, p>=p_tropopause))

dLW_ta_TOA = np.ma.sum(ta_kernel_TOA*dta*pdiff, axis=1)
dLW_ta_surf = np.ma.sum(ta_kernel_surf*dta*pdiff, axis=1)
dLW_ta_TOA_cs = np.ma.sum(ta_kernel_TOA_cs*dta*pdiff, axis=1)
dLW_ta_surf_cs = np.ma.sum(ta_kernel_surf_cs*dta*pdiff, axis=1)

dLW_t_TOA_globalmean = spatial_ave(-dLW_ta_TOA-dLW_ts_TOA, lats)
dLW_t_surf_globalmean = spatial_ave(-dLW_ta_surf-dLW_ts_surf, lats)

f, axarr = plt.subplots(2,1, figsize=(8,6))
axarr[0,].plot(dLW_t_TOA_globalmean)
axarr[0,].set_title(r'$\Delta$ LW$_{TOA,T}$')
axarr[0,].set_ylabel('W m$^{-2}$')
axarr[1,].plot(dLW_t_surf_globalmean)
axarr[1,].set_title(r'$\Delta$ LW$_{surf,T}$')
axarr[1,].set_xlabel('time (month)')
axarr[1,].set_ylabel('W m$^{-2}$')
plt.show()




















#1970-2017 MEAN T_S 
#par = np.arange(-90.,91.,30.)
#mer = np.arange(-180.,181.,60.)
#tsmin = 260
#tsmax = 310
#nlevs = 100
#tslevels = np.linspace(tsmin,tsmax,nlevs)
#
#fig =plt.figure(1, figsize=(8,6))
#ax = fig.gca()
#m = Basemap(projection='moll',lon_0=180,resolution='i')
#m.drawcoastlines(linewidth=0.1)
#m.drawparallels(par, dashes=[1,0.001], labels=[1,0,0,1], linewidth=0.1)
#m.drawmeridians(mer, dashes=[1,0.001], linewidth=0.1)
##ts_anncycle, lonsout = shiftgrid(180,ts_anncycle_regr,lons)
#x, y = m(*np.meshgrid(lons, lats))
#m.contourf(x, y, np.average(ts, axis=0), cmap=plt.cm.RdPu, levels=tslevels)
#m.colorbar(label=r'K', format='%3.1f')
#ax.set_title(r'$T_s$')
#plt.savefig(fout + 'MERRA2_ts_maps_test.pdf')
#plt.close()














