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
from ocean_atmosphere.misc_fns import spatial_ave, calcsatspechum
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, shiftgrid, maskoceans
import cdutil
import matplotlib

matplotlib.rcParams.update({'font.size': 20})
matplotlib.rcParams.update({'figure.figsize': (8,6)})
matplotlib.rcParams.update({'lines.linewidth': 3})
matplotlib.rcParams.update({'legend.fontsize': 20})
matplotlib.rcParams.update({'mathtext.fontset': 'cm'})

fdata = '/Users/cpatrizio/data/MERRA2/'
fout = '/Volumes/GoogleDrive/My Drive/PhD/figures/NAindex/'
fkernels = '/Volumes/GoogleDrive/My Drive/data/cam5-kernels/kernels/'

############## LOAD KERNELS: CAM5 KERNELS FOR MERRA2
print 'loading kernels...'
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
 
############## LOAD DATA
print 'loading data...'
#MERRA-2 DATA, t, qv, ts, and radiative fields
fnamet = glob.glob(fdata + 'MERRA2_t_monthly1980to2017.nc')[0]
fnameqv = glob.glob(fdata + 'MERRA2_qv_monthly1980to2017.nc')[0]
fnametsps = glob.glob(fdata + 'MERRA2_tsps_monthly1980to2017.nc')[0]
fnamerad = glob.glob(fdata + 'MERRA2_rad_monthly1980to2017.nc')[0]

tfile = cdms2.open(fnamet)
qvfile = cdms2.open(fnameqv)
tspsfile = cdms2.open(fnametsps)
radfile = cdms2.open(fnamerad)

tinit = 12*3


#cdutil.setTimeBoundsMonthly(swdwns)
#cdutil.setTimeBoundsMonthly(swnets_cs)
#cdutil.setTimeBoundsMonthly(swdwns_cs)

#LOAD t, ts, qv, AND albedo TO DO CLEAR-SKY ENERGY CLOSURE OF SEASONAL CYCLE
#LOAD all-sky and clear-sky LW, SW radiation TO DO ALL-SKY ENERGY CLOSURE OF SEASONAL CYCLE
ta = tfile['T']
ta = ta[-tinit:-1,:]
cdutil.setTimeBoundsMonthly(ta)
ta= cdutil.ANNUALCYCLE.climatology(ta)

ts = tspsfile['TS']
ts = ts[-tinit:-1,:]
cdutil.setTimeBoundsMonthly(ts)
ts =  cdutil.ANNUALCYCLE.climatology(ts)

qv = qvfile['QV']
qv = qv[-tinit:-1,:]
cdutil.setTimeBoundsMonthly(qv)
qv = cdutil.ANNUALCYCLE.climatology(qv)


#LW_net_surf = radfile['LWGNT']
#LW_net_surf = LW_net_surf[-tinit:-1,:]
#cdutil.setTimeBoundsMonthly(LW_net_surf)
#LW_net_surf = cdutil.ANNUALCYCLE.climatology(LW_net_surf)
LW_net_surf_cs = radfile('LWGNTCLR')
LW_net_surf_cs = LW_net_surf_cs[-tinit:-1,:]
cdutil.setTimeBoundsMonthly(LW_net_surf_cs)
LW_net_surf_cs = cdutil.ANNUALCYCLE.climatology(LW_net_surf_cs) #this is net radiation at surface (positive down)
#LW_net_TOA = radfile['LWTUP']
#LW_net_TOA = LW_net_TOA[-tinit:-1,:]
#cdutil.setTimeBoundsMonthly(LW_net_TOA)
#LW_net_TOA = cdutil.ANNUALCYCLE.climatology(LW_net_TOA)
LW_net_TOA_cs = radfile('LWTUPCLR')
LW_net_TOA_cs = LW_net_TOA_cs[-tinit:-1,:]
cdutil.setTimeBoundsMonthly(LW_net_TOA_cs)
LW_net_TOA_cs = cdutil.ANNUALCYCLE.climatology(LW_net_TOA_cs) #this is upwelling radiation (positive up)

#SW_net_surf = radfile['SWGNT']
#SW_net_surf = SW_net_surf[-tinit:-1,:]
#swnets_cs = radfile('SWGNTCLR')
#swnettoa = radfile('SWTNT')
#swnettoa_cs = radfile('SWTNTCLR')
#SW_dwn_surf = radfile['SWGDN']
#SW_dwn_surf = SW_dwn_surf[-tinit:-1,:]
#swdwns_cs = radfile('SWGDNCLR')
#psvar = tspsfile['PS']

#alb = 1 - SW_net_surf/SW_dwn_surf
#cdutil.setTimeBoundsMonthly(alb)
#alb = cdutil.ANNUALCYCLE.departures(alb)


merralats = ts.getLatitude()[:]
merralons = ts.getLongitude()[:]

#TODO: Get clear-sky radiative variables

#TODO: Get cloud fraction (low, mid, high)

#TODO: Mask land from variables?

#HORIZONTALLY INTERPOLATE TO KERNEL GRID
print 'interpolating to kernel grid...'
ta = ta.regrid(grid, regridTool="esmf", regridMethod = "linear")
ts = ts.regrid(grid, regridTool="esmf", regridMethod = "linear")
qv = qv.regrid(grid, regridTool="esmf", regridMethod = "linear")
LW_net_TOA_cs = LW_net_TOA_cs.regrid(grid, regridTool="esmf", regridMethod = "linear")
LW_net_surf_cs = LW_net_surf_cs.regrid(grid, regridTool="esmf", regridMethod = "linear")
#alb = alb.regrid(grid, regridTool="esmf", regridMethod = "linear")


#COMPUTE MONTHLY ANOMALIES (FROM ANNUAL MEAN)
tabar = np.ma.average(ta, axis=0)
tabar_grid = np.zeros(ta.shape)
tabar_grid[:] = tabar
tsbar = np.ma.average(ts, axis=0)
tsbar_grid = np.zeros(ts.shape)
tsbar_grid[:] = tsbar
qvbar = np.ma.average(qv, axis=0)
qvbar_grid = np.zeros(qv.shape)
qvbar_grid[:] = qvbar
#albbar = np.ma.average(alb,axis=0)
#albbar_grid = np.zeros(alb.shape)

ta_anom = ta - tabar_grid
ts_anom = ts - tsbar_grid
qv_anom = qv - qvbar_grid
#alb_anom = alb - albbar_grid

LW_net_TOA_cs_bar = np.ma.average(LW_net_TOA_cs, axis=0)
LW_net_TOA_cs_bargrid = np.zeros(LW_net_TOA_cs.shape)
LW_net_TOA_cs_bargrid[:] = LW_net_TOA_cs_bar
LW_net_surf_cs_bar = np.ma.average(LW_net_surf_cs, axis=0)
LW_net_surf_cs_bargrid = np.zeros(LW_net_surf_cs.shape)
LW_net_surf_cs_bargrid[:] = LW_net_surf_cs_bar

LW_net_TOA_cs_anom = LW_net_TOA_cs - LW_net_TOA_cs_bargrid
LW_net_surf_cs_anom = LW_net_surf_cs - LW_net_surf_cs_bargrid

############## DECOMPOSE ALL-SKY AND CLEAR-SKY RADIATION INTO DIFFERENT COMPONENTS (t, qv, albedo, and cloud)
print 'decomposing radiation fields...'
print 'temperature component...'
#TEMPERATURE COMPONENT
dts=ts_anom
#dts_globalmean = np.ma.average(spatial_ave(dts,lats))

#ts_kernel_TOA = ftskern('FLNT')
#ts_kernel_surf = ftskern('FLNS')
ts_kernel_TOA_cs = ftskern('FLNTC')
ts_kernel_surf_cs = ftskern('FLNSC')

#dLW_ts_TOA=ts_kernel_TOA*dts  
#dLW_ts_surf = ts_kernel_surf*dts
dLW_ts_TOA_cs = ts_kernel_TOA_cs*dts
dLW_ts_surf_cs = ts_kernel_surf_cs*dts

dta = ta_anom
#ta_kernel_TOA = ftkern('FLNT')
#ta_kernel_surf = ftkern('FLNS')
ta_kernel_TOA_cs = ftkern['FLNTC']
ta_kernel_surf_cs = ftkern['FLNSC']

mask_temp = dta.mask
#mask stratospheric values
dta = np.ma.MaskedArray(dta, mask=np.bitwise_and(mask_temp, p>=p_tropopause))

#dLW_ta_TOA = np.ma.sum(ta_kernel_TOA*dta*pdiff, axis=1)
#dLW_ta_surf = np.ma.sum(ta_kernel_surf*dta*pdiff, axis=1)
dLW_ta_TOA_cs = np.ma.sum(ta_kernel_TOA_cs*dta*pdiff, axis=1)
dLW_ta_surf_cs = np.ma.sum(ta_kernel_surf_cs*dta*pdiff, axis=1)


#WATER VAPOR COMPONENT
print 'water vapor component...'
q_LW_kernel_TOA_cs=fqkern['FLNTC']
q_LW_kernel_surf_cs=fqkern['FLNSC']

dq=qv_anom

q1 = qv
t1 = ta
qs1 = calcsatspechum(t1,p)
qs2 = calcsatspechum(t1+dta,p)
dqsdt = (qs2 - qs1)/dta
rh = q1/qs1
dqdt = rh*dqsdt

q_LW_kernel_TOA_cs=q_LW_kernel_TOA_cs/dqdt
q_LW_kernel_surf_cs=q_LW_kernel_surf_cs/dqdt

mask_temp = dq.mask
#mask stratospheric values
dq = np.ma.MaskedArray(dq, mask=np.bitwise_and(mask_temp, p>=p_tropopause))
#
#% Convolve moisture kernel with change in moisture
dLW_q_TOA_cs=np.ma.sum(q_LW_kernel_TOA_cs*dq*pdiff, axis=1)
dLW_q_surf_cs=np.ma.sum(q_LW_kernel_surf_cs*dq*pdiff, axis=1)


#TODO: Calculate surface albedo component of TOA, ATM, and surface radiation 

#TODO: Compute cloud component of TOA, ATM and surface radiation

############## SELECT REGION FOR ANALYSIS
nlats = len(lats)
lindx = nlats/2
uindx = np.where(lats > 60)[0][0]
reglats = lats[lindx:uindx]

#dLW_t_TOA_globalmean = spatial_ave(-dLW_ta_TOA-dLW_ts_TOA, lats)
#dLW_t_surf_globalmean = spatial_ave(-dLW_ta_surf-dLW_ts_surf, lats)

dLW_t_TOA_cs_globalmean = spatial_ave(-dLW_ta_TOA_cs-dLW_ts_TOA_cs, lats)
dLW_t_surf_cs_globalmean = spatial_ave(-dLW_ta_surf_cs-dLW_ts_surf_cs, lats)

dLW_q_TOA_cs_globalmean = spatial_ave(-dLW_q_TOA_cs, lats)
dLW_q_surf_cs_globalmean = spatial_ave(-dLW_q_surf_cs, lats)

#dLW_t_TOA_NHmean = spatial_ave(-dLW_ta_TOA[:,lindx:uindx,:]-dLW_ts_TOA[:,lindx:uindx,:], reglats)
#dLW_t_surf_NHmean = spatial_ave(-dLW_ta_surf[:,lindx:uindx,:]-dLW_ts_surf[:,lindx:uindx,:], reglats)

dLW_t_TOA_cs_NHmean = spatial_ave(-dLW_ta_TOA_cs[:,lindx:uindx,:]-dLW_ts_TOA_cs[:,lindx:uindx,:], reglats)
dLW_t_surf_cs_NHmean = spatial_ave(-dLW_ta_surf_cs[:,lindx:uindx,:]-dLW_ts_surf_cs[:,lindx:uindx,:], reglats)

dLW_q_TOA_cs_NHmean = spatial_ave(-dLW_q_TOA_cs[:,lindx:uindx,:], reglats)
dLW_q_surf_cs_NHmean = spatial_ave(-dLW_q_surf_cs[:,lindx:uindx,:], reglats)

#calculate clear-sky longwave net radiation using kernel method
#sum of clear-sky longwave temperature component and water vapor component is the net clear-sky longwave radiation 
dLW_net_TOA_cs_globalmean = dLW_t_TOA_cs_globalmean + dLW_q_TOA_cs_globalmean
dLW_net_TOA_cs_NHmean = dLW_t_TOA_cs_NHmean + dLW_q_TOA_cs_NHmean

dLW_net_surf_cs_globalmean = dLW_t_surf_cs_globalmean + dLW_q_surf_cs_globalmean
dLW_net_surf_cs_NHmean = dLW_t_surf_cs_NHmean + dLW_q_surf_cs_NHmean

#calculate clear-sky longwave net radiation from MERRA2
LW_net_TOA_cs_anom_globalmean = spatial_ave(-LW_net_TOA_cs_anom, lats)
LW_net_surf_cs_anom_globalmean = spatial_ave(LW_net_surf_cs_anom, lats)

LW_net_TOA_cs_anom_NHmean = spatial_ave(-LW_net_TOA_cs_anom[:,lindx:uindx,:], reglats)
LW_net_surf_cs_anom_NHmean = spatial_ave(LW_net_surf_cs_anom[:,lindx:uindx,:], reglats)

#TODO: Calculate seasonal cycle of cloud types (low, mid, high)

#TODO: Compute feedbacks for each radiation component (normalize by SST in NH?)

############## PLOTTING
print 'plotting...'
dts_NHmean = spatial_ave(dts[:,lindx:uindx,:], reglats)
p_s = 1000*1e2
dta_vertmean = (1./p_s)*np.ma.sum(dta*pdiff, axis=1)
dq_vertmean = (1./p_s)*np.ma.sum(dq*pdiff, axis=1)
dq_NHmean = spatial_ave(dq_vertmean[:,lindx:uindx,:], reglats)
dta_NHmean = spatial_ave(dta_vertmean[:,lindx:uindx,:], reglats)

#Plot NH averaged t and ts monthly anomalies
f, axarr = plt.subplots(3,1, figsize=(16,10))
axarr[0,].plot(dts_NHmean)
axarr[0,].axhline(0, color='k', alpha=0.5, linewidth=1)
axarr[0,].set_title(r'$T_{s}$')
axarr[0,].set_ylabel('K')
axarr[1,].plot(dta_NHmean)
axarr[1,].axhline(0, color='k', alpha=0.5, linewidth=1)
axarr[1,].set_title(r'$T_a$')
axarr[1,].set_ylabel('K')
axarr[1,].set_xlabel('time (month)')
axarr[2,].plot(dq_NHmean)
axarr[2,].axhline(0, color='k', alpha=0.5, linewidth=1)
axarr[2,].set_xlabel('time (month)')
axarr[2,].set_ylabel('kg/kg')
plt.savefig(fout + 'MERRA2_tqv_seasonalcycle.pdf')
plt.close()


#Plot temperature component of TOA, ATM, and surface all-sky radiation
#f, axarr = plt.subplots(3,1, figsize=(14,16))
#axarr[0,].plot(dLW_t_TOA_NHmean, color='r', label=r'$\Delta LW^{all}_{T}$')
#axarr[0,].axhline(0, color='k', alpha=0.5, linewidth=1)
#axarr[0,].set_title('TOA')
#axarr[0,].set_ylabel('W m$^{-2}$')
#axarr[1,].plot(dLW_t_surf_NHmean, color='r',)
#axarr[1,].axhline(0, color='k', alpha=0.5, linewidth=1)
#axarr[1,].set_title('SURF')
#axarr[1,].set_ylabel('W m$^{-2}$')
#axarr[2,].plot(dLW_t_surf_NHmean-dLW_t_TOA_NHmean,  color='r')
#axarr[2,].axhline(0, color='k', alpha=0.5, linewidth=1)
#axarr[2,].set_title('ATM')
#axarr[2,].set_xlabel('time (month)')
#axarr[2,].set_ylabel('W m$^{-2}$')
#axarr[0,].legend()
#plt.savefig(fout + 'MERRA2_deltaLW_t_monthly.pdf')
#plt.close()

minLW = -40
maxLW = 40

#Plot temperature component of TOA, ATM, and surface clear-sky radiation
f, axarr = plt.subplots(3,1, figsize=(16,18))
axarr[0,].plot(dLW_t_TOA_cs_NHmean,  color='r', label=r'$\Delta LW^{clear}_{T}$')
axarr[0,].plot(dLW_q_TOA_cs_NHmean,  color='g', label=r'$\Delta LW^{clear}_{q}$')
axarr[0,].plot(dLW_net_TOA_cs_NHmean,  color='grey', label=r'$\Delta LW^{clear}_{T} + \Delta LW^{clear}_{q}$')
axarr[0,].plot(LW_net_TOA_cs_anom_NHmean,  color='k', label=r'$\Delta LW^{clear}$')
axarr[0,].axhline(0, color='k', alpha=0.5, linewidth=1)
axarr[0,].set_title(r'TOA')
axarr[0,].set_ylabel('W m$^{-2}$')
axarr[0,].set_ylim(minLW,maxLW)
axarr[1,].plot(dLW_t_surf_cs_NHmean,  color='r', label=r'$\Delta LW^{clear}_{T}$')
axarr[1,].plot(dLW_q_surf_cs_NHmean,  color='g', label=r'$\Delta LW^{clear}_{q}$')
axarr[1,].plot(dLW_net_surf_cs_NHmean,  color='grey', label=r'$\Delta LW^{clear}_{T} + \Delta LW^{clear}_{q}$')
axarr[1,].plot(LW_net_surf_cs_anom_NHmean,  color='k', label=r'$\Delta LW^{clear}$')
axarr[1,].axhline(0, color='k', alpha=0.5, linewidth=1)
axarr[1,].set_title(r'SURF')
axarr[1,].set_ylabel('W m$^{-2}$')
axarr[1,].set_ylim(minLW,maxLW)
axarr[2,].plot(dLW_t_TOA_cs_NHmean-dLW_t_surf_cs_NHmean, color='r', label=r'$\Delta LW^{clear}_{T}$')
axarr[2,].plot(dLW_t_TOA_cs_NHmean-dLW_q_surf_cs_NHmean,  color='g', label=r'$\Delta LW^{clear}_{q}$')
axarr[2,].plot(dLW_net_TOA_cs_NHmean-dLW_net_surf_cs_NHmean,  color='grey', label=r'$\Delta LW^{clear}_{T} + \Delta LW^{clear}_{q}$')
axarr[2,].plot(LW_net_TOA_cs_anom_NHmean-LW_net_surf_cs_anom_NHmean,  color='k', label=r'$\Delta LW^{clear}$')
axarr[2,].axhline(0, color='k', alpha=0.5, linewidth=1)
axarr[2,].set_title(r'ATM')
axarr[2,].set_xlabel('time (month)')
axarr[2,].set_ylabel('W m$^{-2}$')
axarr[2,].set_ylim(minLW,maxLW)
axarr[2,].legend()
plt.savefig(fout + 'MERRA2_deltaLWclear_seasonalcycle.pdf')
plt.close()


#TODO: Plot albedo component of clear-sky radiation

#TODO: Plot water vapor component of clear-sky radiation

#TODO: Plot clear-sky radiation from MERRA-2 on plots of decomposed clear-sky radiation

#TODO: Plot cloud component of all-sky radiation

#TODO: Plot water vapor component of all-sky radiation

#TODO: Plot all-sky radiation from MERRA-2 on plots of decomposed all-sky radiation

#TODO: Plot maps of each radiation component
































