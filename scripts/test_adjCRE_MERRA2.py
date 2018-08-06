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
import MV2 as MV
from thermolib.thermo import r_star

matplotlib.rcParams.update({'font.size': 20})
matplotlib.rcParams.update({'figure.figsize': (8,6)})
matplotlib.rcParams.update({'lines.linewidth': 3})
matplotlib.rcParams.update({'legend.fontsize': 20})
matplotlib.rcParams.update({'mathtext.fontset': 'cm'})

fdata = '/Users/cpatrizio/data/MERRA2/'
fout = '/Volumes/GoogleDrive/My Drive/PhD/figures/NA index/'
fkernels = '/Volumes/GoogleDrive/My Drive/data/cam5-kernels/kernels/'

############## LOAD KERNELS: CAM5 KERNELS FOR MERRA2
print 'loading kernels...'
falbkern = cdms2.open(fkernels + 'alb.kernel.nc')
fqkern = cdms2.open(fkernels + 'q.kernel.MERRA2plev.nc')
ftkern = cdms2.open(fkernels + 't.kernel.MERRA2plev.nc')
#ftkerngw = cdms2.open(fin + 't.kernel.nc')
ftskern = cdms2.open(fkernels + 'ts.kernel.nc')



pdiff = cdms2.open(fdata + 'dp_MERRA2plev.nc')['dp'][:]/100.
pdiff = MV.average(pdiff,axis=0)

#NOTE: plev kernels have units W m^-2 hPa^-1

FLNSt = ftkern('FLNS')
#gw = ftkerngw('gw')

#lats = FLNSt.getLatitude()[:]
#lons = FLNSt.getLongitude()[:]
#t = FLNSt.getTime()[:]
#nlat = lats.size
#nlon = lons.size
#nt = t.size
p_Pa = ftkern('plev')
nplevs = p_Pa.size

#p_tropopause = np.repeat(p_tropopause[np.newaxis,...],nt,axis=0)

#grid = cdms2.createGenericGrid(lats,lons)
 
############## LOAD DATA
print 'loading data...'
#MERRA-2 DATA, t, qv, ts, and radiative fields
#fdeltat = glob.glob(fdata + 'deltat3D_SST4to20_detrT_2x2.npy')[0]
#fdeltaqv = glob.glob(fdata + 'deltaqv3D_SST4to20_detrT_2x2.npy')[0]
#fdeltats = glob.glob(fdata + 'deltaSST_SST4to20_detrT_2x2.npz')[0]
#
##need to load all radiation fields here (LW, SW and clear-sky fields)
#fdeltaSWcs = glob.glob(fdata + 'deltaSWnetsurfcs_SST4to20_detrT_2x2.npy')[0]
#fdeltaLWcs = glob.glob(fdata + 'deltaLWnetsurfcs_SST4to20_detrT_2x2.npy')[0]
#
#
#fdeltaSW = glob.glob(fdata + 'deltaSWnetsurf_SST4to20_detrT_2x2.npy')[0]
#fdeltaLW = glob.glob(fdata + 'deltaLWnetsurf_SST4to20_detrT_2x2.npy')[0]
#
#fdeltaalb = glob.glob(fdata + 'deltaalbedo_SST4to20_detrT_2x2.npy')[0]
#
#fdeltaCREsurf = glob.glob(fdata + 'deltaCREsurf_SST4to20_detrT_2x2.npy')[0]
#
#fdeltaSWCREsurf = glob.glob(fdata + 'deltaSWCREsurf_SST4to20_detrT_2x2.npy')[0]
#
#fdeltaLWCREsurf = glob.glob(fdata + 'deltaLWCREsurf_SST4to20_detrT_2x2.npy')[0]

fdeltat = glob.glob(fdata + 'deltat3D_SST4to20_7.0LP_detrT_2x2.npy')[0]
fdeltaqv = glob.glob(fdata + 'deltaqv3D_SST4to20_7.0LP_detrT_2x2.npy')[0]
fdeltats = glob.glob(fdata + 'deltaSST_SST4to20_7.0LP_detrT_2x2.npz')[0]

#need to load all radiation fields here (LW, SW and clear-sky fields)
fdeltaSWcs = glob.glob(fdata + 'deltaSWnetsurfcs_SST4to20_7.0LP_detrT_2x2.npy')[0]
fdeltaLWcs = glob.glob(fdata + 'deltaLWnetsurfcs_SST4to20_7.0LP_detrT_2x2.npy')[0]
fdeltaalb = glob.glob(fdata + 'deltaalbedo_SST4to20_7.0LP_detrT_2x2.npy')[0]

fdeltaSW = glob.glob(fdata + 'deltaSWnetsurf_SST4to20_7.0LP_detrT_2x2.npy')[0]
fdeltaLW = glob.glob(fdata + 'deltaLWnetsurf_SST4to20_7.0LP_detrT_2x2.npy')[0]

fdeltaCREsurf = glob.glob(fdata + 'deltaCREsurf_SST4to20_7.0LP_detrT_2x2.npy')[0]

fdeltaSWCREsurf = glob.glob(fdata + 'deltaSWCREsurf_SST4to20_7.0LP_detrT_2x2.npy')[0]

fdeltaLWCREsurf = glob.glob(fdata + 'deltaLWCREsurf_SST4to20_7.0LP_detrT_2x2.npy')[0]

fqv3D = glob.glob(fdata + 'MERRA2_qv3D_monthly1980to2017.nc')[0]
fT3D = glob.glob(fdata + 'MERRA2_t3D_monthly1980to2017.nc')[0]

qv3D = cdms2.open(fqv3D)['QV']
T3D = cdms2.open(fT3D)['T']

qv = MV.average(qv3D, axis=0) 
ta = MV.average(T3D, axis=0)


ta_anom = np.load(fdeltat)
qv_anom = np.load(fdeltaqv)*1e-3 #convert to kg/kg
ts_anom = np.load(fdeltats)['sst']
lats = np.load(fdeltats)['lats']
lons = np.load(fdeltats)['lons']
SW_net_surf_cs_anom = np.load(fdeltaSWcs)
LW_net_surf_cs_anom = np.load(fdeltaLWcs)
SW_net_surf_anom = np.load(fdeltaSW)
LW_net_surf_anom = np.load(fdeltaLW)
CRE_net_surf_anom = np.load(fdeltaCREsurf)
SWCRE_net_surf_anom = np.load(fdeltaSWCREsurf)
LWCRE_net_surf_anom = np.load(fdeltaLWCREsurf)
alb_anom = np.load(fdeltaalb)


ta_anom = MV.array(ta_anom, mask=np.isnan(ta_anom))
ts_anom = MV.array(ts_anom, mask=np.isnan(ts_anom))
ts_anom = MV.array(ts_anom, mask=np.abs(ts_anom) > 1e4)
qv_anom = MV.array(qv_anom, mask=np.isnan(qv_anom))
SW_net_surf_cs_anom = MV.array(SW_net_surf_cs_anom, mask=np.isnan(SW_net_surf_cs_anom))
LW_net_surf_cs_anom = MV.array(LW_net_surf_cs_anom, mask=np.isnan(LW_net_surf_cs_anom))
SW_net_surf_anom = MV.array(SW_net_surf_anom, mask=np.isnan(SW_net_surf_anom))
LW_net_surf_anom = MV.array(LW_net_surf_anom, mask=np.isnan(LW_net_surf_anom))
CRE_net_surf_anom = MV.array(CRE_net_surf_anom, mask=np.isnan(CRE_net_surf_anom))
SWCRE_net_surf_anom = MV.array(SWCRE_net_surf_anom, mask=np.isnan(SWCRE_net_surf_anom))
LWCRE_net_surf_anom = MV.array(LWCRE_net_surf_anom, mask=np.isnan(LWCRE_net_surf_anom))
alb_anom = MV.array(alb_anom, mask=np.isnan(alb_anom))


lonbounds = [lons[0],lons[-1]]
latbounds = [lats[0],lats[-1]]



nlat = len(lats)
nlon = len(lons)
#londiff=2
#latdiff=2
#
#lats = np.arange(0, nlat*latdiff,latdiff)
#lons = np.arange(0, nlon*londiff,londiff)

grid = cdms2.createGenericGrid(lats,lons)

p = np.zeros((nplevs, nlat, nlon)).T
p[:,:,:] = p_Pa/100.
p = p.T

x=np.cos(np.deg2rad(lats))
p_tropopause_zonalmean=300-200*x

p_tropopause = np.zeros((nlat, nlon)).T
p_tropopause[:,:] = p_tropopause_zonalmean
p_tropopause = p_tropopause.T
p_tropopause = np.repeat(p_tropopause[np.newaxis,...],nplevs,axis=0)



#LOAD KERNELS AND
#HORIZONTALLY INTERPOLATE KERNEL GRID TO MERRA2 GRID
#ts_kernel_TOA_cs = ftskern('FLNTC')
ts_kernel_surf_cs = ftskern('FLNSC')

#ta_kernel_TOA_cs = ftkern['FLNTC']
ta_kernel_surf_cs = ftkern('FLNSC')

#q_LW_kernel_TOA_cs=fqkern['FLNTC']
q_LW_kernel_surf_cs=fqkern('FLNSC')

q_SW_kernel_surf_cs = fqkern('FSNSC')

alb_kernel_surf_cs = falbkern('FSNSC')


#ts_kernel_TOA_cs = ftskern('FLNTC')
ts_kernel_surf = ftskern('FLNS')

#ta_kernel_TOA_cs = ftkern['FLNTC']
ta_kernel_surf = ftkern('FLNS')

#q_LW_kernel_TOA_cs=fqkern['FLNTC']
q_LW_kernel_surf =fqkern('FLNS')

q_SW_kernel_surf = fqkern('FSNS')

alb_kernel_surf = falbkern('FSNS')

#convert longwave kernels to positive down
ts_kernel_surf_cs = -ts_kernel_surf_cs 
ta_kernel_surf_cs = -ta_kernel_surf_cs 
ts_kernel_surf = -ts_kernel_surf 
ta_kernel_surf = -ta_kernel_surf

q_LW_kernel_surf_cs = -q_LW_kernel_surf_cs
q_LW_kernel_surf = -q_LW_kernel_surf



ts_kernel_surf_cs = MV.average(ts_kernel_surf_cs, axis=0)
ta_kernel_surf_cs = MV.average(ta_kernel_surf_cs, axis=0)
q_LW_kernel_surf_cs = MV.average(q_LW_kernel_surf_cs, axis=0)
q_SW_kernel_surf_cs = MV.average(q_SW_kernel_surf_cs, axis=0)
alb_kernel_surf_cs = MV.average(alb_kernel_surf_cs, axis=0)

ts_kernel_surf = MV.average(ts_kernel_surf, axis=0)
ta_kernel_surf = MV.average(ta_kernel_surf, axis=0)
q_LW_kernel_surf = MV.average(q_LW_kernel_surf, axis=0)
q_SW_kernel_surf = MV.average(q_SW_kernel_surf, axis=0)
alb_kernel_surf = MV.average(alb_kernel_surf, axis=0)


#ts_kernel_surf_cs = ts_kernel_surf_cs.subRegion(latitude=(lats[0], lats[-1]),longitude=(lons[0], lons[-1]))
#ta_kernel_surf_cs = ta_kernel_surf_cs.subRegion(latitude=(lats[0], lats[-1]),longitude=(lons[0], lons[-1]))
#q_LW_kernel_surf_cs = q_LW_kernel_surf_cs.subRegion(latitude=(lats[0], lats[-1]),longitude=(lons[0], lons[-1]))
#q_SW_kernel_surf_cs = q_SW_kernel_surf_cs.subRegion(latitude=(lats[0], lats[-1]),longitude=(lons[0], lons[-1]))
#alb_kernel_surf_cs = alb_kernel_surf_cs.subRegion(latitude=(lats[0], lats[-1]),longitude=(lons[0], lons[-1]))
#pdiff = pdiff.subRegion(latitude=(lats[0], lats[-1]),longitude=(lons[0], lons[-1]))
#ta = ta.subRegion(latitude=(lats[0], lats[-1]),longitude=(lons[0], lons[-1]))
#qv = qv.subRegion(latitude=(lats[0], lats[-1]),longitude=(lons[0], lons[-1]))



print 'interpolating to MERRA2 grid...'
#ta = ta.regrid(grid, regridTool="esmf", regridMethod = "linear")
#ts = ts.regrid(grid, regridTool="esmf", regridMethod = "linear")
#qv = qv.regrid(grid, regridTool="esmf", regridMethod = "linear")
#LW_net_TOA_cs = LW_net_TOA_cs.regrid(grid, regridTool="esmf", regridMethod = "linear")
#LW_net_surf_cs = LW_net_surf_cs.regrid(grid, regridTool="esmf", regridMethod = "linear")
#alb = alb.regrid(grid, regridTool="esmf", regridMethod = "linear")

#ts_kernel_TOA_cs = ts_kernel_TOA_cs.regrid(grid, regridTool="esmf", regridMethod = "linear")
ts_kernel_surf_cs = ts_kernel_surf_cs.regrid(grid, regridTool="esmf", regridMethod = "linear")
ta_kernel_surf_cs = ta_kernel_surf_cs.regrid(grid, regridTool="esmf", regridMethod = "linear")
q_LW_kernel_surf_cs = q_LW_kernel_surf_cs.regrid(grid, regridTool="esmf", regridMethod = "linear")
q_SW_kernel_surf_cs = q_SW_kernel_surf_cs.regrid(grid, regridTool="esmf", regridMethod = "linear")
alb_kernel_surf_cs = alb_kernel_surf_cs.regrid(grid, regridTool="esmf", regridMethod = "linear")
pdiff = pdiff.regrid(grid, regridTool="esmf", regridMethod = "linear")
ta = ta.regrid(grid, regridTool="esmf", regridMethod = "linear")
qv = qv.regrid(grid, regridTool="esmf", regridMethod = "linear")

ts_kernel_surf = ts_kernel_surf.regrid(grid, regridTool="esmf", regridMethod = "linear")
ta_kernel_surf = ta_kernel_surf.regrid(grid, regridTool="esmf", regridMethod = "linear")
q_LW_kernel_surf = q_LW_kernel_surf.regrid(grid, regridTool="esmf", regridMethod = "linear")
q_SW_kernel_surf = q_SW_kernel_surf.regrid(grid, regridTool="esmf", regridMethod = "linear")
alb_kernel_surf = alb_kernel_surf.regrid(grid, regridTool="esmf", regridMethod = "linear")





############## DECOMPOSE ALL-SKY AND CLEAR-SKY RADIATION INTO DIFFERENT COMPONENTS (t, qv, albedo, and cloud)
print 'decomposing radiation fields...'
print 'temperature component...'
#TEMPERATURE COMPONENT
dts=ts_anom
#dts_globalmean = np.ma.average(spatial_ave(dts,lats))



#dLW_ts_TOA=ts_kernel_TOA*dts  
#dLW_ts_surf = ts_kernel_surf*dts
#dLW_ts_TOA_cs = ts_kernel_TOA_cs*dts
dLW_ts_surf_cs = ts_kernel_surf_cs*dts

dLW_ts_surf = ts_kernel_surf*dts

dta = ta_anom


mask_temp = dta.mask
#TODO: MASK STRATOSPHERE?
#dta = np.ma.MaskedArray(dta, mask=np.bitwise_and(mask_temp, p>=p_tropopause))

#dLW_ta_TOA = np.ma.sum(ta_kernel_TOA*dta*pdiff, axis=1)
#dLW_ta_surf = np.ma.sum(ta_kernel_surf*dta*pdiff, axis=1)
#dLW_ta_TOA_cs = np.ma.sum(ta_kernel_TOA_cs*dta*pdiff, axis=0)
dLW_ta_surf_cs = np.ma.sum(ta_kernel_surf_cs*dta*pdiff, axis=0)

dLW_ta_surf = np.ma.sum(ta_kernel_surf*dta*pdiff, axis=0)


#WATER VAPOR COMPONENT
print 'water vapor component...'

dq=qv_anom

#TODO: fix q kernel 
q1 = qv
t1 = ta
#qs1 = r_star(p*1e2,t1)
qs1 = calcsatspechum(t1, p)
qs2 = calcsatspechum(t1+1, p)
#qs1 = np.ma.array(qs1, mask=qs1<0)
#qs2 = r_star(p*1e2,t1+dta)
#qs2 = np.ma.array(qs2, mask=qs2<0)
#dqsdt = (qs2 - qs1)/dta
dqsdt = (qs2 - qs1)
rh = q1/qs1
dqdt = rh*dqsdt

#q_LW_kernel_TOA_cs=q_LW_kernel_TOA_cs/dqdt
q_LW_kernel_surf_cs=q_LW_kernel_surf_cs/dqdt

q_SW_kernel_surf_cs=q_SW_kernel_surf_cs/dqdt

q_LW_kernel_surf=q_LW_kernel_surf/dqdt

q_SW_kernel_surf=q_SW_kernel_surf/dqdt

mask_temp = dq.mask
#TODO: mask stratospheric values
#dq = np.ma.MaskedArray(dq, mask=np.bitwise_and(mask_temp, p>=p_tropopause))
#
#% Convolve moisture kernel with change in moisture
#dLW_q_TOA_cs=np.ma.sum(q_LW_kernel_TOA_cs*dq*pdiff, axis=1)
dLW_q_surf_cs=np.ma.sum(q_LW_kernel_surf_cs*dq*pdiff, axis=0)

dSW_q_surf_cs=np.ma.sum(q_SW_kernel_surf_cs*dq*pdiff, axis=0)

dLW_q_surf=np.ma.sum(q_LW_kernel_surf*dq*pdiff, axis=0)

dSW_q_surf=np.ma.sum(q_SW_kernel_surf*dq*pdiff, axis=0)

print 'albedo component'

dalb = alb_anom

dSW_alb_surf_cs = alb_kernel_surf_cs*dalb

dSW_alb_surf = alb_kernel_surf*dalb


#TODO: Compute cloud component of surface radiation

dLW_t_surf = dLW_ta_surf + dLW_ts_surf
dLW_t_surf_cs = dLW_ta_surf_cs + dLW_ts_surf_cs

dLWCRE_adj = LWCRE_net_surf_anom - (dLW_t_surf - dLW_t_surf_cs) - (dLW_q_surf - dLW_q_surf_cs)
dSWCRE_adj = SWCRE_net_surf_anom - (dSW_q_surf - dSW_q_surf_cs) - (dSW_alb_surf - dSW_alb_surf_cs)


############## SELECT REGION FOR ANALYSIS
lonbounds = [305.,335.]
#latbounds = [0,70.]
#
#si = np.argmin(np.abs(lats - latbounds[0]))
#ni = np.argmin(np.abs(lats - latbounds[1]))
wi = np.argmin(np.abs(lons - lonbounds[0]))
ei = np.argmin(np.abs(lons - lonbounds[1]))

#calculate clear-sky temperature, moisture, and albedo components of net surface radiation
dLW_t_surf_cs_zmean = MV.average(dLW_ta_surf_cs[:,wi:ei]+dLW_ts_surf_cs[:,wi:ei],axis=1)

dLW_q_surf_cs_zmean = MV.average(dLW_q_surf_cs[:,wi:ei],axis=1)

dSW_q_surf_cs_zmean = MV.average(dSW_q_surf_cs[:,wi:ei],axis=1)

dSW_alb_cs_zmean = MV.average(dSW_alb_surf_cs[:,wi:ei],axis=1)
 
dLW_net_surf_cs_zmean = dLW_t_surf_cs_zmean + dLW_q_surf_cs_zmean

#calculate clear-sky net radiation from MERRA2

LW_net_surf_cs_anom_zmean = np.ma.average(LW_net_surf_cs_anom[:,wi:ei],axis=1)

SW_net_surf_cs_anom_zmean = np.ma.average(SW_net_surf_cs_anom[:,wi:ei],axis=1)

dSW_net_surf_cs_zmean =  dSW_q_surf_cs_zmean + dSW_alb_cs_zmean

dLWCRE_adj_zmean = np.ma.average(dLWCRE_adj[:,wi:ei], axis=1)
dSWCRE_adj_zmean = np.ma.average(dSWCRE_adj[:,wi:ei], axis=1)



dLW_t_surf_zmean = MV.average(dLW_ta_surf[:,wi:ei]+dLW_ts_surf[:,wi:ei],axis=1)

dLW_q_surf_zmean = MV.average(dLW_q_surf[:,wi:ei],axis=1)

dSW_q_surf_zmean = MV.average(dSW_q_surf[:,wi:ei],axis=1)

dSW_alb_zmean = MV.average(dSW_alb_surf[:,wi:ei],axis=1)
 
dLW_net_surf_zmean = dLW_t_surf_zmean + dLW_q_surf_zmean + dLWCRE_adj_zmean

#calculate clear-sky net radiation from MERRA2

LW_net_surf_anom_zmean = np.ma.average(LW_net_surf_anom[:,wi:ei],axis=1)

SW_net_surf_anom_zmean = np.ma.average(SW_net_surf_anom[:,wi:ei],axis=1)

dSW_net_surf_zmean =  dSW_q_surf_zmean + dSW_alb_zmean + dSWCRE_adj_zmean


#TODO: Compute feedbacks cloud component of radiation (divide by SST pattern?)

############## PLOTTING
print 'plotting...'


latsplot = lats

#Plot surface clear-sky radiation decomposition
f, axarr = plt.subplots(2,1, figsize=(14,12))
axarr[0,].plot(dLW_t_surf_cs_zmean.getValue(), latsplot,  color='C0', label=r'$\Delta LW^{clear}_{T}$')
axarr[0,].plot(dLW_q_surf_cs_zmean.getValue(), latsplot,  color='C1', label=r'$\Delta LW^{clear}_{q}$')
axarr[0,].plot(dLW_net_surf_cs_zmean.getValue(), latsplot,  color='grey', label=r'sum')
axarr[0,].plot(LW_net_surf_cs_anom_zmean, latsplot,  color='k', label=r'$\Delta LW^{clear}$')
axarr[0,].axhline(0, color='k', alpha=0.5, linewidth=1)
axarr[0,].set_xlabel(r'Clear-sky Net Surface Longwave (W/m$^{2}$)')
axarr[0,].legend(loc='upper right')
axarr[0,].axvline(0, linewidth=1, color='k')
axarr[0,].set_xlim(-0.8,1.4)
axarr[0,].set_ylabel('Latitude (degrees)')
axarr[1,].plot(dSW_q_surf_cs_zmean.getValue(), latsplot, color='C1', label=r'$\Delta SW^{clear}_{q}$')
axarr[1,].plot(dSW_alb_cs_zmean.getValue(), latsplot, color='C3', label=r'$\Delta SW^{clear}_{\alpha}$')
axarr[1,].plot(dSW_net_surf_cs_zmean.getValue(), latsplot, color='grey', label=r'sum')
axarr[1,].plot(SW_net_surf_cs_anom_zmean, latsplot, color='k', label=r'$\Delta SW^{clear}$')
axarr[1,].axhline(0, color='k', alpha=0.5, linewidth=1)
#axarr[1,].set_title(r'Clear-sky Shortwave')
axarr[1,].set_xlabel(r'Clear-sky Net Surface Shortwave (W/m$^{2}$)')
axarr[1,].set_ylabel('Latitude (degrees)')
axarr[1,].legend(loc='upper right')
axarr[1,].axvline(0, linewidth=1, color='k')
axarr[1,].set_xlim(-0.8,1.4)
plt.savefig(fout + 'MERRA2_SST5to20_clearsky_kerneldecomp.pdf')
plt.close()

#TODO: Calculate Adjusted CRE (shortwave and longwave)

f, axarr = plt.subplots(2,1, figsize=(14,12))
axarr[0,].plot(dLW_t_surf_zmean.getValue(), latsplot,  color='C0', label=r'$\Delta LW^{all}_{T}$')
axarr[0,].plot(dLW_q_surf_zmean.getValue(), latsplot,  color='C1', label=r'$\Delta LW^{all}_{q}$')
axarr[0,].plot(dLWCRE_adj_zmean, latsplot,  color='C4', label=r'$\Delta LW_{cloud}$')
axarr[0,].plot(dLW_net_surf_zmean.getValue(), latsplot,  color='grey', label=r'sum')
axarr[0,].plot(LW_net_surf_anom_zmean, latsplot,  color='k', label=r'$\Delta LW^{all}$')
axarr[0,].axhline(0, color='k', alpha=0.5, linewidth=1)
axarr[0,].set_xlabel(r'Net Surface Longwave (W/m$^{2}$)')
axarr[0,].legend(loc='upper right')
axarr[0,].axvline(0, linewidth=1, color='k')
axarr[0,].set_xlim(-2,2)
axarr[0,].set_ylabel('Latitude (degrees)')
axarr[1,].plot(dSW_q_surf_zmean.getValue(), latsplot, color='C1', label=r'$\Delta SW^{all}_{q}$')
axarr[1,].plot(dSW_alb_zmean.getValue(), latsplot, color='C3', label=r'$\Delta SW^{all}_{\alpha}$')
axarr[1,].plot(dSWCRE_adj_zmean, latsplot,  color='C4', label=r'$\Delta SW_{cloud}$')
axarr[1,].plot(dSW_net_surf_zmean.getValue(), latsplot, color='grey', label=r'sum')
axarr[1,].plot(SW_net_surf_anom_zmean, latsplot, color='k', label=r'$\Delta SW^{all}$')
axarr[1,].axhline(0, color='k', alpha=0.5, linewidth=1)
#axarr[1,].set_title(r'Clear-sky Shortwave')
axarr[1,].set_xlabel(r'Net Surface Shortwave (W/m$^{2}$)')
axarr[1,].set_ylabel('Latitude (degrees)')
axarr[1,].legend(loc='upper right')
axarr[1,].axvline(0, linewidth=1, color='k')
axarr[1,].set_xlim(-2,2)
plt.savefig(fout + 'MERRA2_SST5to20_allsky_kerneldecomp.pdf')
plt.close()




































