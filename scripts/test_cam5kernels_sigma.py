#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  4 10:41:53 2018

@author: cpatrizio
"""

import sys
sys.path.append("/Users/cpatrizio/repos/")
import cdms2 as cdms2
import glob
import numpy as np
import MV2 as MV2
from AMO.misc_fns import spatial_ave
fin = '/Users/cpatrizio/repos/cam5-kernels/kernels/'
fdata = '/Users/cpatrizio/repos/cam5-kernels/demodata/'

#model_names =['MRI-CGCM3']
#
#print 'model:', model_names[0]
#
#if model_names[0] == 'MIROC5':
#    year1='21'
#    year2='22'
#else:
#    year1='18'
#    year2='19'




fbasefields = glob.glob(fdata + 'basefields.nc')[0]
#fchangefields3D = glob.glob(fdata + 'changefields.plev.nc')[0]
fchangefields = glob.glob(fdata + 'changefields.nc')[0]



#NCL SCRIPTS TO RUN BEFORE USING KERNELS:
#convert kernels to input data pressure levels:
#q_kernel_to_plev.ncl
#t_kernel_to_plev.ncl
#calculate pressure thickness of each layer:
#calcdp_plev.ncl

basefields = cdms2.open(fbasefields)
#changefields3D = cdms2.open(fchangefields3D)
changefields = cdms2.open(fchangefields)
#pdiff = cdms2.open(fdata + 'dp_plev.nc')['dp'][:]/100.
p = cdms2.open(fdata + 'p_sigma.nc')['pmid'][:]/100.
#pdiff_test = changefields3D.variables['pdiff'][:]

falbkern = cdms2.open(fin + 'alb.kernel.nc')
fqkern = cdms2.open(fin + 'q.kernel.nc')
ftkern = cdms2.open(fin + 't.kernel.nc')
#ftkerngw = cdms2.open(fin + 't.kernel.nc')
ftskern = cdms2.open(fin + 'ts.kernel.nc')

#NOTE: plev kernels have units W m^-2 hPa^-1

FLNSt = ftkern('FLNS')
#gw = ftkerngw('gw')

lats = FLNSt.getLatitude()[:]
lons = FLNSt.getLongitude()[:]
t = FLNSt.getTime()[:]
nlat = lats.size
nlon = lons.size
nt = t.size
#p_Pa = ftkern('plev')
nplevs = p.shape[1]
#p = np.zeros((nplevs, nlat, nlon)).T
#p[:,:,:] = p_Pa/100.
#p = p.T
#p = np.repeat(p[np.newaxis,...], nt, axis=0)

#crude tropopause
x=np.cos(np.deg2rad(lats))
p_tropopause_zonalmean=300-200*x

p_tropopause = np.zeros((nlat, nlon)).T
p_tropopause[:,:] = p_tropopause_zonalmean
p_tropopause = p_tropopause.T
p_tropopause = np.repeat(p_tropopause[np.newaxis,...],nplevs,axis=0)
p_tropopause = np.repeat(p_tropopause[np.newaxis,...],nt,axis=0)

#TEMPERATURE FEEDBACKS
dts=changefields('ts')
dts_globalmean = MV2.average(spatial_ave(dts,lats))

ts_kernel_TOA = ftskern('FLNT')
ts_kernel_surf = ftskern('FLNS')
ts_kernel_TOA_cs = ftskern('FLNTC')
ts_kernel_surf_cs = ftskern('FLNSC')

dLW_ts_TOA=ts_kernel_TOA*dts  
dLW_ts_surf = ts_kernel_surf*dts
dLW_ts_TOA_cs = ts_kernel_TOA_cs*dts
dLW_ts_surf_cs = ts_kernel_surf_cs*dts

dta = changefields('temp')
ta_kernel_TOA = ftkern('FLNT')
ta_kernel_surf = ftkern('FLNS')
ta_kernel_TOA_cs = ftkern['FLNTC']
ta_kernel_surf_cs = ftkern['FLNSC']

mask_temp = dta.mask
#mask stratospheric values
dta = np.ma.MaskedArray(dta, mask=np.bitwise_and(mask_temp, p>=p_tropopause))

dLW_ta_TOA = np.ma.sum(ta_kernel_TOA*dta, axis=1)
dLW_ta_surf = np.ma.sum(ta_kernel_surf*dta, axis=1)
dLW_ta_TOA_cs = np.ma.sum(ta_kernel_TOA_cs*dta, axis=1)
dLW_ta_surf_cs = np.ma.sum(ta_kernel_surf_cs*dta, axis=1)

dLW_t_TOA_globalmean = MV2.average(spatial_ave(-dLW_ta_TOA-dLW_ts_TOA, lats))
dLW_t_surf_globalmean = MV2.average(spatial_ave(-dLW_ta_surf-dLW_ts_surf, lats))

t_feedback_TOA = dLW_t_TOA_globalmean/dts_globalmean
t_feedback_surf = dLW_t_surf_globalmean/dts_globalmean

print 'TOA temperature feedback (W m^-2 K^-1)', t_feedback_TOA
print 'Surface temperature feedback (W m^-2 K^-1)', t_feedback_surf
print 'ATM temperature feedback (W m^-2 K^-1)', t_feedback_TOA - t_feedback_surf

#ALBEDO FEEDBACKS
SW_sfc_net_1 = basefields('FSNS')
SW_sfc_down_1 = basefields('FSDS')
SW_sfc_net_2 = changefields('FSNS') + SW_sfc_net_1
SW_sfc_down_2 = changefields('FSDS') + SW_sfc_down_1

alb1 = 1 - SW_sfc_net_1/SW_sfc_down_1
alb2 = 1 - SW_sfc_net_2/SW_sfc_down_2

dalb=(alb2-alb1)*100;

alb_kernel_TOA=falbkern('FSNT')
alb_kernel_surf = falbkern('FSNS')
alb_kernel_TOA_cs = falbkern('FSNTC') 
alb_kernel_surf_cs = falbkern('FSNSC')

dSW_alb_TOA=alb_kernel_TOA*dalb
dSW_alb_surf = alb_kernel_surf*dalb
dSW_alb_TOA_cs = alb_kernel_TOA_cs*dalb
dSW_alb_surf_cs = alb_kernel_surf_cs*dalb

dSW_alb_TOA_globalmean=MV2.average(spatial_ave(dSW_alb_TOA,lats))
dSW_alb_surf_globalmean=MV2.average(spatial_ave(dSW_alb_surf,lats))

alb_feedback_TOA=dSW_alb_TOA_globalmean/dts_globalmean
alb_feedback_surf=dSW_alb_surf_globalmean/dts_globalmean
print '---------------------------------'
print 'TOA surface albedo feedback (W m^-2 K^-1)', alb_feedback_TOA
print 'Surface surface albedo feedback (W m^-2 K^-1)', alb_feedback_surf
print 'ATM surface albedo feedback (W m^-2 K^-1)', alb_feedback_TOA - alb_feedback_surf


#WATER VAPOR FEEDBACKS

#% Calculate the change in moisture for 1 K warming at constant relative humidity. 
#% Run the accompanying NCL script with your input files, or
#% implement here.                                                                                                                          
#dq1k=ncread('dq1k.plev.nc','dq1k');


#BETTER TO CALCULATE DQ1K WITH MY OWN SCRIPTS
fdq1k = cdms2.open(glob.glob(fdata + 'dq1k.nc')[0])
dq1k = fdq1k('dq1k')
#
#% Read kernels
q_LW_kernel_TOA=fqkern('FLNT')
q_SW_kernel_TOA=fqkern('FSNT')
q_LW_kernel_surf=fqkern('FLNS')
q_SW_kernel_surf=fqkern('FSNS')

q_LW_kernel_TOA_cs=fqkern['FLNTC']
q_SW_kernel_TOA_cs=fqkern['FSNTC']
q_LW_kernel_surf_cs=fqkern['FLNSC']
q_SW_kernel_surf_cs=fqkern['FSNSC']

#% Normalize kernels by the change in moisture for 1 K warming at
# constant RH
q_LW_kernel_TOA=q_LW_kernel_TOA/dq1k
q_SW_kernel_TOA=q_SW_kernel_TOA/dq1k
q_LW_kernel_surf=q_LW_kernel_surf/dq1k
q_SW_kernel_surf=q_SW_kernel_surf/dq1k

q_LW_kernel_TOA_cs=q_LW_kernel_TOA_cs/dq1k
q_SW_kernel_TOA_cs=q_SW_kernel_TOA_cs/dq1k
q_LW_kernel_surf_cs=q_LW_kernel_surf_cs/dq1k
q_SW_kernel_surf_cs=q_SW_kernel_surf_cs/dq1k


#
#% Read the change in moisture
dq=changefields['Q'][:]

mask_temp = dq.mask
#mask stratospheric values
dq = np.ma.MaskedArray(dq, mask=np.bitwise_and(mask_temp, p>=p_tropopause))
#
#% Convolve moisture kernel with change in moisture
dLW_q_TOA=np.ma.sum(q_LW_kernel_TOA*dq, axis=1)
dSW_q_TOA=np.ma.sum(q_SW_kernel_TOA*dq, axis=1)
dLW_q_surf=np.ma.sum(q_LW_kernel_surf*dq, axis=1)
dSW_q_surf=np.ma.sum(q_SW_kernel_surf*dq, axis=1)

dLW_q_TOA_cs=np.ma.sum(q_LW_kernel_TOA_cs*dq, axis=1)
dSW_q_TOA_cs=np.ma.sum(q_SW_kernel_TOA_cs*dq, axis=1)
dLW_q_surf_cs=np.ma.sum(q_LW_kernel_surf_cs*dq, axis=1)
dSW_q_surf_cs=np.ma.sum(q_SW_kernel_surf_cs*dq, axis=1)

#% Add the LW and SW responses. Note the sign convention difference
#% between LW and SW!
dR_q_TOA_globalmean=MV2.average(spatial_ave(-dLW_q_TOA+dSW_q_TOA,lats))
dR_q_surf_globalmean=MV2.average(spatial_ave(-dLW_q_surf+dSW_q_surf,lats))                   
#
#% Divide by the global annual mean surface warming (units: W/m2/K)
q_feedback_TOA=dR_q_TOA_globalmean/dts_globalmean
q_feedback_surf=dR_q_surf_globalmean/dts_globalmean

print '---------------------------------'
print 'TOA water vapor feedback (W m^-2 K^-1)', q_feedback_TOA
print 'Surface water vapor feedback (W m^-2 K^-1)', q_feedback_surf
print 'ATM water vapor feedback (W m^-2 K^-1)', q_feedback_TOA - q_feedback_surf

#%%%%% Change in Cloud Radiative Effect (CRE) 
d_sw_TOA=changefields('FSNT')
d_lw_TOA=changefields('FLNT')
d_sw_TOA_cs=changefields('FSNTC')
d_lw_TOA_cs=changefields('FLNTC')

d_sw_surf=changefields('FSNS')
d_lw_surf=changefields('FLNS')
d_sw_surf_cs=changefields('FSNSC')
d_lw_surf_cs=changefields('FLNSC')

d_cre_sw_TOA=d_sw_TOA_cs-d_sw_TOA
d_cre_lw_TOA=d_lw_TOA_cs-d_lw_TOA

d_cre_sw_surf=d_sw_surf_cs-d_sw_surf
d_cre_lw_surf=d_lw_surf_cs-d_lw_surf
#
#
#%%%% Cloud masking of radiative forcing
ghgfile=glob.glob(fdata + 'ghg.forcing.nc')[0]
ghgfields = cdms2.open(ghgfile)

sw_TOA=ghgfields('FSNT')
lw_TOA=ghgfields('FLNT')
sw_TOA_cs=ghgfields('FSNTC')
lw_TOA_cs = ghgfields('FLNTC')
ghg_sw_TOA=sw_TOA_cs-sw_TOA
ghg_lw_TOA=lw_TOA_cs-lw_TOA

sw_surf=ghgfields('FSNS')
lw_surf=ghgfields('FLNS')
sw_surf_cs=ghgfields('FSNSC')
lw_surf_cs = ghgfields('FLNSC')
ghg_sw_surf=sw_surf_cs-sw_surf
ghg_lw_surf=lw_surf_cs-lw_surf

aerosolfile=glob.glob(fdata + 'aerosol.forcing.nc')[0]
aerosolfields = cdms2.open(aerosolfile)
sw_TOA=aerosolfields('FSNT')
sw_TOA_cs=aerosolfields('FSNTC')
lw_TOA=aerosolfields('FLNT')
lw_TOA_cs=aerosolfields('FLNTC')
aerosol_sw_TOA=sw_TOA_cs-sw_TOA
aerosol_lw_TOA=lw_TOA_cs-lw_TOA

sw_surf=aerosolfields('FSNT')
lw_surf=aerosolfields('FLNS')
sw_surf_cs=aerosolfields('FSNSC')
lw_surf_cs=aerosolfields('FLNSC')
aerosol_lw_surf=lw_surf_cs-lw_surf
aerosol_sw_surf = sw_surf_cs - sw_surf

cloud_masking_of_forcing_sw_TOA=aerosol_sw_TOA+ghg_sw_TOA
cloud_masking_of_forcing_lw_TOA=aerosol_lw_TOA+ghg_lw_TOA
cloud_masking_of_forcing_sw_surf=aerosol_sw_surf+ghg_sw_surf
cloud_masking_of_forcing_lw_surf=aerosol_lw_surf+ghg_lw_surf

#%%%%%% Cloud feedback. 
#%%% CRE + cloud masking of radiative forcing + corrections for each feedback
#
dLW_cloud_TOA=-d_cre_lw_TOA+cloud_masking_of_forcing_lw_TOA+(dLW_q_TOA_cs-dLW_q_TOA)+(dLW_ta_TOA_cs-dLW_ta_TOA)+(dLW_ts_TOA_cs-dLW_ts_TOA)
dSW_cloud_TOA=-d_cre_sw_TOA+cloud_masking_of_forcing_sw_TOA+(dSW_q_TOA_cs-dSW_q_TOA)+(dSW_alb_TOA_cs-dSW_alb_TOA)
dLW_cloud_surf=-d_cre_lw_surf+cloud_masking_of_forcing_lw_surf+(dLW_q_surf_cs-dLW_q_surf)+(dLW_ta_surf_cs-dLW_ta_surf)+(dLW_ts_surf_cs-dLW_ts_surf)
dSW_cloud_surf=-d_cre_sw_surf+cloud_masking_of_forcing_sw_surf+(dSW_q_surf_cs-dSW_q_surf)+(dSW_alb_surf_cs-dSW_alb_surf)

#%Take global and annual averages
dLW_cloud_TOA_globalmean=MV2.average(spatial_ave(-dLW_cloud_TOA, lats))
dSW_cloud_TOA_globalmean=MV2.average(spatial_ave(dSW_cloud_TOA, lats))
dLW_cloud_surf_globalmean=MV2.average(spatial_ave(-dLW_cloud_surf, lats))
dSW_cloud_surf_globalmean=MV2.average(spatial_ave(dSW_cloud_surf, lats))
#
#%Divide by global, annual mean temperature change to get W/m2/K
lw_cloud_feedback_TOA=dLW_cloud_TOA_globalmean/dts_globalmean
sw_cloud_feedback_TOA=dSW_cloud_TOA_globalmean/dts_globalmean
lw_cloud_feedback_surf=dLW_cloud_surf_globalmean/dts_globalmean
sw_cloud_feedback_surf=dSW_cloud_surf_globalmean/dts_globalmean

print '---------------------------------'
print 'TOA LW cloud feedback (W m^-2 K^-1)', lw_cloud_feedback_TOA
print 'Surface LW cloud feedback (W m^-2 K^-1)', lw_cloud_feedback_surf
print 'ATM LW cloud feedback (W m^-2 K^-1)', lw_cloud_feedback_TOA - lw_cloud_feedback_surf
print '---------------------------------'
print 'TOA SW cloud feedback (W m^-2 K^-1)', sw_cloud_feedback_TOA
print 'Surface SW cloud feedback (W m^-2 K^-1)', sw_cloud_feedback_surf
print 'ATM SW cloud feedback (W m^-2 K^-1)', sw_cloud_feedback_TOA - sw_cloud_feedback_surf

#
#disp(['LW Cloud Feedback: ' num2str(lw_cloud_feedback) ' W m^-2 K^-1'])
#disp(['SW Cloud Feedback: ' num2str(sw_cloud_feedback) ' W m^-2 K^-1'])








