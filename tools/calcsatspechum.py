#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 14:44:17 2018

@author: cpatrizio
"""

import numpy as np

def calcsatspechum(t,p):

## T is temperature in Kelvins, P is pressure in hPa 

  ## Formulae from Buck (1981):
  es = (1.0007+(3.46e-6*p))*6.1121.*np.exp(17.502*(t-273.15)/(240.97+(t-273.15)))
  wsl = .622*es./(p-es); %# saturation mixing ratio wrt liquid water (g/kg)
  
  es = (1.0003+(4.18e-6*p))*6.1115.*np.exp(22.452*(t-273.15)/(272.55+(t-273.15)))

  wsi = .622*es./(p-es) %# saturation mixing ratio wrt ice (g/kg)
  
  ws = wsl
  ws[t<273.15]=wsi[t<273.15]
 
  return qs=ws/(1+ws)