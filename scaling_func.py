#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  3 15:14:48 2020

@author: wieka20
"""
import numpy as np
from netCDF4 import Dataset

def upscaling(lat_h,lon_h,lat_c,lon_c):
  r_sel=[]
  for i in range(len(lat_h)):
    dlat_min=1000.
    row_min=0
    for j in range(len(lat_c)):
      dlat=abs(lat_h[i]-lat_c[j])
      if dlat<dlat_min:
        dlat_min=1*dlat
        row_min=j      
    r_sel.append(row_min)
  k_sel=[]
  for i in range(len(lon_h)):
    dlon_min=1000.
    kol_min=0
    for j in range(len(lon_c)):
      dlon=abs(lon_h[i]-lon_c[j])
      if dlon<dlon_min:
        dlon_min=1*dlon
        kol_min=j      
    k_sel.append(kol_min)
  r_sel=np.asarray(r_sel) #closest coarse res rows corresponding to the high res lats of basin
  k_sel=np.asarray(k_sel)
  
  return r_sel,k_sel

def select_region(lats,lons,llx,urx,lly,ury):
    rows_sel=[]
    for i in range(len(lats)):
        if lats[i]>=lly and lats[i]<=ury:
            rows_sel.append(i)
    kols_sel=[]
    for j in range(len(lons)):
        if lons[j]>=llx and lons[j]<=urx:
            kols_sel.append(j)
    rows_sel=np.asarray(rows_sel)
    kols_sel=np.asarray(kols_sel)
    return rows_sel,kols_sel
