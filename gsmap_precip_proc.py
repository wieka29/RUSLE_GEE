#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 12 16:17:52 2020

@author: wieka20
"""
import gdal
import numpy as np
from netCDF4 import Dataset

year=2012
months=np.array(['01','02','03','04','05','06','07','08','09','10','11','12'])
fname='/home/wieka20/Dokumente/Adj_RUSLE/Data_input_test/erosivity'
# derive rows kols of precip file
file_in= '%s/PrM01D1.tif' % (fname)
ds = gdal.Open(file_in)
band = ds.GetRasterBand(1)
P = band.ReadAsArray()
rows=P.shape[0]
kols=P.shape[1]
band=None
ds=None

P_monthly=[]
wetdays=np.zeros((rows,kols))
P_wetdays=np.zeros((rows,kols))
for m in range(len(months)): 
    if m==0 or m==2 or m==4 or m==6 or m==7 or m==9 or m==11:
        dd=31
    elif m==1:
        dd=28
    else:
        dd=30
    days=np.arange(1,dd+1)
    P_daily=[]
    for d in days:
        # read daily precip data from gsmap
        file_in= '%s/PrM%sD%i.tif' % (fname,months[m],d)
        ds = gdal.Open(file_in)
        band = ds.GetRasterBand(1)
        P = band.ReadAsArray()
        wetdays[P>=1.]+=1
        P_wetdays[P>=1.]+=P[P>=1.]
        band=None
        ds=None
        P_daily.append(P)
    P_daily=np.asarray(P_daily).reshape(dd,rows,kols)
    P_monthly.append(np.nansum(P_daily,axis=0)) #monthly sum
P_monthly=np.asarray(P_monthly).reshape(12,rows,kols)
P_yearly=np.nansum(P_monthly,axis=0) #yearly sum
SDII=P_wetdays/wetdays
SDII[np.isfinite(SDII)==False]=0.
SDII[np.isnan(SDII)==True]=0.
#output
output = Dataset('%s/P_total_%i.nc' % (fname,year),'w')
output.createDimension('lat',rows)
output.createDimension('lon',kols)
output.createVariable('P','d',('lat','lon',))
output.variables['P'][:] = P_yearly
output.close() 

output = Dataset('%s/SDII_%i.nc' % (fname,year),'w')
output.createDimension('lat',rows)
output.createDimension('lon',kols)
output.createVariable('SDII','d',('lat','lon',))
output.variables['SDII'][:] = SDII
output.close() 