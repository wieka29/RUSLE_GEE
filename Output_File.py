#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  6 18:21:53 2020

@author: wieka20
"""
import numpy as np
from netCDF4 import Dataset
from Adj_RUSLE import Eros,LS_factor,R_factor,C_factor,K_factor,P_factor
from Input_File import input_grav,input_lc
import gdal
from netCDF4 import Dataset
#import rasterio
#from rasterio import Affine as A
#from rasterio.warp import calculate_default_transform,reproject, Resampling
#select coordinates test region
llx = 106.;lly = 36.;urx = 108.;ury = 38.
#select years for ndvi
year_start=2012;year_end=2017
#select work directory
fname='/home/wieka20/Dokumente/Adj_RUSLE/'
# read high resolution coordinates
file_in = '%s/Data_input_test/ndvi/NDVI_China_30m_test_1_2012-2017.tif' % (fname)
ds = gdal.Open(file_in)
width = ds.RasterXSize
height = ds.RasterYSize
gt = ds.GetGeoTransform()
minx = gt[0]
miny = gt[3] + width*gt[4] + height*gt[5] 
maxx = gt[0] + width*gt[1] + height*gt[2]
maxy = gt[3] 
lat_h=np.arange(miny,maxy,gt[5]*-1)
lon_h=np.arange(minx,maxx,gt[1])
band=None
ds=None
#derive RUSLe factors
LS_out=LS_factor(fname,llx,urx,lly,ury)
K_out=K_factor(fname,llx,urx,lly,ury,lat_h,lon_h)
R_out=R_factor(fname,llx,urx,lly,ury,lat_h,lon_h)
C_out=C_factor(fname,year_start,year_end,llx,urx,lly,ury)
P_out=P_factor(fname,llx,urx,lly,ury)
gravel=input_grav(fname,llx,urx,lly,ury,lat_h,lon_h)
treefrac,grassfrac,cropfrac=input_lc(fname,llx,urx,lly,ury)
#reproject C, K, R factors from WGS84 to Asia_Lambert_Conformal_Conic (EPSG:102012) on local machine 
#reprojecting using rasterio on local machine (Reprojecting these rasters within python) does not work somehow
#output factors to netcdf
output = Dataset('%s/Data_output_test/C_factor_China_30m_degree.nc' % (fname),'w')
output.createDimension('lat',C_out.shape[0])
output.createDimension('lon',C_out.shape[1])
output.createVariable('C_out','d',('lat','lon',))
output.variables['C_out'][:] = C_out
output.close()
output = Dataset('%s/Data_output_test/K_factor_China_30m_degree.nc' % (fname),'w')
output.createDimension('lat',K_out.shape[0])
output.createDimension('lon',K_out.shape[1])
output.createVariable('K_out','d',('lat','lon',))
output.variables['K_out'][:] = K_out
output.close()
output = Dataset('%s/Data_output_test/R_factor_China_30m_degree.nc' % (fname),'w')
output.createDimension('lat',R_out.shape[0])
output.createDimension('lon',R_out.shape[1])
output.createVariable('R_out','d',('lat','lon',))
output.variables['R_out'][:] = R_out
output.close()
output = Dataset('%s/Data_output_test/treefrac_China_30m_degree.nc' % (fname),'w')
output.createDimension('lat',treefrac.shape[0])
output.createDimension('lon',treefrac.shape[1])
output.createVariable('treefrac','d',('lat','lon',))
output.variables['treefrac'][:] = treefrac
output.close()
output = Dataset('%s/Data_output_test/grassfrac_China_30m_degree.nc' % (fname),'w')
output.createDimension('lat',grassfrac.shape[0])
output.createDimension('lon',grassfrac.shape[1])
output.createVariable('grassfrac','d',('lat','lon',))
output.variables['grassfrac'][:] = grassfrac
output.close()
output = Dataset('%s/Data_output_test/cropfrac_China_30m_degree.nc' % (fname),'w')
output.createDimension('lat',cropfrac.shape[0])
output.createDimension('lon',cropfrac.shape[1])
output.createVariable('cropfrac','d',('lat','lon',))
output.variables['cropfrac'][:] = cropfrac
output.close()
output = Dataset('%s/Data_output_test/gravel_China_30m_degree.nc' % (fname),'w')
output.createDimension('lat',gravel.shape[0])
output.createDimension('lon',gravel.shape[1])
output.createVariable('grav','d',('lat','lon',))
output.variables['grav'][:] = gravel
output.close()
#on local machine do:
'''
cdo setgrid,sourcegrid.txt R_factor_China_30m_degree.nc R_factor_China_30m_degree_grid.nc
gdal_translate -of 'GTiff' R_factor_China_30m_degree_grid.nc R_factor_China_30m_degree.tif
gdalwarp -s_srs epsg:4326 -t_srs epsg:102012 -r near -of gtiff -tr 30 -30 R_factor_China_30m_degree.tif R_factor_China_30m.tif
'''
#read in factors again
file_in_C = '%s/Data_output_test/C_factor_China_30m.tif' % (fname)
ds = gdal.Open(file_in_C)
band = ds.GetRasterBand(1)
C_out = band.ReadAsArray()
C_out[C_out>100.]=0.
band=None
ds=None
file_in_R = '%s/Data_output_test/R_factor_China_30m.tif' % (fname)
ds = gdal.Open(file_in_R)
band = ds.GetRasterBand(1)
R_out = band.ReadAsArray()
R_out[R_out>1000000.]=0.
band=None
ds=None
file_in_K = '%s/Data_output_test/K_factor_China_30m.tif' % (fname)
ds = gdal.Open(file_in_K)
band = ds.GetRasterBand(1)
K_out = band.ReadAsArray()
K_out[K_out>100.]=0.
band=None
ds=None
file_in_tree = '%s/Data_output_test/treefrac_China_30m.tif' % (fname)
ds = gdal.Open(file_in_tree)
band = ds.GetRasterBand(1)
treefrac = band.ReadAsArray()
treefrac[treefrac>100.]=0.
band=None
ds=None
file_in_crop = '%s/Data_output_test/cropfrac_China_30m.tif' % (fname)
ds = gdal.Open(file_in_crop)
band = ds.GetRasterBand(1)
cropfrac = band.ReadAsArray()
cropfrac[cropfrac>100.]=0.
band=None
ds=None
file_in_grass = '%s/Data_output_test/grassfrac_China_30m.tif' % (fname)
ds = gdal.Open(file_in_grass)
band = ds.GetRasterBand(1)
grassfrac = band.ReadAsArray()
grassfrac[grassfrac>100.]=0.
band=None
ds=None
file_in_grav = '%s/Data_output_test/gravel_China_30m.tif' % (fname)
ds = gdal.Open(file_in_grav)
band = ds.GetRasterBand(1)
gravel = band.ReadAsArray()
gravel[gravel>100.]=0.
band=None
ds=None
#This function calculates the annual average soil erosion rate in t/ha/yr at 30m spatial resolution
E=LS_out*K_out[:,:-1]*R_out[:,:-1]*C_out[:,:-1]*P_out
E_CG=1.*E;E_rest=1.*E
E_CG[(cropfrac[:,:-1]+grassfrac[:,:-1])<1.]=0.
E_CG[gravel[:,:-1]>12.]=E_CG[gravel[:,:-1]>12.]-(0.8*E_CG[gravel[:,:-1]>12.])
E_rest[(cropfrac[:,:-1]+grassfrac[:,:-1])>=1.]=0.
E_rest[gravel[:,:-1]>30.]=E_rest[gravel[:,:-1]>30.]-(0.3*E_rest[gravel[:,:-1]>30.])
E_out=E_CG+E_rest #yearly average soil erosion rates in tonne/ha/year
E_out[E_out>100.]=100. #limit soil erosion to 100t/ha/y
#output erosion
output = Dataset('%s/Data_output_test/Erosion_China_30m.nc' % (fname),'w')
output.createDimension('lat',E_out.shape[0])
output.createDimension('lon',E_out.shape[1])
output.createVariable('E','d',('lat','lon',))
output.variables['E'][:] = E_out
output.close()
