# -*- coding: utf-8 -*-
import numpy as np
from netCDF4 import Dataset
import sys
import os.path
import gdal
#import ee 
#ee.Initialize()
#import richdem as rd

def input_LS(fname,llx,urx,lly,ury):
    
    '''
    Derive slope, aspect, and contributing area from a 30m DEM on GEE platform 
    Here I only derived slope and aspect from GEE and exported that as Gtiff images    
    '''
    
    ##################GEE part########################
    #region=ee.Geometry.Rectangle([llx,lly,urx,ury])
    #image = ee.Image('USGS/SRTMGL1_003') #30m resolution DEM SRTM void-filled
    ##print(image.getInfo())
    #Project the image to Asia_Lambert_Conformal_Conic
    #wkt = 'PROJCS["Asia_Lambert_Conformal_Conic", \
    #    GEOGCS["GCS_WGS_1984", \
    #    DATUM["WGS_1984", \
    #    SPHEROID["WGS_1984",6378137,298.257223563]], \
    #    PRIMEM["Greenwich",0], \
    #    UNIT["Degree",0.017453292519943295]], \
    #    PROJECTION["Lambert_Conformal_Conic_2SP"], \
    #    PARAMETER["False_Easting",0], \
    #    PARAMETER["False_Northing",0],\
    #    PARAMETER["Central_Meridian",105],\
    #    PARAMETER["Standard_Parallel_1",30],\
    #    PARAMETER["Standard_Parallel_2",62],\
    #    PARAMETER["Latitude_Of_Origin",0],\
    #    UNIT["Meter",1], AUTHORITY["EPSG","102012"]]'

    #proj_Asia_Lambert = ee.Projection(wkt);
    #image_Asia_Lambert = image.reproject(crs=proj_Asia_Lambert,scale= 30);

    #elevation = image_Asia_Lambert.select('elevation')
    #slope = ee.Terrain.slope(elevation)
    #aspect = ee.Terrain.aspect(elevation)
    ##################################################
    
    # read slope data
    file_in_slope = '%s/Data_input_test/topo/Slope_China_30m_test.tif' % (fname)
    # file_out_slope= '%s/topo/Slope_China_30m_xxx.nc' % (fname)
    # !gdal_translate -of 'netCDF' file_in_slope file_out_slope
    ds = gdal.Open(file_in_slope)
    band = ds.GetRasterBand(1)
    slope = band.ReadAsArray()
    band=None
    ds=None
    # read aspect data
    file_in_aspect = '%s/Data_input_test/topo/Aspect_China_30m_test.tif' % (fname)
    ds = gdal.Open(file_in_aspect)
    band = ds.GetRasterBand(1)
    aspect = band.ReadAsArray()
    band=None
    ds=None
    # read contributing area data
    file_in_A_contr = '%s/Data_input_test/topo/A_contr_China_30m_test.tif' % (fname)
    ds = gdal.Open(file_in_A_contr)
    band = ds.GetRasterBand(1)
    A_contr = band.ReadAsArray()
    band=None
    ds=None
    
    return slope,aspect,A_contr

def input_K(fname,llx,urx,lly,ury):
    '''
    Read the soil input data for the calculation of the K factor at 1 km spatial resolution from GSDE (Shangguan et al., 2014)
    Extract data for test region
    '''
    
    from scaling_func import select_region
    data = Dataset('%s/Data_input_test/soil/SAND1.nc' % (fname),'r')   
    lats=data.variables['lat'][:]
    lons=data.variables['lon'][:]
    rows_sel,kols_sel=select_region(lats,lons,llx,urx,lly,ury) 
    lat_c=lats[rows_sel[0]:rows_sel[-1]+1]
    lon_c=lons[kols_sel[0]:kols_sel[-1]+1]
    
    sand = data.variables['SAND'][3,rows_sel[0]:rows_sel[-1]+1,kols_sel[0]:kols_sel[-1]+1] #select topsoil layer for all, units % 
    sand_t = sand.astype(float)
    sand_t[sand_t<0.] = np.nan
    sand_t = sand_t*0.01 # now fraction
    data.close()

    data = Dataset('%s/Data_input_test/soil/SILT1.nc' % (fname) ,'r')
    silt = data.variables['SILT'][3,rows_sel[0]:rows_sel[-1]+1,kols_sel[0]:kols_sel[-1]+1] 
    silt_t = silt.astype(float)
    silt_t[silt_t<0.] = np.nan
    silt_t = silt_t*0.01 
    data.close()

    data = Dataset('%s/Data_input_test/soil/CLAY1.nc' % (fname),'r')
    clay = data.variables['CLAY'][3,rows_sel[0]:rows_sel[-1]+1,kols_sel[0]:kols_sel[-1]+1] 
    clay_t = clay.astype(float)
    clay_t[clay_t<0.] = np.nan
    clay_t = clay_t*0.01 
    data.close()

    data = Dataset('%s/Data_input_test/soil/OC1.nc' % (fname),'r')
    oc = data.variables['OC'][3,rows_sel[0]:rows_sel[-1]+1,kols_sel[0]:kols_sel[-1]+1] 
    oc_t = oc.astype(float)
    oc_t[oc_t<0.] = np.nan
    oc_t = oc_t*0.01 # to have % weight
    oc_t = oc_t*0.01 # to have fractions!
    data.close()
    
    #soil type classification: Taxo USDA 2014 from soil grids at 1km;
    file_in = '%s/Data_input_test/soil/TAXOUSDA_China_1km_test.tif' % (fname)
    ds = gdal.Open(file_in)
    band = ds.GetRasterBand(1)
    su_code90 = band.ReadAsArray()
    band=None
    ds=None

    return sand_t,silt_t,clay_t,oc_t,su_code90,lat_c,lon_c

def input_grav(fname,llx,urx,lly,ury,lat_h,lon_h):
    '''
    Read the gravel data at 1 km spatial resolution from GSDE (Shangguan et al., 2014) 
    then upscale to 30m for Eros function in Adj_RUSLE.py
    Extract data for test region
    '''
    
    from scaling_func import select_region,upscaling
    
    data = Dataset('%s/Data_input_test/soil/GRAV1.nc' % (fname),'r')    
    lats=data.variables['lat'][:]
    lons=data.variables['lon'][:]
    rows_sel,kols_sel=select_region(lats,lons,llx,urx,lly,ury)    
    lat_c=lats[rows_sel[0]:rows_sel[-1]+1]
    lon_c=lons[kols_sel[0]:kols_sel[-1]+1]
    
    gravel = data.variables['GRAV'][3,rows_sel[0]:rows_sel[-1]+1,kols_sel[0]:kols_sel[-1]+1] # units %
    gravel=gravel.astype(float)
    gravel[gravel<0.] = 0.
    gravel[np.isnan(gravel)==True] = 0.
    data.close()
    #rescale gravel to 30m
    r_sel,k_sel=upscaling(lat_h,lon_h,lat_c,lon_c)
    gravel_out1=gravel[r_sel[:]];gravel_out=gravel_out1[:,k_sel[:]]
    
    return gravel_out

def input_lc(fname,llx,urx,lly,ury):  
    '''
    Derives landcover at 100m res on GEE platform 
    Outputs landcover data test region as Gtiff images
    Landcover fractions are then rescaled to 30m for calculation of the C factor 
    and as input for the Eros function in Adj_RUSLE.py
       
    '''
    from scaling_func import upscaling
    
    ##################GEE part########################
    #region=ee.Geometry.Rectangle([llx,lly,urx,ury])

    ##landcover data
    #imagecol=ee.ImageCollection("COPERNICUS/Landcover/100m/Proba-V/Global").filterDate('2015-01-01')
    #imagereg=imagecol.clip(region)
    #crops = imagereg.select(['crops-coverfraction'])
    #tree= imagereg.select(['tree-coverfraction'])
    #grass=imagereg.select(['grass-coverfraction'])
    #baresoil=imagereg.select(['baresoil-coverfraction'])
    #urban=imagereg.select(['urban-coverfraction'])

    ###################################################
    
    # read ndvi data to derive high-res lat and lon (lat_h,lon_h)
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
    # read landcover data
    file_in_tree = '%s/Data_input_test/landcover/Tree_China_100m_test_2015.tif' % (fname)
    ds = gdal.Open(file_in_tree)
    band = ds.GetRasterBand(1)
    treefrac_c = band.ReadAsArray()*0.01
    width = ds.RasterXSize
    height = ds.RasterYSize
    gt = ds.GetGeoTransform()
    minx = gt[0]
    miny = gt[3] + width*gt[4] + height*gt[5] 
    maxx = gt[0] + width*gt[1] + height*gt[2]
    maxy = gt[3] 
    lat_c=np.arange(miny,maxy,gt[5]*-1)
    lon_c=np.arange(minx,maxx,gt[1])
    band=None
    ds=None
    #rescale tree to 30m
    r_sel,k_sel=upscaling(lat_h,lon_h,lat_c,lon_c)
    treefrac1=treefrac_c[r_sel[:]];treefrac=treefrac1[:,k_sel[:]]
    
    file_in_crop = '%s/Data_input_test/landcover/Crops_China_100m_test_2015.tif' % (fname)
    ds = gdal.Open(file_in_crop)
    band = ds.GetRasterBand(1)
    cropfrac_c = band.ReadAsArray()*0.01
    band=None
    ds=None
    #rescale crop to 30m
    cropfrac1=cropfrac_c[r_sel[:]];cropfrac=cropfrac1[:,k_sel[:]]
    
    file_in_grass = '%s/Data_input_test/landcover/Grass_China_100m_test_2015.tif' % (fname)
    ds = gdal.Open(file_in_grass)
    band = ds.GetRasterBand(1)
    grassfrac_c = band.ReadAsArray()*0.01
    band=None
    ds=None
    #rescale grass to 30m
    grassfrac1=grassfrac_c[r_sel[:]];grassfrac=grassfrac1[:,k_sel[:]]
    
    return treefrac,grassfrac,cropfrac

def input_C(fname,year_start,year_end,llx,urx,lly,ury):  
    '''
    Derives ndvi (monthly) at 30m resolution on GEE platform 
    Outputs ndvi test region as Gtiff images       
    '''
    from scaling_func import upscaling
    
    ##################GEE part########################
    #region=ee.Geometry.Rectangle([llx,lly,urx,ury])
    #months=np.array(['01','02','03','04','05','06','07','08','09','10','11','12'])
    #enddate=np.array(['31','28','31','30','31','30','31','31','30','31','30','31'])

    ##ndvi data
    #for m in range(len(months)): 
    #imagecol = ee.ImageCollection('LANDSAT/LC08/C01/T1_TOA')    
    #img = imagecol.filter(ee.Filter.calendarRange(%i,%i,'year') %(year_start,year_end))\
    #    .filter(ee.Filter.calendarRange(m+1,m+1,'month'))
    #imagereg=img.filterBounds(region)
    ##print('spatialFiltered', imagereg.getInfo())
    #median = imagereg.median()
    #ndvi = median.normalizedDifference(['B5', 'B4']).rename('NDVI')
    ###################################################
    
    # read ndvi data
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
    ndvi_all=[]
    for m in range(1,13):
        file_in = '%s/Data_input_test/ndvi/NDVI_China_30m_test_%i_2012-2017.tif' % (fname,m)
        ds = gdal.Open(file_in)
        band = ds.GetRasterBand(1)
        ndvi = band.ReadAsArray()
        band=None
        ds=None
        ndvi_all.append(ndvi)
    ndvi_all=np.asarray(ndvi_all).reshape(12,ndvi.shape[0],ndvi.shape[1])
 
    return ndvi_all

def input_R(fname,llx,urx,lly,ury):
    
    '''
    Derives readily available climate and elevation data for calculation of erosivity
    Total yearly precipitation is derived from GSMaP on GEE 
    SDII is calculated from the daily precip exported from GEE on local machine
    high-resolution (0.005deg) precip data is the interpolated to 10km with bilinear method using cdo's
    '''
    
    from scaling_func import select_region
 
    ##################GEE part########################
    #months=np.array(['01','02','03','04','05','06','07','08','09','10','11','12'])
    #for m in range(len(months)): 
    #    if m==0 or m==2 or m==4 or m==6 or m==7 or m==9 or m==11:
    #        days=np.arange(1,32)
    #    elif m==1:
    #        days=np.arange(1,29)
    #    else:
    #        days=np.arange(1,31)
    #    #loop over days
    #    inidate = ee.Date.fromYMD(2012,m+1,1)
    #    if m<11:
    #        enddate = ee.Date.fromYMD(2012,m+2,1) 
    #    else:
    #        enddate = ee.Date.fromYMD(2013,1,1) 
    #    #Import GSMaP data
    #    gsmap = ee.ImageCollection('JAXA/GPM_L3/GSMaP/v6/reanalysis')\
    #        .filterDate(inidate, enddate).select('hourlyPrecipRateGC');
    #    gsmap_reg=gsmap.filterBounds(region)
    #    for d in days:
    #        if d<10:
    #            day=ee.Date('2012-%s-0%i' % (months[m],d))
    #            print(day)
    #        else:
    #            day=ee.Date('2012-%s-%i' % (months[m],d))
    #            print(day)
    #        day_collection = gsmap_reg.filterDate(day, day.advance(1 , 'day'))
    #        sum_d = ee.Image(day_collection.sum())
    # 
    #        task_config = { 'description' : 'ImagetoDriveM%sD%i' %(months[m],d), 'folder' : 'GEE', \
    #                       'fileNamePrefix' : 'PrM%sD%i' %(months[m],d)\
    #                           , 'dimensions' : 400, 'region' : geometry }
    #
    #        task = ee.batch.Export.image.toDrive(sum_d, **task_config)
    #
    #        task.start()
    ##################################################
    
    #input readily calculated coarse res precip data rescaled to 10km resolution
    data = Dataset('%s/Data_input_test/erosivity/climate_classification_5m.nc' % (fname),'r') #koppen-geiger climate classification 
    lats=data.variables['lat'][:]
    lons=data.variables['lon'][:]
    rows_sel,kols_sel=select_region(lats,lons,llx,urx,lly,ury)    
    lat_c=lats[rows_sel[0]:rows_sel[-1]+1]
    lon_c=lons[kols_sel[0]:kols_sel[-1]+1]
    k = data.variables['climate'][rows_sel[0]:rows_sel[-1]+1,kols_sel[0]:kols_sel[-1]+1]
    k[k<0.]=np.nan
    k = np.ravel(k)
    data.close()
    data = Dataset('%s/Data_input_test/erosivity/P_total_2012_5m.nc' % (fname),'r') #total yearly 50 year mean precipitation data
    precip = data.variables['P'][:]
    precip[precip<=0.]=np.nan
    precip = np.ravel(precip)
    data.close()
    # SDII (simple precipitation index data) is the sum of precip on wet days (precip>1mm) divided by the number of wet days in a period
    data = Dataset('%s/Data_input_test/erosivity/SDII_2012_5m.nc' % (fname),'r') 
    sdii = data.variables['SDII'][:]
    sdii[sdii<0.]=np.nan
    sdii = np.ravel(sdii)
    data.close()
    data = Dataset('%s/Data_input_test/erosivity/etopo5_5m.nc' % (fname),'r') #ETOPO 5 arcmin elevation data
    elev = data.variables['elev'][rows_sel[0]:rows_sel[-1]+1,kols_sel[0]:kols_sel[-1]+1]
    elev[elev<0.]=-9999.
    elev = np.ravel(elev)
    data.close()
    # remove too small elevation and precip values and zeros
    z = []
    for i in range(len(elev)):
        if 0.<=elev[i]<0.001:
            z.append(0.001)
        elif elev[i]<0.:
            z.append(np.nan)
        else:
            z.append(elev[i])
    z = np.asarray(z)
    p = []
    for i in range(len(precip)):
        if 0.<precip[i]<1.:
            p.append(1.)
        elif precip[i]<0.:
            p.append(np.nan)
        else:
            p.append(precip[i])
    p = np.asarray(p)
    z = np.ma.masked_array(z,[np.isnan(x) for x in z])
    p = np.ma.masked_array(p,[np.isnan(x) for x in p])
    k = np.ma.masked_array(k,[np.isnan(x) for x in k])
  
    return z,p,k,sdii,rows_sel,kols_sel,lat_c,lon_c

def input_P(fname,llx,urx,lly,ury):
    
    '''
    Here we read in the high-resolution terracing map of develop by Le and his team
    The original terracing map is available on GEE at 30m resolution
    '''
    
    ################GEE Part################################
    #geometry = ee.Geometry.Rectangle([llx,lly,urx,ury])
    #geometry = geometry['coordinates'][0]
    #region=ee.Geometry.Rectangle([llx,lly,urx,ury])
    #
    #ChinaTerrace = ee.Image("users/cbw/TerraceResults/China_terrace_1");
    #
    #wkt = 'PROJCS["Asia_Lambert_Conformal_Conic", \
    #    GEOGCS["GCS_WGS_1984", \
    #    DATUM["WGS_1984", \
    #    SPHEROID["WGS_1984",6378137,298.257223563]], \
    #    PRIMEM["Greenwich",0], \
    #    UNIT["Degree",0.017453292519943295]], \
    #    PROJECTION["Lambert_Conformal_Conic_2SP"], \
    #    PARAMETER["False_Easting",0], \
    #    PARAMETER["False_Northing",0],\
    #    PARAMETER["Central_Meridian",105],\
    #    PARAMETER["Standard_Parallel_1",30],\
    #    PARAMETER["Standard_Parallel_2",62],\
    #    PARAMETER["Latitude_Of_Origin",0],\
    #    UNIT["Meter",1], AUTHORITY["EPSG","102012"]]'
    #
    #proj_Asia_Lambert = ee.Projection(wkt);
    #image_Asia_Lambert = ChinaTerrace.reproject(crs=proj_Asia_Lambert,scale= 30);
    #
    #imagereg=image_Asia_Lambert.clip(region)
    #
    #task_config = { 'description' : 'imageToDriveTer', 'folder' : 'GEE', \
    #               'fileNamePrefix' : 'Chinaterrace', 'scale' : 30, 'region' : geometry }
    #
    #task_terrace = ee.batch.Export.image.toDrive(imagereg, **task_config)
    #
    #task_terrace.start()
    ######################################
    
    # read terrace data
    file_in = '%s/Data_input_test/Terracing/Terracing_China_30m_test.tif' % (fname)
    ds = gdal.Open(file_in)
    band = ds.GetRasterBand(1)
    gt =ds.GetGeoTransform()
    pixelSizeX = gt[1]
    pixelSizeY = -gt[5]
    terrace_frac = band.ReadAsArray()/(pixelSizeX*pixelSizeY)
    band=None
    ds=None

    return terrace_frac