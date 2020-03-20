# -*- coding: utf-8 -*-
"""
test GEE
"""
import ee
import numpy as np

ee.Initialize()

image = ee.Image('USGS/SRTMGL1_003') #30m resolution DEM SRTM void-filled
#print(image.getInfo())
#image.projection().nominalScale().getInfo()
#Project the image to Asia_Lambert_Conformal_Conic
wkt = 'PROJCS["Asia_Lambert_Conformal_Conic", \
    GEOGCS["GCS_WGS_1984", \
    DATUM["WGS_1984", \
    SPHEROID["WGS_1984",6378137,298.257223563]], \
    PRIMEM["Greenwich",0], \
    UNIT["Degree",0.017453292519943295]], \
    PROJECTION["Lambert_Conformal_Conic_2SP"], \
    PARAMETER["False_Easting",0], \
    PARAMETER["False_Northing",0],\
    PARAMETER["Central_Meridian",105],\
    PARAMETER["Standard_Parallel_1",30],\
    PARAMETER["Standard_Parallel_2",62],\
    PARAMETER["Latitude_Of_Origin",0],\
    UNIT["Meter",1], AUTHORITY["EPSG","102012"]]'

proj_Asia_Lambert = ee.Projection(wkt);
image_Asia_Lambert = image.reproject(crs=proj_Asia_Lambert,scale= 30);

elevation = image_Asia_Lambert.select('elevation')
slope = ee.Terrain.slope(elevation)
aspect = ee.Terrain.aspect(elevation)


'''
Test some calculations

xy = ee.Geometry.Point([86.9250, 27.9881])
elev = image.sample(xy, 30).first().get('elevation').getInfo()
print('Mount Everest elevation (m):', elev)

slope_sel=slope.sample(xy,30).first().getInfo()
print('Mount Everest slope (degree):', slope_sel)

region=ee.Geometry.Rectangle([86.76,34.69,88.5,35.52]) #minLon, minLat, maxLon, maxLat
affine = [0.0002777777777777778, 0, -180.0001388888889, 0, -0.0002777777777777778, 60.00013888888889]

meanElevation = elevation.reduceRegion(reducer=ee.Reducer.mean(),\
                                       geometry=region,crs='EPSG:4326',\
                                       crsTransform=affine,\
                                       maxPixels=1e9)
print('Mean elevation:', meanElevation.getInfo())
'''

# Export image/data
llx = 106.
lly = 36.
urx = 108.
ury = 38.
geometry = ee.Geometry.Rectangle([llx,lly,urx,ury])
geometry = geometry['coordinates'][0]

task_config_elev = { 'description' : 'imageToDriveDEM', 'folder' : 'GEE', \
               'fileNamePrefix' : 'DEMTest', 'scale' : 30, 'region' : geometry }

task_elev = ee.batch.Export.image.toDrive(elevation, **task_config_elev)

task_elev.start()

task_config_slope = { 'description' : 'imageToDriveSlope', 'folder' : 'GEE', \
               'fileNamePrefix' : 'SlopeTest', 'scale' : 30, 'region' : geometry }

task_slope = ee.batch.Export.image.toDrive(slope, **task_config_slope)

task_slope.start()

task_config_aspect = { 'description' : 'imageToDriveAspect', 'folder' : 'GEE', \
               'fileNamePrefix' : 'AspectTest', 'scale' : 30, 'region' : geometry }

task_aspect = ee.batch.Export.image.toDrive(slope, **task_config_aspect)

task_aspect.start()

"""
Test NDVI calculation
"""
region=ee.Geometry.Rectangle([llx,lly,urx,ury])
months=np.array(['01','02','03','04','05','06','07','08','09','10','11','12'])
enddate=np.array(['31','28','31','30','31','30','31','31','30','31','30','31'])
for m in range(len(months)): 
  imagecol = ee.ImageCollection('LANDSAT/LC08/C01/T1_TOA')
  #filterDate('2017-%s-01' % (months[m]), '2017-%s-%s' % (months[m],enddate[m]))     
  img = imagecol.filter(ee.Filter.calendarRange(2012,2017,'year'))\
      .filter(ee.Filter.calendarRange(m+1,m+1,'month'))
  imagereg=img.filterBounds(region)
  #print('spatialFiltered', imagereg.getInfo())
  median = imagereg.median()
  ndvi = median.normalizedDifference(['B5', 'B4']).rename('NDVI')
  #export ndvi image
  task_config_ndvi = { 'description' : 'imageToDrivendvi%s' % (months[m]), 'folder' : 'GEE', \
               'fileNamePrefix' : 'ndviTest%s' % (months[m]), 'scale' : 30, 'region' : geometry }

  task_ndvi = ee.batch.Export.image.toDrive(ndvi, **task_config_ndvi)

  task_ndvi.start()
  
"""
Derive landcover
"""
imagecol=ee.ImageCollection("COPERNICUS/Landcover/100m/Proba-V/Global").filterDate('2015-01-01').first()
imagereg=imagecol.clip(region)
crops = imagereg.select(['crops-coverfraction'])
tree= imagereg.select(['tree-coverfraction'])
grass=imagereg.select(['grass-coverfraction'])
baresoil=imagereg.select(['bare-coverfraction'])
urban=imagereg.select(['urban-coverfraction'])
#discrete=imagereg.select(['discrete_classification'])
#palette = ee.List(discrete.get('discrete_classification_class_palette'))



task_config_tree = { 'description' : 'imageToDriveTree', 'folder' : 'GEE', \
               'fileNamePrefix' : 'TreeTest', 'scale' : 100, 'region' : geometry }

task_tree = ee.batch.Export.image.toDrive(tree, **task_config_tree)

task_tree.start()

task_config_crops = { 'description' : 'imageToDriveCrops', 'folder' : 'GEE', \
               'fileNamePrefix' : 'CropsTest', 'scale' : 100, 'region' : geometry }

task_crops = ee.batch.Export.image.toDrive(crops, **task_config_crops)

task_crops.start()

task_config_grass = { 'description' : 'imageToDriveGrass', 'folder' : 'GEE', \
               'fileNamePrefix' : 'GrassTest', 'scale' : 100, 'region' : geometry }

task_grass = ee.batch.Export.image.toDrive(grass, **task_config_grass)

task_grass.start()

task_config_baresoil = { 'description' : 'imageToDrivebaresoil', 'folder' : 'GEE', \
               'fileNamePrefix' : 'baresoilTest', 'scale' : 100, 'region' : geometry }

task_baresoil = ee.batch.Export.image.toDrive(baresoil, **task_config_baresoil)

task_baresoil.start()

task_config_urban = { 'description' : 'imageToDriveurban', 'folder' : 'GEE', \
               'fileNamePrefix' : 'urbanTest', 'scale' : 100, 'region' : geometry }

task_urban = ee.batch.Export.image.toDrive(urban, **task_config_urban)

task_urban.start()

"""
Derive precipitation 
"""
llx = 106.
lly = 36.
urx = 108.
ury = 38.
geometry = ee.Geometry.Rectangle([llx,lly,urx,ury])
geometry = geometry['coordinates'][0]
region=ee.Geometry.Rectangle([llx,lly,urx,ury])

'''
#test daily precipitation (working code)
inidate = ee.Date.fromYMD(2012,12,1)
enddate = ee.Date.fromYMD(2013,1,1)
#Import GSMaP data
gsmap = ee.ImageCollection('JAXA/GPM_L3/GSMaP/v6/reanalysis')\
    .filterDate(inidate, enddate).select('hourlyPrecipRateGC');
gsmap_reg=gsmap.filterBounds(region)
#lapse = ee.List.sequence(0, difdate.subtract(1))
#inidate = ee.Date('2014-01-01')
#days = ee.List.sequence(1, 365)
#Filter the collection in one day (24 images)
day=ee.Date('2012-12-31')
day_collection = gsmap_reg.filterDate(day, day.advance(1, 'day'))
sum_d = ee.Image(day_collection.sum())
#month_collection = gsmap_reg.filterDate(inidate,enddate)
#Get the sum of all 24 images into one Image
#sum_m = ee.Image(month_collection.sum())

task_config = { 'description' : 'ImagetoDriveM12D31', 'folder' : 'GEE', \
               'fileNamePrefix' : 'PrM12D31', 'dimensions' : 400, 'region' : geometry }

task = ee.batch.Export.image.toDrive(sum_d, **task_config)

task.start()


'''

'''
# The code below does not work

#test
inidate = ee.Date.fromYMD(2014,12,1)
enddate = ee.Date.fromYMD(2014,12,31)

#Difference between start and end in days 
difdate = enddate.difference(inidate, 'day')

def function(day):
    return inidate.advance(day, 'day')
listdates=lapse.map(function(day))

#Iterate over the list of dates
def function2(day):
    #cast
    day = ee.Date(day)

    #Filter the collection in one day (24 images)
    day_collection = gsmap.filterDate(day, day.advance(1, 'day'))

    #Get the sum of all 24 images into one Image
    sum = ee.Image(day_collection.sum())

    value = sum.reduceRegion(ee.Reducer.first(), geometry, 30).get('hourlyPrecipRateGC')
    
    return value
'''

months=np.array(['01','02','03','04','05','06','07','08','09','10','11','12'])
for m in range(len(months)): 
    if m==0 or m==2 or m==4 or m==6 or m==7 or m==9 or m==11:
        days=np.arange(1,32)
    elif m==1:
        days=np.arange(1,29)
    else:
        days=np.arange(1,30)
    #loop over days
    inidate = ee.Date.fromYMD(2012,m+1,1)
    if m<11:
        enddate = ee.Date.fromYMD(2012,m+2,1) 
    else:
        enddate = ee.Date.fromYMD(2013,1,1) 
    #Import GSMaP data
    gsmap = ee.ImageCollection('JAXA/GPM_L3/GSMaP/v6/reanalysis')\
        .filterDate(inidate, enddate).select('hourlyPrecipRateGC');
    gsmap_reg=gsmap.filterBounds(region)
    for d in days:
        if d<10:
            day=ee.Date('2012-%s-0%i' % (months[m],d))
            print(day)
        else:
            day=ee.Date('2012-%s-%i' % (months[m],d))
            print(day)
        day_collection = gsmap_reg.filterDate(day, day.advance(1 , 'day'))
        sum_d = ee.Image(day_collection.sum())
    
        task_config = { 'description' : 'ImagetoDriveM%sD%i' %(months[m],d), 'folder' : 'GEE', \
                       'fileNamePrefix' : 'PrM%sD%i' %(months[m],d)\
                           , 'dimensions' : 400, 'region' : geometry }

        task = ee.batch.Export.image.toDrive(sum_d, **task_config)

        task.start()
        
"""
Export terracing map
"""
llx = 106.
lly = 36.
urx = 108.
ury = 38.
geometry = ee.Geometry.Rectangle([llx,lly,urx,ury])
geometry = geometry['coordinates'][0]
region=ee.Geometry.Rectangle([llx,lly,urx,ury])

ChinaTerrace = ee.Image("users/cbw/TerraceResults/China_terrace_1");

wkt = 'PROJCS["Asia_Lambert_Conformal_Conic", \
    GEOGCS["GCS_WGS_1984", \
    DATUM["WGS_1984", \
    SPHEROID["WGS_1984",6378137,298.257223563]], \
    PRIMEM["Greenwich",0], \
    UNIT["Degree",0.017453292519943295]], \
    PROJECTION["Lambert_Conformal_Conic_2SP"], \
    PARAMETER["False_Easting",0], \
    PARAMETER["False_Northing",0],\
    PARAMETER["Central_Meridian",105],\
    PARAMETER["Standard_Parallel_1",30],\
    PARAMETER["Standard_Parallel_2",62],\
    PARAMETER["Latitude_Of_Origin",0],\
    UNIT["Meter",1], AUTHORITY["EPSG","102012"]]'

proj_Asia_Lambert = ee.Projection(wkt);
image_Asia_Lambert = ChinaTerrace.reproject(crs=proj_Asia_Lambert,scale= 30);

imagereg=image_Asia_Lambert.clip(region)

task_config = { 'description' : 'imageToDriveTer', 'folder' : 'GEE', \
               'fileNamePrefix' : 'Chinaterrace', 'scale' : 30, 'region' : geometry }
    
task_terrace = ee.batch.Export.image.toDrive(imagereg, **task_config)

task_elev.start()