from osgeo import osr,gdal
from datetime import datetime
from dateutil import tz

import numpy as np
from numpy import *

from math import sin, cos, sqrt, atan2, radians

from flask import Flask
from flask import jsonify

app = Flask(__name__)

def GetLocalTime(utcdatetime, localtz):
	localdatetime = utcdatetime.replace(tzinfo=tz.UTC).astimezone(localtz)

	return str(localdatetime.isoformat())

@app.route("/vegetation-cover")
def vegetationcover():
	# TODO: aceitar como parametro
	filename = "319567_2331703_2016-12-07_0c0b-20161207T151953Z.tif"

	filenameinfo = GetFileNameInfo(filename);

	# aqui assumimos que o indice ndvi acima de 0.4 se considera vegetacao verde
	# esse indice pode ser alterado dependendo do nivel de exatidao necessaria
	coverageIdx = 0.4
	vegCover = GetVegetationCover(filename, coverageIdx)

	area = GetArea(filename)
	centroid = GetCentroid(filename)

	localtz = tz.gettz("America/Sao_Paulo")

	# nao consegui achar um metadado no arquivo que me mostra o timestamp/timezoneId
	# portanto, fiz o calculo manualmente
	capdatetimeutc = datetime(2016, 12, 7, 15, 19, 53, 000, tzinfo=tz.UTC)

	localTime = GetLocalTime(capdatetimeutc, localtz)

	result = {
		"filename": filenameinfo,
		"cover": vegCover,
		"area": area,
		"centroid":centroid,
		"local_time": localTime
	}

	return jsonify(result)

def GetFileNameInfo(filename):
	dataset = gdal.Open(filename)

	return dataset.GetDescription()

def GetVegetationCover(filename, coverageIdx):
	dataset = gdal.Open(filename)

	red = dataset.GetRasterBand(3).ReadAsArray();
	nir = dataset.GetRasterBand(4).ReadAsArray();

	red = red.astype(np.float64)
	nir = nir.astype(np.float64)

	ndvi = CalculateNDVI(red, nir)

	avg = 0;
	count = 0;
	for x in ndvi:
		for i in x:
			count = count + 1
			if i >= coverageIdx:
				avg = avg + 1

	perc = float((avg * 100) / count) / 100;

	dataset = None;

	return round(perc, 1)

def GetCentroid(filename):
	dataset = gdal.Open(filename)
	old_cs= osr.SpatialReference()
	old_cs.ImportFromWkt(dataset.GetProjectionRef())

	# create the new coordinate system
	wgs84_wkt = """
	GEOGCS["WGS 84",
		DATUM["WGS_1984",
			SPHEROID["WGS 84",6378137,298.257223563,
				AUTHORITY["EPSG","7030"]],
			AUTHORITY["EPSG","6326"]],
		PRIMEM["Greenwich",0,
			AUTHORITY["EPSG","8901"]],
		UNIT["degree",0.01745329251994328,
			AUTHORITY["EPSG","9122"]],
		AUTHORITY["EPSG","4326"]]"""
	new_cs = osr.SpatialReference()
	new_cs .ImportFromWkt(wgs84_wkt)

	# create a transform object to convert between coordinate systems
	transform = osr.CoordinateTransformation(old_cs,new_cs)

	# get the point to transform, pixel (0,0) in this case
	width = dataset.RasterXSize
	height = dataset.RasterYSize
	gt = dataset.GetGeoTransform()
	minx = gt[0]
	miny = gt[3] + width*gt[4] + height*gt[5]
	maxx = gt[0] + width*gt[1] + height*gt[2]
	maxy = gt[3]

	mid = MidpointEuclidean(minx, miny, maxx, maxy)

	centroid = transform.TransformPoint(mid[0],mid[1])

	dataset = None;

	return {
		"type": "Point",
    	"coordinates": [centroid[1], centroid[0]]
	}

def GetArea(filename):
	# approximate radius of earth in km
	R = 6371.0

	dataset = gdal.Open(filename)
	old_cs= osr.SpatialReference()
	old_cs.ImportFromWkt(dataset.GetProjectionRef())

	# create the new coordinate system
	wgs84_wkt = """
	GEOGCS["WGS 84",
	    DATUM["WGS_1984",
	        SPHEROID["WGS 84",6378137,298.257223563,
	            AUTHORITY["EPSG","7030"]],
	        AUTHORITY["EPSG","6326"]],
	    PRIMEM["Greenwich",0,
	        AUTHORITY["EPSG","8901"]],
	    UNIT["degree",0.01745329251994328,
	        AUTHORITY["EPSG","9122"]],
	    AUTHORITY["EPSG","4326"]]"""
	new_cs = osr.SpatialReference()
	new_cs .ImportFromWkt(wgs84_wkt)

	# create a transform object to convert between coordinate systems
	transform = osr.CoordinateTransformation(old_cs,new_cs)

	# get the point to transform, pixel (0,0) in this case
	width = dataset.RasterXSize
	height = dataset.RasterYSize
	gt = dataset.GetGeoTransform()
	minx = gt[0]
	miny = gt[3] + width*gt[4] + height*gt[5]
	maxx = gt[0] + width*gt[1] + height*gt[2]
	maxy = gt[3]

	minbase = transform.TransformPoint(minx,miny)
	maxbase = transform.TransformPoint(maxx,miny)

	lat1 = radians(minbase[1])
	lon1 = radians(minbase[0])
	lat2 = radians(maxbase[1])
	lon2 = radians(maxbase[0])

	dlon = lon2 - lon1
	dlat = lat2 - lat1

	a = sin(dlat / 2)**2 + cos(lat1) * cos(lat2) * sin(dlon / 2)**2
	c = 2 * atan2(sqrt(a), sqrt(1 - a))

	baseDistance = R * c

	minheight = transform.TransformPoint(maxx,miny)
	maxheight = transform.TransformPoint(maxx,maxy)


	lat1 = radians(minheight[1])
	lon1 = radians(minheight[0])
	lat2 = radians(maxheight[1])
	lon2 = radians(maxheight[0])

	dlon = lon2 - lon1
	dlat = lat2 - lat1

	a = sin(dlat / 2)**2 + cos(lat1) * cos(lat2) * sin(dlon / 2)**2
	c = 2 * atan2(sqrt(a), sqrt(1 - a))

	heightDistance = R * c
	area = float(baseDistance) * float(heightDistance)

	dataset = None;

	return "%.2f" % area

def MidpointEuclidean(x1,y1,x2,y2):
	dist_x = abs(x1-x2) / 2.
	dist_y = abs(y1-y2) / 2.
	res_x = x1 - dist_x if x1 > x2 else x2 - dist_x
	res_y = y1 - dist_y if y1 > y2 else y2 - dist_y

	return res_x, res_y

def CalculateNDVI(red, nir):
    return (nir - red) / (nir + red)

if __name__ == "__main__":
    app.run()
