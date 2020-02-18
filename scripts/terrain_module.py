"""
Terrain profile module.

Written by Ed Oughton and Tom Russell

June 2019

"""
import configparser
import glob

import os
import math
from functools import partial, lru_cache

import fiona
import pyproj
import rasterio
import numpy as np
from shapely.ops import transform
from shapely.geometry import Point, LineString, mapping
from rasterstats import zonal_stats
import configparser
from collections import OrderedDict

CONFIG = configparser.ConfigParser()
CONFIG.read(os.path.join(os.path.dirname(__file__), 'script_config.ini'))
BASE_PATH = CONFIG['file_locations']['base_path']


def terrain_area(dem_folder, transmitter, cell_range, current_crs):
    """
    This module takes a single set of point coordinates for a site
    along with an estimate of the cell range. The irregular terrain
    parameter is returned.

    Parameters
    ----------
    dem_folder : string
        Folder path to the available Digital Elevation Model tiles.
    transmitter : dict
        Geojson. Must be in WGS84 / EPSG: 4326
    cell_range : int
        Radius of cell area in meters.
    current_crs : string
        The coordinate reference system the transmitter coordinates
        are in (must be WGS84 / EPSG: 4326).

    Returns
    -------
    Inter-decile range : int
        The terrain irregularity parameter.

    """
    point_geometry = Point(transmitter['geometry']['coordinates'])

    wgs84_to_projected = partial(
        pyproj.transform,
        pyproj.Proj(init=current_crs),
        pyproj.Proj(init='EPSG:3857'))

    #convert line geometry to projected to allow buffer in meters
    point_geometry = transform(wgs84_to_projected, point_geometry)

    cell_area_projected = point_geometry.buffer(cell_range)

    projected_to_wgs84 = partial(
        pyproj.transform,
        pyproj.Proj(init='EPSG:3857'),
        pyproj.Proj(init=current_crs))

    cell_area_unprojected = transform(projected_to_wgs84, cell_area_projected)

    stats = zonal_stats([cell_area_unprojected],
                        os.path.join(BASE_PATH,'ASTGTM2_N51W001_dem.tif'),
                        add_stats={'interdecile_range':interdecile_range})

    return stats[0]['interdecile_range']


def interdecile_range(x):
    q90, q10 = np.percentile(x, [90, 10])
    return int(round(q90 - q10, 0))


def terrain_p2p(dem_folder, line, current_crs):
    """
    This module takes a set of point coordinates and returns
    the elevation profile.

    Line : dict
        Geojson. Must be in WGS84 / EPSG: 4326

    """
    extents = load_extents(dem_folder)

    line_geometry = LineString(line['geometry']['coordinates'])

    ll_to_osgb = partial(
        pyproj.transform,
        pyproj.Proj(init=current_crs),
        pyproj.Proj(init='EPSG:3857'))

    #convert line geometry to projected
    line_geometry = transform(ll_to_osgb, line_geometry)

    distance = int(line_geometry.length)

    increment = determine_distance_increment(distance)

    distance_km = distance / 1e3

    x = []
    y = []
    elevation_profile = []
    osgb_to_ll = partial(
        pyproj.transform,
        pyproj.Proj(init='EPSG:3857'),
        pyproj.Proj(init=current_crs))

    points = []

    for currentdistance  in range(0, distance, increment):
        point_old_crs = line_geometry.interpolate(currentdistance)
        point_new_crs = transform(osgb_to_ll, point_old_crs)
        xp, yp = point_new_crs.x, point_new_crs.y
        x.append(xp)
        y.append(yp)
        tile_path = get_tile_path_for_point(extents, xp, yp)
        z = get_value_from_dem_tile(tile_path, xp, yp)

        elevation_profile.append(z)

        points.append({
            'type': 'Feature',
            'geometry': mapping(point_old_crs),
            'properties': {
                'elevation': float(z),
                }
            })

    return elevation_profile, distance_km, points


def load_extents(dem_folder):
    """
    Check the extent of each DEM tile, save to dict for future reference.

    """
    extents = {}
    for tile_path in glob.glob(os.path.join(dem_folder, "*.tif")):
        dataset = rasterio.open(tile_path)
        # print("Extent of", tile_path, tuple(dataset.bounds))
        extents[tuple(dataset.bounds)] = tile_path
    return extents


def get_tile_path_for_point(extents, x, y):
    for (left, bottom, right, top), path in extents.items():
        if x >= left and x <= right and y <= top and y >= bottom:
            return path
    raise ValueError("No tile includes x {} y {}".format(x, y))


def get_value_from_dem_tile(tile_path, x, y):
    """
    Read all tile extents, load value from relevant tile.

    """
    dataset = rasterio.open(tile_path)
    row, col = dataset.index(x, y)
    band = dataset.read(1)
    dataset.close()
    return band[row, col]


def determine_distance_increment(distance):
    """
    Longley-Rice Irregular Terrain Model is limited to only 600
    elevation points, so this function ensures this number is not
    passed.

    """
    if distance >= 60000:
        return int(distance / 100)
    elif distance >= 30000:
        return int(distance / 50)
    elif distance >= 10000:
        return int(distance / 25)
    elif distance >= 1000:
        return int(distance / 10)
    else:
        return int(distance / 2)
