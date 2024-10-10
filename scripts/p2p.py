"""
Point to Point prediction mode runner.

Referred to as qkpfl in the original Fortran codebase.

Written by Ed Oughton

June 2019

"""
import configparser
import json
import os
import csv
import math
import pickle
import time

from create_random import *

import numpy as np
from functools import partial
from collections import OrderedDict

import fiona
import rasterio
from fiona.crs import from_epsg
from pyproj import Transformer
from rasterio.warp import Resampling
from rasterio.warp import calculate_default_transform, reproject
from shapely import Point
from shapely.geometry import LineString, mapping
from shapely.ops import transform

from itmlogic.misc.qerfi import qerfi
from itmlogic.preparatory_subroutines.qlrpfl import qlrpfl
from itmlogic.statistics.avar import avar
from terrain_module import terrain_p2p, terrain_p2p_wgs, terrain_p2p_linear

# #set up file paths
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_PATH = os.path.abspath(os.path.join(SCRIPT_DIR, '..'))
DATA_FOLDER = os.path.join(BASE_PATH, 'data')
DATA_PROCESSED = os.path.join(DATA_FOLDER, 'processed')

RESULTS = os.path.join(BASE_PATH, 'results')

DEM_FOLDER = os.path.join(DATA_FOLDER)
DIRECTORY_SHAPES = os.path.join(DATA_PROCESSED, 'shapes')



def reproject_crs(input_file, target_crs):
    # Open the source file (in WGS84, EPSG:4326)
    with rasterio.open(input_file) as src:
    # with rasterio.open(input_file) as src:
        # Calculate the transform and dimensions for the new projection
        transform, width, height = calculate_default_transform(
            src.crs, target_crs, src.width, src.height, *src.bounds
        )

        # Define metadata for the reprojected file
        kwargs = src.meta.copy()
        kwargs.update({
            'crs': target_crs,
            'transform': transform,
            'width': width,
            'height': height
        })

        # Create a new array to hold the reprojected data
        data_reprojected = np.empty((height, width), dtype=src.meta['dtype'])

        # Reproject the data
        reproject(
            source=rasterio.band(src, 1),
            destination=data_reprojected,
            src_transform=src.transform,
            src_crs=src.crs,
            dst_transform=transform,
            dst_crs=target_crs,
            resampling=Resampling.nearest
        )
        return data_reprojected, transform


def itmlogic_p2p(main_user_defined_parameters, surface_profile_m):
    """
    Run itmlogic in point to point (p2p) prediction mode.

    Parameters
    ----------
    main_user_defined_parameters : dict
        User defined parameters.
    surface_profile_m : list
        Contains surface profile measurements in meters.

    Returns
    -------
    output : list of dicts
        Contains model output results.

    """
    prop = main_user_defined_parameters

    #DEFINE ENVIRONMENTAL PARAMETERS
    # Terrain relative permittivity
    prop['eps'] = 15

    # Terrain conductivity (S/m)
    prop['sgm']   = 0.005

    # Climate selection (1=equatorial,
    # 2=continental subtropical, 3=maritime subtropical,
    # 4=desert, 5=continental temperate,
    # 6=maritime temperate overland,
    # 7=maritime temperate, oversea (5 is the default)
    prop['klim']  =   5

    # Surface refractivity (N-units): also controls effective Earth radius
    prop['ens0']  =   314

    #DEFINE STATISTICAL PARAMETERS
    # Confidence  levels for predictions
    # qc = [50, 90, 10]
    qc = [50]

    # Reliability levels for predictions
    # qr = [1, 10, 50, 90, 99]
    qr = [50]

    # Number of points describing profile -1
    pfl = []
    pfl.append(len(surface_profile_m) - 1)
    pfl.append(0)

    for profile in surface_profile_m:
        pfl.append(profile)

    # Refractivity scaling ens=ens0*exp(-zsys/9460.)
    # (Average system elev above sea level)
    zsys = 0

    # Note also defaults to a continental temperate climate

    # Setup some intermediate quantities
    # Initial values for AVAR control parameter: LVAR=0 for quantile change,
    # 1 for dist change, 2 for HE change, 3 for WN change, 4 for MDVAR change,
    # 5 for KLIM change
    prop['lvar'] = 5

    # Inverse Earth radius
    prop['gma']  = 157E-9

    # Conversion factor to db
    db = 8.685890

    #Number of confidence intervals requested
    nc = len(qc)

    #Number of reliability intervals requested
    nr = len(qr)

    #Length of profile in km
    dkm = prop['d']

    #Profile range step, select option here to define range step from profile
    #length and # of points
    xkm = 0

    #If DKM set <=0, find DKM by mutiplying the profile step by number of
    #points (not used here)
    if dkm <= 0:
        dkm = xkm * pfl[0]

    #If XKM is <=0, define range step by taking the profile length/number
    #of points in profile
    if xkm <= 0:

        xkm = dkm // pfl[0]

        #Range step in meters stored in PFL(2)
        pfl[1] = dkm * 1000 / pfl[0]

        #Store profile in prop variable
        prop['pfl'] = pfl
        #Zero out error flag
        prop['kwx'] = 0
        #Initialize omega_n quantity
        prop['wn'] = prop['fmhz'] / 47.7
        #Initialize refractive index properties
        prop['ens'] = prop['ens0']

    #Scale this appropriately if zsys set by user
    if zsys != 0:
        prop['ens'] = prop['ens'] * math.exp(-zsys / 9460)

    #Include refraction in the effective Earth curvature parameter
    prop['gme'] = prop['gma'] * (1 - 0.04665 * math.exp(prop['ens'] / 179.3))

    #Set surface impedance Zq parameter
    zq = complex(prop['eps'], 376.62 * prop['sgm'] / prop['wn'])

    #Set Z parameter (h pol)
    prop['zgnd'] = np.sqrt(zq - 1)

    #Set Z parameter (v pol)
    if prop['ipol'] != 0:
        prop['zgnd'] = prop['zgnd'] / zq

    #Flag to tell qlrpfl to set prop.klim=prop.klimx and set lvar to initialize avar routine
    prop['klimx'] = 0

    #Flag to tell qlrpfl to use prop.mdvar=prop.mdvarx and set lvar to initialize avar routine
    prop['mdvarx'] = 11

    #Convert requested reliability levels into arguments of standard normal distribution
    zr = qerfi([x / 100 for x in qr])
    #Convert requested confidence levels into arguments of standard normal distribution
    zc = qerfi([x / 100 for x in qc])

    #Initialization routine for point-to-point mode that sets additional parameters
    #of prop structure
    prop = qlrpfl(prop)

    ## Here HE = effective antenna heights, DL = horizon distances,
    ## THE = horizon elevation angles
    ## MDVAR = mode of variability calculation: 0=single message mode,
    ## 1=accidental mode, 2=mobile mode, 3 =broadcast mode, +10 =point-to-point,
    ## +20=interference

    #Free space loss in db
    fs = db * np.log(2 * prop['wn'] * prop['dist'])

    #Used to classify path based on comparison of current distance to computed
    #line-of-site distance
    q = prop['dist'] - prop['dlsa']

    #Scaling used for this classification
    q = max(q - 0.5 * pfl[1], 0) - max(-q - 0.5 * pfl[1], 0)

    #Report dominant propagation type predicted by model according to parameters
    #obtained from qlrpfl
    # if q < 0:
    #     print('Line of sight path')
    # elif q == 0:
    #     print('Single horizon path')
    # else:
    #     print('Double-horizon path')
    # if prop['dist'] <= prop['dlsa']:
    #     print('Diffraction is the dominant mode')
    # elif prop['dist'] > prop['dx']:
    #     print('Tropospheric scatter is the dominant mode')
    #
    # print('Estimated quantiles of basic transmission loss (db)')
    # print('Free space value {} db'.format(str(fs)))

    # print('Confidence levels {}, {}, {}'.format(
    #     str(qc[0]), str(qc[1]), str(qc[2])))

    # # Confidence  levels for predictions
    # qc = [50, 90, 10]
    #
    # # Reliability levels for predictions
    # qr = [1, 10, 50, 90, 99]

    output = []
    for jr in range(0, (nr)):
        for jc in range(0, nc):
            #Compute corrections to free space loss based on requested confidence
            #and reliability quantities
            avar1, prop = avar(zr[jr], 0, zc[jc], prop)
            output.append({
                'distance_km': prop['d'],
                'reliability_level_%': qr[jr],
                'confidence_level_%': qc[jc],
                'propagation_loss_dB': fs + avar1 #Add free space loss and correction
                })

    return output


def csv_writer(data, directory, filename):
    """
    Write data to a CSV file path.

    Parameters
    ----------
    data : list of dicts
        Data to be written.
    directory : string
        Folder to write the results to.
    filename : string
        Name of the file to write.

    """
    # Create path
    if not os.path.exists(directory):
        os.makedirs(directory)

    fieldnames = []
    for name, value in data[0].items():
        fieldnames.append(name)

    with open(os.path.join(directory, filename), 'w') as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames, lineterminator = '\n')
        writer.writeheader()
        writer.writerows(data)


def straight_line_from_points(a, b):
    """
    Generate a geojson LineString object from two geojson points.

    Parameters
    ----------
    a : geojson
        Point A
    b : geojson
        Point B

    Returns
    -------
    line : geojson
        A geojson LineString object.

    """
    line = {
        'type': 'Feature',
        'geometry': {
            'type': 'LineString',
            'coordinates': [
                (
                    a['geometry']['coordinates'][0],
                    a['geometry']['coordinates'][1]
                ),
                (
                    b['geometry']['coordinates'][0],
                    b['geometry']['coordinates'][1]
                ),
            ]
        },
        'properties': {
            'id': 'terrain path'
        }
    }

    return line


if __name__ == '__main__':
    #Set coordinate reference systems
    old_crs = 'EPSG:4326'

    main_user_defined_parameters = {}
    main_user_defined_parameters['fmhz']  =  12450.0
    main_user_defined_parameters['hg'] = [143.9, 8.5]
    main_user_defined_parameters['ipol'] = 0

    #Create new geojson for Crystal Palace radio transmitter
    transmitter = {
        'type': 'Feature',
        'geometry': {
            'type': 'Point',
            'coordinates': (-0.07491679518573545, 51.42413477117786)
        },
        'properties': {
            'id': 'Crystal Palace radio transmitter'
        }
    }

    transmitter_point = Point(-0.07491679518573545, 51.42413477117786)
    transmitter_point_linear = Point(533941.9879948243, 171216.04150133877)

    # Load the receiver data from the JSON file
    with open('receivers.json', 'r') as f:
        receivers_json = json.load(f)

    with open("receivers.pkl", 'rb') as f:
        receivers = pickle.load(f)

    print(receivers.json[0])
    print(receivers.wgs[0])
    print(receivers.linear[0])

    dem_string = os.path.join(DEM_FOLDER, 'Copernicus_DSM_30_N51_00_W001_00_DEM.tif')
    dem_wgs, transform_wgs = reproject_crs(dem_string, 'EPSG:4326')
    dem_linear, transform_linear = reproject_crs(dem_string, 'EPSG:27700')

    output = []
    # Iterate over each receiver
    for i in range(10):
        line = straight_line_from_points(transmitter, receivers.json[i])

        #Run terrain modules
        start_time = time.time()
        measured_terrain_profile, distance_km, _ = terrain_p2p(dem_string, line)
        print(f"terrain_old took {time.time() - start_time} seconds")

        start_time = time.time()
        measured_terrain_profile_wgs, distance_km_wgs, _ = terrain_p2p_wgs(dem_wgs, transform_wgs, transmitter_point, receivers.wgs[i])
        print(f"terrain_new took {time.time() - start_time} seconds")

        start_time = time.time()
        measured_terrain_profile_linear, distance_km_linear, _ = terrain_p2p_linear(dem_linear, transform_linear, transmitter_point_linear, receivers.linear[i])
        print(f"terrain_linear took {time.time() - start_time} seconds")

        main_user_defined_parameters['d'] = distance_km

        #Run model and get output
        output.extend(itmlogic_p2p(main_user_defined_parameters, measured_terrain_profile))
        output.extend(itmlogic_p2p(main_user_defined_parameters, measured_terrain_profile_wgs))

        start_time = time.time()
        output.extend(itmlogic_p2p(main_user_defined_parameters, measured_terrain_profile_linear))
        print(f"itmlogic_p2p took {time.time() - start_time} seconds")
        output.append({})

    #Write results to .csv
    csv_writer(output, RESULTS, 'p2p_results.csv')


    print(f'Completed run in {time.time() - start_time}')
