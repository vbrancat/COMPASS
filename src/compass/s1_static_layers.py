import re
import time

import isce3
import journal

from osgeo import gdal
from compass.utils.geo_runconfig import GeoRunConfig
from compass.utils.helpers import (bursts_grouping_generator, get_module_name,
                                   get_time_delta_str)
from compass.utils.yaml_argparse import YamlArgparse

file_name_los_east = 'los_east'
file_name_los_north = 'los_north'
file_name_local_incidence = 'local_incidence_angle'
file_name_layover = 'layover_shadow_mask'
file_name_x = 'x'
file_name_y = 'y'
file_name_z = 'z'


def _make_rdr2geo_cfg(yaml_runconfig_str):
    '''
    Make a rdr2geo specific runconfig with latitude, longitude, and height
    layers enabled for static layer product generation while preserving all
    other rdr2geo config settings

    Parameters
    ----------
    yaml_runconfig_str: str
        Workflow runconfig as a string

    Returns
    -------
    rdr2geo_cfg: dict
        Dictionary with rdr2geo longitude, latitude, and height layers
        enabled. All other rdr2geo parameters are from *yaml_runconfig_str*
    '''
    # If any of the requisite layers are false, make them true in yaml cfg str
    for layer in ['latitude', 'longitude', 'incidence_angle']:
        re.sub(f'compute_{layer}:\s+[Ff]alse', f'compute_{layer}: true',
               yaml_runconfig_str)

    # Load a GeoRunConfig from modified yaml cfg string
    rdr2geo_cfg = GeoRunConfig.load_from_yaml(yaml_runconfig_str,
                                              workflow_name='s1_cslc_geo')

    return rdr2geo_cfg


def run(cfg: GeoRunConfig):
    '''
    Run static layers workflow (i.e., generate static layers,
    geocode them, create product HDF5) with user-defined
    args stored in dictionary runconfig *cfg*

    Parameters
    ---------
    cfg: GeoRunConfig
        GeoRunConfig object with user runconfig options
    '''

    module_name = get_module_name(__file__)
    info_channel = journal.info(f"{module_name}.run")
    info_channel.log(f"Starting {module_name} burst")

    # Extract which layers to generate
    rdr2geo_cfg = cfg.rdr2geo_params

    # Start tracking processing time
    t_start = time.perf_counter()

    for burst_id, bursts in bursts_grouping_generator(cfg.bursts):
        burst = bursts[0]

        date_str = burst.sensing_start.strftime("%Y%m%d")

        info_channel.log(f'Starting geocoding of {burst_id} for {date_str}')

        # Generate required static layers
        # rdr2geo_cfg = _make_rdr2geo_cfg(cfg.yaml_string)
        # s1_rdr2geo.run(rdr2geo_cfg, burst, save_in_scratch=True)
        # s1_geocode_metadata.run(cfg, burst, fetch_from_scratch=True)

        # Generate static layers using Gustavo's get_radar_grid
        geogrid = cfg.geogrids[burst_id]
        dem_raster = isce3.io.Raster(cfg.dem)
        radar_grid = burst.as_isce3_radargrid()
        native_doppler = burst.doppler.lut2d
        grid_doppler = isce3.core.LUT2d()

        date_str = burst.sensing_start.strftime("%Y%m%d")
        burst_id = str(burst.burst_id)

        # init output directory in product_path
        burst_id_date_key = (burst_id, date_str)

        out_paths = cfg.output_paths[burst_id_date_key]
        output_path = out_paths.output_directory
        output_path = out_paths.scratch_directory

        topo_output = {file_name_x: (rdr2geo_cfg.compute_longitude, gdal.GDT_Float64),
                       file_name_y: (rdr2geo_cfg.compute_latitude, gdal.GDT_Float64),
                       file_name_z: (rdr2geo_cfg.compute_height, gdal.GDT_Float64),
                       file_name_layover: (
                       cfg.rdr2geo_params.compute_layover_shadow_mask,
                       gdal.GDT_Byte),
                       file_name_local_incidence: (
                       rdr2geo_cfg.compute_local_incidence_angle,
                       gdal.GDT_Float32),
                       file_name_los_east: (
                       rdr2geo_cfg.compute_ground_to_sat_east, gdal.GDT_Float32),
                       file_name_los_north: (
                       rdr2geo_cfg.compute_ground_to_sat_north, gdal.GDT_Float32),
                       }
        raster_list = [
            isce3.io.Raster(f'{output_path}/{fname}.geo', geogrid.width,
                            geogrid.length, 1, dtype, 'ENVI')
            if enabled else None
            for fname, (enabled, dtype) in topo_output.items()]

        (x_raster, y_raster, z_raster, layover_shadow_raster,
         local_incident_angle_raster, los_east_raster,
         los_north_raster) = raster_list

        isce3.geogrid.get_radar_grid(radar_grid.lookside, burst.wavelength,
                                     dem_raster, geogrid, burst.orbit,
                                     native_doppler, grid_doppler,
                                     local_incidence_angle_raster=local_incident_angle_raster,
                                     lo_unit_vector_x_raster=los_east_raster,
                                     los_unit_vector_y_raster=los_north_raster)

    dt = get_time_delta_str(t_start)
    info_channel.log(f"{module_name} burst successfully ran in {dt} (hr:min:sec)")


def main():
    """Create the CLI and run the static layers workflow"""
    # load arguments from command line
    parser = YamlArgparse()

    # Get a runconfig dict from command line arguments
    cfg = GeoRunConfig.load_from_yaml(parser.run_config_path,
                                      workflow_name='s1_cslc_geo')

    run(cfg)

if __name__ == "__main__":
    main()
