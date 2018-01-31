# Copyright (c) 2018, Vienna University of Technology (TU Wien), Department of
# Geodesy and Geoinformation (GEO).
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#
# The views and conclusions contained in the software and documentation are
# those of the authors and should not be interpreted as representing official
# policies, either expressed or implied, of the FreeBSD Project.


"""
Code for converting raster images to TPS grids.
"""

import os
import errno

from osgeo import osr

import pytileproj.gdalport as gdalport
import pytileproj.geometry as geometry


def warp2tiledgeotiff(TPS, image, output_dir,
                      gdal_path=None, subgrid_ids=None,
                      accurate_boundary=True, e7_folder=True, ftiles=None,
                      roi=None, coverland=True, outshortname=None,
                      withtilenameprefix=False, withtilenamesuffix=True,
                      compress=True, compresstype="LZW", resampling_type="near",
                      overwrite=False, image_nodata=None, tile_nodata=None,
                      tiledtiff=True, blocksize=512):
    """warp the image to tiles in target TPS.


    Parameters
    ----------
    grid : TiledProjectionSystem
        a compatible TPS object, e.g. a Equi7Grid(500) object.
    image : string
        Image file path.
    ouput_dir : string
        Output directory path.
    e7_folder: boolean
        if set (default), the output data will be stored in equi7 folder structure
    gdal_path : string
        Gdal utilities location.
    subgrid_ids : list
        Only resample to the specified continents,
        default is to resample to all continents.
    roi : OGRGeometry
        Region of interest.
        The roi will beignored if ftiles keyword is provided
    ftiles : list of tile names
        full name of tiles to which data should be resampled
    coverland : bool
        Only resample to the tiles that cover land.
    outshortname : string
        The short name will be main part of the output tile name.
    withtilenameprefix : bool
        Prepend the tile name in the tile file name.
    withtilenamesuffix : bool
        Append the tile name in the tile file name.
    compress : bool
        The output tiles is compressed or not.
    resampling_type : string
        Resampling method.
    overwrite : bool
        Overwrite the old tile or not.
    image_nodata : double
        The nodata value of input images.
    tile_nodata : double
        The nodata value of tile.
    tiledtiff : bool
        Set to yes for tiled geotiff output
    blocksize: integer
        sets tile width, in x and y direction

    """
    # get the output shortname if not provided

    grid = TPS

    if not outshortname:
        outshortname = os.path.splitext(os.path.basename(image))[0]

    # find overlapping tiles
    if ftiles is None:
        if roi is None:
            if accurate_boundary:
                try:
                    roi_geom = gdalport.retrieve_raster_boundary(image,
                                                                 gdal_path=gdal_path,
                                                                 nodata=image_nodata)
                except Exception as e:
                    print("retrieve_raster_boundary failed:", str(e))
                    roi_geom = None
            else:
                roi_geom = None
            if roi_geom:
                ftiles = grid.search_tiles_in_geo_roi(geom_area=roi_geom,
                                                      subgrid_ids=subgrid_ids,
                                                      coverland=coverland,
                                                      gdal_path=gdal_path)
            else:
                ds = gdalport.open_image(image)
                img_extent = ds.get_extent()
                geo_extent = geometry.extent2polygon(img_extent)
                geo_sr = osr.SpatialReference()
                geo_sr.ImportFromWkt(ds.projection())
                geo_extent.AssignSpatialReference(geo_sr)
                ftiles = grid.search_tiles_in_geo_roi(geom_area=geo_extent,
                                                      subgrid_ids=subgrid_ids,
                                                      coverland=coverland,
                                                      gdal_path=gdal_path)
        else:
            ftiles = grid.search_tiles_in_geo_roi(geom_area=roi,
                                                  subgrid_ids=subgrid_ids,
                                                  coverland=coverland,
                                                  gdal_path=gdal_path)
    else:
        if type(ftiles) != list:
            ftiles = [ftiles]

    # keep the full path of the output(resampled) files
    dst_file_names = []

    # resample for each tile sequentially
    for ftile in ftiles:
        # create grid folder
        if e7_folder:
            grid_folder = "EQUI7_{}".format(ftile[0:6])
            tile_path = os.path.join(output_dir, grid_folder, ftile[7:])
            if not os.path.exists(tile_path):
                makedirs(tile_path)
        else:
            tile_path = output_dir

        # make output filename
        outbasename = outshortname
        if withtilenameprefix:
            outbasename = "_".join((ftile, outbasename))
        if withtilenamesuffix:
            outbasename = "_".join((outbasename, ftile))
        filename = os.path.join(tile_path, "".join((outbasename, ".tif")))

        # get grid tile object
        # TODO Generalise! Now only the Equi7Grid case is consideres (2 digits
        # for subgrid)
        gridtile = getattr(grid, ftile[0:2]).tilesys.create_tile(ftile)

        # prepare options for gdalwarp
        options = {'-t_srs': gridtile.core.projection.wkt, '-of': 'GTiff',
                   '-r': resampling_type,
                   '-te': " ".join(map(str, gridtile.limits_m())),
                   '-tr': "{} -{}".format(grid.core.res, grid.core.res)}

        options["-co"] = list()
        if compress:
            options["-co"].append("COMPRESS={0}".format(compresstype))
        if image_nodata != None:
            options["-srcnodata"] = image_nodata
        if tile_nodata != None:
            options["-dstnodata"] = tile_nodata
        if overwrite:
            options["-overwrite"] = " "
        if tiledtiff:
            options["-co"].append("TILED=YES")
            options["-co"].append("BLOCKXSIZE={0}".format(blocksize))
            options["-co"].append("BLOCKYSIZE={0}".format(blocksize))

        # call gdalwarp for resampling
        succeed, _ = gdalport.call_gdal_util('gdalwarp', src_files=image,
                                             dst_file=filename, gdal_path=gdal_path,
                                             options=options)

        # full path to the output(resampled) files
        if succeed:
            dst_file_names.extend([filename])

    return dst_file_names


def makedirs(dir_path):
    """ make directory safely for parallelism
    """
    try:
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise


def mkdir(name):
    """ make directory safely for parallelism
    """
    try:
        if not os.path.exists(name):
            os.mkdir(name)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise
