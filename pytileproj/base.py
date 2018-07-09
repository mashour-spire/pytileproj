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
Code for Tiled Projection Systems.
"""

import abc
import itertools

import numpy as np
from osgeo import osr

import pytileproj.geometry as geometry
import pyproj


class TPSCoreProperty(object):

    """
    Class holding information needed at every level of `TiledProjectionSystem`,
    the alltime-valid "core properties".
    With this, core parameters are everywhere accessible via the same name.
    """

    def __init__(self, tag, projection, sampling, tiletype,
                 tile_xsize_m, tile_ysize_m):
        """
        Initialises a TPSCoreProperty.

        Parameters
        ----------
        tag : str
            identifier of the object holding the TPSCoreProperty
            e.g. 'EU' or 'Equi7'.
        projection : None or TPSProjection
            A TPSProjection() holding info on the spatial reference
        sampling : int
            the grid sampling = size of pixels; in metres.
        tiletype : str
            tilecode (related the tile size of the grid)
        tile_xsize_m : int
            tile size in x direction defined for the grid's sampling
        tile_ysize_m :
            tile size in y direction defined for the grid's sampling
        """

        self.tag = tag
        self.projection = projection
        self.sampling = sampling
        self.tiletype = tiletype
        self.tile_xsize_m = tile_xsize_m
        self.tile_ysize_m = tile_ysize_m


class TPSProjection():

    """
    Projection class holding and translating the definitions of a projection when initialising.
    """

    def __init__(self, epsg=None, proj4=None, wkt=None):
        """
        Initialises a TPSProjection().

        Parameters
        ----------
        epsg : inte
            The EPSG-code of the spatial reference.
            As from http://www.epsg-registry.org
            Not all reference do have a EPSG code.
        proj4 : str
            The proj4-string defining the spatial reference.
        wkt : str
            The wkt-string (well-know-text) defining the spatial reference.

        Notes
        -----
        Either one of epsg, proj4, or wkt must be given.
        """

        checker = {epsg, proj4, wkt}
        checker.discard(None)
        if len(checker) == 0:
            raise ValueError('Projection is not defined!')

        if len(checker) != 1:
            raise ValueError('Projection is defined ambiguously!')

        spref = osr.SpatialReference()

        if epsg is not None:
            spref.ImportFromEPSG(epsg)
            self.osr_spref = spref
            self.proj4 = spref.ExportToProj4()
            self.wkt = spref.ExportToWkt()
            self.epsg = epsg

        if proj4 is not None:
            spref.ImportFromProj4(proj4)
            self.osr_spref = spref
            self.proj4 = proj4
            self.wkt = spref.ExportToWkt()
            self.epsg = self.extract_epsg(self.wkt)

        if wkt is not None:
            spref.ImportFromWkt(wkt)
            self.osr_spref = spref
            self.proj4 = spref.ExportToProj4()
            self.wkt = wkt
            self.epsg = self.extract_epsg(self.wkt)


    def extract_epsg(self, wkt):
        """
        Checks if the WKT contains an EPSG code for the spatial reference,
        and returns it, if found.

        Parameters
        ----------
        wkt : string
            The wkt-string (well-know-text) defining the spatial reference.

        Returns
        -------
        epsg : integer, None
            the EPSG code of the spatial reference (if found). Else: None
        """

        pos_last_code = wkt.rfind('EPSG')
        pos_end = len(wkt)
        if pos_end - pos_last_code < 16:
            epsg = int(wkt[pos_last_code + 7:pos_last_code + 11])
        else:
            epsg = None

        return epsg


class TiledProjectionSystem(object):

    __metaclass__ = abc.ABCMeta

    # placeholders for static data defining the grid
    # static attribute
    _static_data = None
    # sub grid IDs
    _static_subgrid_ids = ['SG']
    # supported tile widths (linked to grid sampling)
    _static_tilecodes = ['T1']
    # supported grid spacing ( = the pixel sampling)
    _static_sampling = [1]

    def __init__(self, sampling, tag='TPS'):
        """
        Initialises a TiledProjectionSystem().

        Parameters
        ----------
        sampling : int
            the grid sampling = size of pixels; in metres.
        tag : str
            identifier of the object holding the TPSCoreProperty
            e.g. 'EU' or 'Equi7'.
        """

        tiletype = self.get_tiletype(sampling)
        tile_xsize_m, tile_ysize_m = self.get_tilesize(sampling)

        self.core = TPSCoreProperty(
            tag, None, sampling, tiletype, tile_xsize_m, tile_ysize_m)

        self.subgrids = self.define_subgrids()

    def __getattr__(self, item):
        '''
        short link for items of subgrids and core
        '''
        if item in self.subgrids:
            return self.subgrids[item]
        elif item in self.core.__dict__:
            return self.core.__dict__[item]
        else:
            return self.__dict__[item]


    @abc.abstractmethod
    def define_subgrids(self):
        pass


    def locate_geometry_in_subgrids(self, geom):
        """
        finds overlapping subgrids of given geometry.

        Attributes
        ----------
        geom : OGRGeometry
            a geometry to be located

        Returns
        -------
        list of TiledProjection()
            all subgrids that overlap with geom
        """

        covering_subgrid = list()
        for x in self.subgrids.keys():
            geometry.intersect_geometry(geom, self.subgrids.get(x).polygon_geog)
            if geom.Intersects(self.subgrids.get(x).polygon_geog):
                covering_subgrid.append(x)

        return covering_subgrid


    def lonlat2xy(self, lon, lat, subgrid=None):
        """
        converts latitude and longitude coordinates to TPS grid coordinates

        Parameters
        ----------
        lon : list of numbers
            longitude coordinates
        lat : list of numbers
            latitude coordinates
        subgrid : str
            optional: acronym / subgrid ID to search within (speeding up)

        Returns
        -------
        subgrid : str
            subgrid ID
        x, y : list of float
            TPS grid coordinates
        """

        if subgrid is None:
            vfunc = np.vectorize(self._lonlat2xy)
            return vfunc(lon, lat)
        else:
            return self._lonlat2xy_subgrid(lon, lat, subgrid)


    def _lonlat2xy(self, lon, lat):
        """
        finds overlapping subgrids of a given point in lon-lat-space
        and computes the projected coordinates.

        Parameters
        ----------
        lon : number
            longitude coordinate
        lat : number
            latitude coordinate

        Returns
        -------
        subgrid : str
            subgrid ID
        x, y : float
            TPS grid coordinates
        """

        # create point geometry
        lonlatprojection = TPSProjection(epsg=4326)
        point_geom = geometry.create_point_geom(lon, lat,
                                                lonlatprojection.osr_spref)

        # search for co-locating subgrid
        subgrid = self.locate_geometry_in_subgrids(point_geom)[0]

        x, y, = geometry.uv2xy(lon, lat,
                               lonlatprojection.osr_spref,
                               self.subgrids[subgrid].core.projection.osr_spref)

        return np.full_like(x, subgrid, dtype=(np.str, len(subgrid))), x, y


    def _lonlat2xy_subgrid(self, lon, lat, subgrid):
        """
        computes the projected coordinates in given subgrid.

        Parameters
        ----------
        lon : number
            longitude coordinate
        lat : number
            latitude coordinate
        subgrid : str
            acronym / subgrid ID to search within (speeding up)

        Returns
        -------
        subgrid : str
            acronym / subgrid ID
        x, y : int
            TPS grid coordinates
        """

        # set up spatial references
        p_grid = pyproj.Proj(self.subgrids[subgrid].core.projection.proj4)
        p_geo = pyproj.Proj(init="EPSG:4326")

        x, y, = pyproj.transform(p_geo, p_grid, lon, lat)

        return subgrid, x, y


    @abc.abstractmethod
    def create_tile(self, name):
        pass


    @abc.abstractmethod
    def get_tiletype(self, sampling):
        pass


    @abc.abstractmethod
    def get_tilesize(self, sampling):
        pass


    def search_tiles_in_roi(self,
                            geom_area=None,
                            extent=None,
                            osr_spref=None,
                            subgrid_ids=None,
                            coverland=False):
        """
        Search the tiles of the grid which intersect by the given area.

        Parameters
        ----------
        geom_area : geometry
            a polygon or multipolygon geometery object representing the ROI
        extent : list
            It is a list of coordinates representing either
                a) the rectangle-region-of-interest in the format of
                    [xmin, ymin, xmax, ymax]
                b) the tuple-list of points-of-intererst in the format of
                    [(x1, y1), (x2, y2), ...]
        osr_spref : OGRSpatialReference
            spatial reference of input coordinates in extent
        sgrid_ids : string or list of strings
            subgrid IDs, e.g. specifying over which continent
            you want to search.
            Default value is None for searching all subgrids.
        coverland : Boolean
            option to search for tiles covering land at any point in the tile

        Returns
        -------
        list
            return a list of  the overlapped tiles' name.
            If not found, return empty list.
        """

        # check input grids
        if subgrid_ids is None:
            subgrid_ids = self.subgrids.keys()
        if isinstance(subgrid_ids, str):
            subgrid_ids = [subgrid_ids]
        subgrid_ids = [x.upper() for x in subgrid_ids]
        if set(subgrid_ids).issubset(set(self.subgrids.keys())):
            subgrid_ids = list(subgrid_ids)
        else:
            raise ValueError("Invalid argument: grid must one of [ %s ]." %
                             " ".join(self.subgrids.keys()))

        if not geom_area and not extent:
            print("Error: either geom or extent should be given as the ROI.")
            return list()

        # obtain the geometry of ROI
        if not geom_area:
            if osr_spref is None:
                projection = TPSProjection(epsg=4326)
                osr_spref = projection.osr_spref
            geom_area = geometry.extent2polygon(extent, osr_spref)

        # load lat-lon spatial reference as the default
        geo_sr = TPSProjection(epsg=4326).osr_spref

        geom_sr = geom_area.GetSpatialReference()
        if geom_sr is None:
            geom_area.AssignSpatialReference(geo_sr)
        elif not geom_sr.IsSame(geo_sr):
            projected = geom_area.GetSpatialReference().IsProjected()
            if projected == 0:
                max_segment = 0.5
            elif projected == 1:
                max_segment = 50000
            else:
                raise Warning('Please check unit of geometry '
                              'before reprojection!')
            geom_area = geometry.transform_geometry(geom_area, geo_sr,
                                                    segment=max_segment)

        # intersect the given grid ids and the overlapped ids
        overlapped_grids = self.locate_geometry_in_subgrids(geom_area)
        subgrid_ids = list(set(subgrid_ids) & set(overlapped_grids))

        # finding tiles
        overlapped_tiles = list()
        for sgrid_id in subgrid_ids:
            overlapped_tiles.extend(
                self.subgrids[sgrid_id].search_tiles_in_geometry(geom_area,
                                                        coverland=coverland))
        return overlapped_tiles


class TiledProjection(object):

    """
    Class holding the projection and tiling definition of a
    tiled projection space.

    Parameters
    ----------
    Projection : Projection()
        A Projection object defining the spatial reference.
    tile_definition: TilingSystem()
        A TilingSystem object defining the tiling system.
        If None, the whole space is one single tile.
    """

    __metaclass__ = abc.ABCMeta

    staticdata = None

    def __init__(self, core, polygon_geog, tilingsystem=None):
        """
        Initialises a TiledProjection().

        Parameters
        ----------
        core : TPSCoreProperty
            defines core parameters of the (sub-) grid
        polygon_geog : OGRGeometry
            geometry defining the extent/outline of the subgrid.
            if not given, a single global subgrid is assigned to the grid.
        tilingsystem : TilingSystem
            optional; an instance of TilingSystem()
            if not given, a single global tile is assigned to the grid.
        """

        self.core = core
        self.polygon_geog = geometry.segmentize_geometry(polygon_geog,
                                                         segment=0.5)
        self.polygon_proj = geometry.transform_geometry(
            self.polygon_geog, self.core.projection.osr_spref)
        self.bbox_geog = geometry.get_geom_boundaries(
            self.polygon_geog, rounding=self.core.sampling / 1000000.0)
        self.bbox_proj = geometry.get_geom_boundaries(
            self.polygon_proj, rounding=self.core.sampling)

        if tilingsystem is None:
            tilingsystem = GlobalTile(self.core, 'TG', self.get_bbox_proj())
        self.tilesys = tilingsystem

    def __getattr__(self, item):
        '''
        short link for items of core
        '''
        if item in self.core.__dict__:
            return self.core.__dict__[item]
        else:
            return self.__dict__[item]


    def get_bbox_geog(self):
        """
        Returns the limits of the subgrid in the lon-lat-space.

        Returns
        -------
        tuple
            boundind box of subgrid
            as (lonmin, lonmax, latmin, latmax)
        """
        bbox = self.polygon_geog.GetEnvelope()
        return bbox


    def get_bbox_proj(self):
        """
        Returns the limits of the subgrid in the  rojected space.

        Returns
        -------
        tuple
            boundind box of subgrid
            as (xmin, xmax, ymin, ymax)
        """
        bbox = self.polygon_proj.GetEnvelope()
        return bbox


    def xy2lonlat(self, x, y):
        """
        Converts projected coordinates to longitude and latitude coordinates
        Parameters

        ----------
        x : number or list of numbers
            projected x coordinate(s) in metres
        y : number or list of numbers
            projected y coordinate(s) in metres

        Returns
        -------
        lon : float or list of floats
            longitude coordinate(s)
        lat : float or list of floats
            latitude coordinate(s)

        """
        # set up spatial references
        p_grid = pyproj.Proj(self.core.projection.proj4)
        p_geo = pyproj.Proj(init="EPSG:4326")

        lon, lat = pyproj.transform(p_grid, p_geo, x, y)

        return lon, lat


    def search_tiles_in_geometry(self, geom, coverland=True):
        """
        Search the tiles which are overlapping with the subgrid

        Parameters
        ----------
        geom : OGRGeometry
            A polygon geometry representing the region of interest.
        coverland : Boolean
            option to search for tiles covering land at any point in the tile

        Returns
        -------
        overlapped_tiles : list
            Return a list of the overlapped tiles' name.
            If not found, return empty list.
        """

        overlapped_tiles = list()

        # check if geom intersects subgrid
        if geom.Intersects(self.polygon_geog):
            # get intersect area with subgrid in latlon
            intersect = geom.Intersection(self.polygon_geog)
        else:
            return overlapped_tiles

        # get spatial reference of subgrid in grid projection
        grid_sr = self.projection.osr_spref

        # transform intersection geometry back to the spatial reference system
        # of the subgrid.
        # segmentise for high precision during reprojection.
        projected = intersect.GetSpatialReference().IsProjected()
        if projected == 0:
            max_segment = 0.5
        elif projected == 1:
            max_segment = 50000
        else:
            raise Warning('Please check unit of geometry before reprojection!')
        intersect = geometry.transform_geometry(intersect, grid_sr,
                                                segment=max_segment)

        # get envelope of the Geometry and cal the bounding tile of the
        envelope = intersect.GetEnvelope()
        x_min = int(envelope[0]) // self.core.tile_xsize_m \
            * self.core.tile_xsize_m
        x_max = (int(envelope[1]) // self.core.tile_xsize_m + 1) \
            * self.core.tile_xsize_m
        y_min = int(envelope[2]) // self.core.tile_ysize_m * \
            self.core.tile_ysize_m
        y_max = (int(envelope[3]) // self.core.tile_ysize_m + 1) * \
            self.core.tile_ysize_m

        # make sure x_min and y_min greater or equal 0
        x_min = 0 if x_min < 0 else x_min
        y_min = 0 if y_min < 0 else y_min

        # get overlapped tiles
        xr = np.arange(
            x_min, x_max + self.core.tile_xsize_m, self.core.tile_xsize_m)
        yr = np.arange(
            y_min, y_max + self.core.tile_ysize_m, self.core.tile_ysize_m)

        for x, y in itertools.product(xr, yr):
            geom_tile = geometry.extent2polygon(
                (x, y, x + self.core.tile_xsize_m,
                 y + self.core.tile_xsize_m), grid_sr)
            if geom_tile.Intersects(intersect):
                ftile = self.tilesys.point2tilename(x, y)
                if not coverland or self.tilesys.check_tile_covers_land(ftile):
                    overlapped_tiles.append(ftile)

        return overlapped_tiles


class TilingSystem(object):

    """
    Class defining the tiling system and providing methods for queries and handling.

    Parameters (BBM: init(stuff))
    ----------
    projection : :py:class:`Projection`
        A Projection object defining the spatial reference.
    tile_definition: TilingSystem
        A TilingSystem object defining the tiling system.
        If None, the whole space is one single tile.

    Attributes (BBM: stuff that needs to be explained)
    ----------
    extent_geog:
    """

    __metaclass__ = abc.ABCMeta

    def __init__(self, core, polygon_geog, x0, y0):
        """
        Initialises an TilingSystem class for a specified subgrid.

        Parameters
        ----------
        core : TPSCoreProperty
            defines core parameters of the (sub-) grid
        polygon_geog : OGRGeometry
            geometry defining the extent/outline of the subgrid
        x0 : int
            lower-left x (right) coordinates of the subgrid
        y0 : int
            lower-left y (up) coordinates of the subgrid
        """

        self.core = core
        self.x0 = x0
        self.y0 = y0
        self.xstep = self.core.tile_xsize_m
        self.ystep = self.core.tile_ysize_m
        self.polygon_proj = geometry.transform_geometry(
            polygon_geog, self.core.projection.osr_spref)
        self.bbox_proj = geometry.get_geom_boundaries(
            self.polygon_proj, rounding=self.core.sampling)

    def __getattr__(self, item):
        '''
        short link for items of core
        '''
        if item in self.core.__dict__:
            return self.core.__dict__[item]
        else:
            return self.__dict__[item]


    @abc.abstractmethod
    def create_tile(self, name=None, x=None, y=None):
        """
        Returns a Tile object of the grid.

        Parameters
        ----------
        name : str
            name of the tile
        x : int
            x (right) coordinate of a pixel located in the desired tile
            must to given together with y
        y : int
            y (up) coordinate of a pixel located in the desired tile
            must to given together with x

        Returns
        -------
        Tile
            object containing info of the specified tile.

        Notes
        -----
        either name, or x and y, must be given.
        """
        return


    def round_xy2lowerleft(self, x, y):
        """
        Returns the lower-left coordinates of the tile in which the point,
        defined by x and y coordinates (in metres), is located.

        Parameters
        ----------
        x : int
            x (right) coordinate in the desired tile
            must to given together with y
        y : int
            y (up) coordinate in the desired tile
            must to given together with x

        Returns
        -------
        llx, lly: int
            lower-left coordinates of the tile
        """

        llx = x // self.core.tile_xsize_m * self.core.tile_xsize_m
        lly = y // self.core.tile_ysize_m * self.core.tile_ysize_m
        return llx, lly


    @abc.abstractmethod
    def point2tilename(self, x, y):
        """
        Returns the name string of an Equi7Tile in which the point,
        defined by x and y coordinates (in metres), is located.

        Parameters
        ----------
        x : int
            x (right) coordinate in the desired tile
            must to given together with y
        y : int
            y (up) coordinate in the desired tile
            must to given together with x

        Returns
        -------
        str
            the tilename

        """
        return


    @abc.abstractmethod
    def _encode_tilename(self, llx, lly):
        """
        Encodes a tilename defined by the lower-left coordinates of the tile,
        using inherent information

        Parameters
        ----------
        llx : int
            Lower-left x coordinate.
        lly : int
            Lower-left y coordinate.

        Returns
        -------
        str
            the tilename
        """
        return


    @abc.abstractmethod
    def decode_tilename(self, tilename):
        """
        Returns the information assigned to the tilename

        Parameters
        ----------
        tilename : str
            the tilename

        Returns
        -------
        various
            features of the tiles
        """
        a = None
        return a


    def identify_tiles_overlapping_xybbox(self, bbox):
        """Light-weight routine that returns
           the name of tiles overlapping the bounding box.

        Parameters
        ----------
        bbox : list
            list of projected coordinates limiting the bounding box.
            scheme: [xmin, ymin, xmax, ymax]

        Return
        ------
        tilenames : list
            list of tilenames overlapping the bounding box
        """

        xmin, ymin, xmax, ymax = bbox
        if (xmin >= xmax) or (ymin >= ymax):
            raise ValueError("Check order of coordinates of bbox! "
                             "Scheme: [xmin, ymin, xmax, ymax]")

        tsize_x = self.core.tile_xsize_m
        factor_x = tsize_x
        tsize_y = self.core.tile_ysize_m
        factor_y = tsize_y

        llxs = list(
            range(xmin // tsize_x * factor_x, xmax // tsize_x * factor_x + 1,
                  factor_x))
        llys = list(
            range(ymin // tsize_y * factor_y, ymax // tsize_y * factor_y + 1,
                  factor_y))
        tx, ty = np.meshgrid(llxs, llys)
        tx = tx.flatten()
        ty = ty.flatten()

        tilenames = list()
        for i, _ in enumerate(tx):
            tilenames.append(
                self._encode_tilename(tx[i], ty[i]))

        return tilenames


    def create_tiles_overlapping_xybbox(self, bbox):
        """Light-weight routine that returns
           the name of tiles intersecting the bounding box.

        Parameters
        ----------
        bbox : list of numbers
            list of projected coordinates limiting the bounding box.
            scheme: [xmin, ymin, xmax, ymax]

        Return
        ------
        tiles : list of Tiles()
            list of Tiles() intersecting the bounding box,
            with .active_subset_px() holding indices of the tile that cover
            the bounding box.
        """

        tilenames = self.identify_tiles_overlapping_xybbox(bbox)
        tiles = list()

        for t in tilenames:

            tile = self.create_tile(name=t)
            le, te, re, be = tile.active_subset_px
            extent = tile.limits_m()

            # left_edge
            if extent[0] <= bbox[0]:
                le = (bbox[0] - extent[0]) / tile.core.sampling
            # top_edge
            if extent[1] <= bbox[1]:
                te = (bbox[1] - extent[1]) / tile.core.sampling
            # right_edge
            if extent[2] > bbox[2]:
                re = (
                    bbox[2] - extent[2] + self.core.tile_xsize_m) / tile.core.sampling
            # bottom_edge
            if extent[3] > bbox[3]:
                be = (
                    bbox[3] - extent[3] + self.core.tile_ysize_m) / tile.core.sampling

            # subset holding indices of the tile that cover the bounding box.
            tile.active_subset_px = le, te, re, be

            tiles.append(tile)

        return tiles


class Tile(object):
    """
    A tile in the TiledProjectedSystem, holding characteristics of the tile.
    """

    __metaclass__ = abc.ABCMeta

    def __init__(self, core, name, xll, yll):
        """
        Initialises a Tile().

        Parameters
        ----------
        core : TPSCoreProperty
            defines core parameters of the tile
        name : str
            name of the tile
        llx : int
            lower-left x (right) coordinate of the tile
        lly : int
            lower-left y (up) coordinate of the tile
        """

        self.core = core
        self.name = name
        self.typename = core.tiletype
        self.llx = xll
        self.lly = yll
        self.x_size_px = self.core.tile_xsize_m / self.core.sampling
        self.y_size_px = self.core.tile_ysize_m / self.core.sampling
        self._subset_px = (0, 0, self.x_size_px, self.y_size_px)

    def __getattr__(self, item):
        '''
        short link for items of core
        '''
        if item in self.core.__dict__:
            return self.core.__dict__[item]
        else:
            return self.__dict__[item]


    def shape_px(self):
        """
        Returns the shape of the pixel array

        Returns
        -------
        tuple
            shape of the tile's pixel array as (samples_x, samples_y)
        """

        return (self.x_size_px, self.y_size_px)


    def limits_m(self):
        """
        returns the limits of the tile in projected coordinates (in metres)

        Returns
        -------
        tuple
            limits in the terms of (xmin, ymin, xmax, ymax)
        """

        return (self.llx, self.lly,
                self.llx + self.core.tile_xsize_m, self.lly + self.core.tile_ysize_m)


    @property
    def active_subset_px(self):
        """
        holds indices of the active_subset_px-of-interest

        Returns
        -------
        tuple
            active subset as
            (xmin, ymin, xmax, ymax) =
            (left edge, top edge, right edge, bottom edge)
        """

        return self._subset_px


    @active_subset_px.setter
    def active_subset_px(self, limits):
        """
        changes the indices of the active_subset_px-of-interest,
        mostly to a smaller extent, for efficient reading

        Parameters
        ----------
        limits : tuple
            the limits of subsets as
            (xmin, ymin, xmax, ymax) =
            (left edge, top edge, right edge, bottom edge)

        """

        string = ['xmin', 'ymin', 'xmax', 'ymax']
        if len(limits) != 4:
            raise ValueError('Limits are not properly set!')

        _max = [self.x_size_px, self.y_size_px, self.x_size_px, self.y_size_px]

        for l, limit in enumerate(limits):
            if (limit < 0) or (limit > _max[l]):
                raise ValueError('{} is out of bounds!'.format(string[l]))

        xmin, ymin, xmax, ymax = limits

        if xmin >= xmax:
            raise ValueError('xmin >= xmax!')
        if ymin >= ymax:
            raise ValueError('ymin >= ymax!')

        self._subset_px = limits


    def geotransform(self):
        """
        returns the GDAL geotransform list

        Returns
        -------
        list
            a list contain the geotransform elements (no rotation specified)
            as (llx, x pixel spacing, 0, lly, 0, y pixel spacing)
        """

        geot = [self.llx, self.core.sampling, 0,
                self.lly + self.core.tile_ysize_m, 0, -self.core.sampling]

        return geot


    def ij2xy(self, i, j):
        """
        Returns the projected coordinates of a tile pixel in the TilingSystem
        for a given pixel pair defined by column and row

        Parameters
        ----------
        i : number
            pixel row number
        j : number
            pixel collumn number

        Returns
        -------
        x : number
            x coordinate in the projection
        y : number
            y coordinate in the projection
        """

        gt = self.geotransform()

        x = gt[0] + i * gt[1] + j * gt[2]
        y = gt[3] + i * gt[4] + j * gt[5]

        if self.core.sampling <= 1.0:
            precision = len(str(int(1.0 / self.core.sampling))) + 1
            return round(x, precision), round(y, precision)
        else:
            return x, y


    def xy2ij(self, x, y):
        """
        returns the column and row number (i, j)
        of a projection coordinate pair (x, y)

        Parameters
        ----------
        x : number
            x coordinate in the projection
        y : number
            y coordinate in the projection

        Returns
        -------
        i : integer
            pixel row number
        j : integer
            pixel column number
        """

        gt = self.geotransform()

        # TODO: check if 1) round-to-nearest-int or 2) round-down-to-int
        i = int(round(-1.0 * (gt[2] * gt[3] - gt[0] * gt[5] + gt[5] * x - gt[2] * y) /
                      (gt[2] * gt[4] - gt[1] * gt[5])))
        j = int(round(-1.0 * (-1 * gt[1] * gt[3] + gt[0] * gt[4] - gt[4] * x + gt[1] * y) /
                      (gt[2] * gt[4] - gt[1] * gt[5])))

        return i, j


    def get_geotags(self):
        """
        returns the geotags for given tile used as geo-information for GDAL

        Returns
        -------
        geotags : dict
            dict containing the geotransform and the spatial reference in WKT
            format
        """

        geotags = {'geotransform': self.geotransform(),
                   'spatialreference': self.core.projection.wkt}

        return geotags


class GlobalTile(Tile):

    __metaclass__ = abc.ABCMeta

    def __init__(self, core, name, bbox_polygon_proj):
        """
        Initialising a GlobalTile(), covering the whole extent of the subgrid

        Parameters
        ----------
        core : TPSCoreProperty
            defines core parameters of the tile/subgrid
        name : str
            defining the name of the GlobalTile()
        bbox_polygon_proj : tuple
            limits in projection spacve of the subgrid
            as (xmin, xmax, ymin, ymax)
        """

        super(GlobalTile, self).__init__(core, name, 0, 0)
        self.typename = 'TG'
        self.core.tiletype = self.typename
        self.core.tile_xsize_m = np.int((np.floor(bbox_polygon_proj[1] /
                                  self.core.sampling) * self.core.sampling) - \
                                 (np.ceil(bbox_polygon_proj[0] /
                                  self.core.sampling) * self.core.sampling))
        self.core.tile_ysize_m = np.int((np.floor(bbox_polygon_proj[3] /
                                  self.core.sampling) * self.core.sampling) - \
                                 (np.ceil(bbox_polygon_proj[2] /
                                  self.core.sampling) * self.core.sampling))
        self.x_size_px = self.core.tile_xsize_m / self.core.sampling
        self.y_size_px = self.core.tile_ysize_m / self.core.sampling
        self._subset_px = (0, 0, self.x_size_px, self.y_size_px)
