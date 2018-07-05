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
Code for osgeo geometry operations.
"""

import numpy as np

from osgeo import ogr
from osgeo import osr

def uv2xy(u, v, src_ref, dst_ref):
    """
    wrapper; reprojects a pair of point coordinates
    Parameters
    ----------
    u : number
        input coordinate ("Rechtswert")
    v : number
        input coordinate ("Hochwert")
    src_ref : SpatialReference
        osgeo spatial reference defining the input u, v coordinates
    dst_ref : SpatialReference
        osgeo spatial reference defining the output x, y coordinates

    Returns
    -------
    x : number
        output coordinate ("Rechtswert")
    y : : number
        output coordinate ("Hochwert")
    """
    tx = osr.CoordinateTransformation(src_ref, dst_ref)
    x, y, _ = tx.TransformPoint(u, v)
    return x, y


def create_multipoint_geom(u, v, osr_spref):
    """
    wrapper; creates multipoint geometry in given projection
    Parameters
    ----------
    u : list of numbers
        input coordinates ("Rechtswert")
    v : list of numbers
        input coordinates ("Hochwert")
    osr_spref : OGRSpatialReference
        spatial reference of given coordinates

    Returns
    -------
    OGRGeometry
        a geometry holding all points defined by (u, v)

    """
    point_geom = ogr.Geometry(ogr.wkbMultiPoint)
    point_geom.AssignSpatialReference(osr_spref)
    for p, _ in enumerate(u):
        point = ogr.Geometry(ogr.wkbPoint)
        point.AddPoint(u[p], v[p], 0)
        point_geom.AddGeometry(point)

    return point_geom


def create_point_geom(u, v, osr_spref):
    """
    wrapper; creates single point geometry in given projection

    Parameters
    ----------
    u : list of numbers
        input coordinates ("Rechtswert")
    v : list of numbers
        input coordinates ("Hochwert")
    osr_spref : OGRSpatialReference
        spatial reference of given coordinates

    Return
    ------
    OGRGeometry
        a geometry holding point defined by (u, v)

    """
    point_geom = ogr.Geometry(ogr.wkbPoint)
    point_geom.AddPoint(u, v)
    point_geom.AssignSpatialReference(osr_spref)

    return point_geom

def create_geometry_from_wkt(wkt_multipolygon, epsg=4326, segment=None):
    """
    return extent geometry from multipolygon defined by wkt string

    Parameters
    ----------
    wkt_multipolygon : string
        WKT text containing points of geometry (e.g. polygon)
    epsg : int
        EPSG code of spatial reference of the points.
    segment : float
        for precision: distance of longest segment of the geometry polygon
        in units of input osr_spref (defined by epsg)

    Return
    ------
    OGRGeometry
        a geometry holding the multipolygon and spatial reference

    """
    geom = ogr.CreateGeometryFromWkt(wkt_multipolygon)
    geo_sr = osr.SpatialReference()
    geo_sr.SetWellKnownGeogCS("EPSG:{}".format(str(epsg)))
    geom.AssignSpatialReference(geo_sr)

    # modify the geometry such it has no segment longer then the given distance
    if segment is not None:
        geom.Segmentize(segment)

    return geom


def open_geometry(fname, format="shapefile"):
    '''
    opens a geometry from a vector file.

    fname : string
        full path of the output file name
    format : string
        format name. currently only shape file is supported
    def __getitem__(key):
        return geoms[key]
    '''
    driver = ogr.GetDriverByName("ESRI Shapefile")
    ds = driver.Open(fname, 0)
    feature = ds.GetLayer(0).GetFeature(0)
    geom = feature.GetGeometryRef()

    out = geom.Clone()
    ds, feature, geom = None, None, None
    return out


def write_geometry(geom, fname, format="shapefile", segment=None):
    """
    writes a geometry to a vector file.

    parameters
    ----------
    geom : Geometry
        geometry object
    fname : string
        full path of the output file name
    format : string
        format name. currently only shape file is supported
    segment : float
        for precision: distance of longest segment of the geometry polygon
        in units of input osr_spref
    """

    drv = ogr.GetDriverByName("ESRI Shapefile")
    dst_ds = drv.CreateDataSource(fname)
    srs = geom.GetSpatialReference()

    dst_layer = dst_ds.CreateLayer("out", srs=srs)
    fd = ogr.FieldDefn('DN', ogr.OFTInteger)
    dst_layer.CreateField(fd)

    feature = ogr.Feature(dst_layer.GetLayerDefn())
    feature.SetField("DN", 1)

    # modify the geometry such it has no segment longer then the given distance
    if segment is not None:
        geom.Segmentize(segment)

    feature.SetGeometry(geom)
    dst_layer.CreateFeature(feature)

    dst_ds, feature, geom = None, None, None
    return


def transform_geometry(geometry, osr_spref, segment=0.5):
    """
    return extent geometry

    Parameters
    ----------
    geometry : OGRGeomtery
        geometry object
    osr_spref : OGRSpatialReference
        spatial reference to what the geometry should be transformed to
    segment : float
        for precision: distance in units of input osr_spref of longest
        segment of the geometry polygon

    Return
    ------
    OGRGeomtery
        a geometry represented in the target spatial reference

    """
    # to count the number of points of the polygon
    # print(geometry_out.GetGeometryRef(0).GetGeometryRef(0).GetPointCount())

    geometry_out = geometry.Clone()

    # modify the geometry such it has no segment longer then the given distance
    if segment is not None:
        geometry_out.Segmentize(segment)

    # transform geometry to new spatial reference system.
    geometry_out.TransformTo(osr_spref)

    geometry = None
    return geometry_out


def get_geom_boundaries(geometry, rounding=1.0):
    """
    returns the envelope of the geometry

    Parameters
    ----------
    geometry : Geometry
        geometry object
    rounding : float
        precision
    Returns
    -------
    list of numbers
        rounded coordinates of geometry-envelope

    """
    limits = geometry.GetEnvelope()
    limits = [int(x / rounding) * rounding for x in limits]
    return limits


def extent2polygon(extent, osr_spref, segment=None):
    """create a polygon geometry from extent.

    extent : list
        list of coordinates representing either
            a) the rectangle-region-of-interest in the format of
                [xmin, ymin, xmax, ymax]
            b) the list of points-of-interest in the format of
                [(x1, y1), (x2, y2), ...]
    osr_spref : OGRSpatialReference
        spatial reference of the coordinates in extent
    segment : float
        for precision: distance of longest segment of the geometry polygon
        in units of input osr_spref
    """

    if isinstance(extent[0], tuple):
        geom_area = _points2geometry(extent, osr_spref)
    else:
        area = [(extent[0], extent[1]),
                ((extent[0] + extent[2]) / 2., extent[1]),
                (extent[2], extent[1]),
                (extent[2], (extent[1] + extent[3]) / 2.),
                (extent[2], extent[3]),
                ((extent[0] + extent[2]) / 2., extent[3]),
                (extent[0], extent[3]),
                (extent[0], (extent[1] + extent[3]) / 2.)]

        edge = ogr.Geometry(ogr.wkbLinearRing)
        [edge.AddPoint(np.double(x), np.double(y)) for x, y in area]
        edge.CloseRings()

        # modify the geometry such it has no segment longer then the given distance
        if segment is not None:
            edge.Segmentize(segment)

        geom_area = ogr.Geometry(ogr.wkbPolygon)
        geom_area.AddGeometry(edge)

        geom_area.AssignSpatialReference(osr_spref)

    return geom_area


def _points2geometry(coords, osr_spref):
    """create a point geometry from a list of coordinate-tuples

    coords : list of tuples
        point-coords in terms of [(x1, y1), (x2, y2), ...]
    osr_spref : OGRSpatialReference
        spatial reference of the coordinates in extent

    """
    u = []
    v = []
    for co in coords:
        if len(co) == 2:
            u.append(co[0])
            v.append(co[1])

    return create_multipoint_geom(u, v, osr_spref)

