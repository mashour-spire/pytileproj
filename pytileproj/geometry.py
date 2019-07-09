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

from shapely.geometry import Polygon
from shapely.geometry import LineString
from shapely.ops import linemerge
from shapely.ops import unary_union
from shapely.ops import polygonize


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

    Returns
    -------
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
    wkt_multipolygon : str
        WKT text containing points of geometry (e.g. polygon)
    epsg : int
        EPSG code of spatial reference of the points.
    segment : float
        for precision: distance of longest segment of the geometry polygon
        in units of input osr_spref (defined by epsg)

    Returns
    -------
    OGRGeometry
        a geometry holding the multipolygon and spatial reference

    """
    geom = ogr.CreateGeometryFromWkt(wkt_multipolygon)
    geo_sr = osr.SpatialReference()
    geo_sr.SetWellKnownGeogCS("EPSG:{}".format(str(epsg)))
    geom.AssignSpatialReference(geo_sr)

    # modify the geometry such it has no segment longer then the given distance
    if segment is not None:
        geom = segmentize_geometry(geom, segment=segment)

    return geom


def open_geometry(fname, feature=0, format='shapefile'):
    '''
    opens a geometry from a vector file.

    fname : str
        full path of the output file name
    feature : int
        optional; order number of feature that should be returned;
        default is 0 for the first feature/object/geometry.
    format : str
        optional; format name. currently only shapefile is supported
    def __getitem__(key):
        return geoms[key]

    Returns
    -------
    OGRGeometry
        a geometry from the input file as indexed by the feature number
    '''

    if format == 'shapefile':
        drivername = "ESRI Shapefile"

    driver = ogr.GetDriverByName(drivername)
    ds = driver.Open(fname, 0)
    feature = ds.GetLayer(0).GetFeature(feature)
    geom = feature.GetGeometryRef()

    out = geom.Clone()
    ds, feature, geom = None, None, None
    return out


def write_geometry(geom, fname, format="shapefile", segment=None):
    """
    writes a geometry to a vector file.

    parameters
    ----------
    geom : OGRGeometry
        geometry object
    fname : str
        full path of the output file name
    format : str
        optional; format name. currently only shapefile is supported
    segment : float
        for precision: distance of longest segment of the geometry polygon
        in units of input osr_spref
    """

    if format == 'shapefile':
        drivername = "ESRI Shapefile"

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
        geom = segmentize_geometry(geom, segment=segment)

    feature.SetGeometry(geom)
    dst_layer.CreateFeature(feature)

    dst_ds, feature, geom = None, None, None
    return


def transform_geometry(geometry, osr_spref, segment=None):
    """
    return extent geometry

    Parameters
    ----------
    geometry : OGRGeometry
        geometry object
    osr_spref : OGRSpatialReference
        spatial reference to what the geometry should be transformed to
    segment : float
        for precision: distance in units of input osr_spref of longest
        segment of the geometry polygon

    Returns
    -------
    OGRGeometry
        a geometry represented in the target spatial reference

    """
    # to count the number of points of the polygon
    # print(geometry_out.GetGeometryRef(0).GetGeometryRef(0).GetPointCount())

    geometry_out = geometry.Clone()

    # modify the geometry such it has no segment longer then the given distance
    if segment is not None:
        geometry_out = segmentize_geometry(geometry_out, segment=segment)

    # transform geometry to new spatial reference system.
    geometry_out.TransformTo(osr_spref)

    geometry = None
    return geometry_out


def segmentize_geometry(geometry, segment=0.5):
    """
    segmentizes the lines of a geometry

    Parameters
    ----------
    geometry : OGRGeometry
        geometry object
    segment : float
        for precision: distance in units of input osr_spref of longest
        segment of the geometry polygon

    Returns
    -------
    OGRGeometry
        a congruent geometry realised by more vertices along its shape
    """

    geometry_out = geometry.Clone()

    geometry_out.Segmentize(segment)

    geometry = None
    return geometry_out


def intersect_geometry(geometry1, geometry2):
    """
    returns the intersection of two point or polygon geometries

    Parameters
    ----------
    geometry1, geometry2 : OGRGeometry
        geometry objects

    Returns
    -------
    intersection : OGRGeometry
        a geometry representing the intersection area
    """

    geometry1c=geometry1.Clone()
    geometry2c=geometry2.Clone()
    geometry1 = None
    geometry2 = None

    intersection = geometry1c.Intersection(geometry2c)

    return intersection


def check_lonlat_intersection(geometry1, geometry2):
    """
    checks if two polygon geometries intersect in lonlat space.
    geometry1 is split at the antimeridian
    (i.e. the 180 degree dateline)

    Parameters
    ----------
    geometry1 : OGRGeometry
        polygon geometry object in lonlat space
        is split by the antimeridian
    geometry2 : OGRGeometry
        polygon geometry object
        should be the large one

    Returns
    -------
    boolean
        does geometry1 intersect with geometry2?
    """

    area = get_lonlat_intersection(geometry1, geometry2)

    if area.Area() == 0.0:
        return False
    if area.Area() != 0.0:
        return True


def get_lonlat_intersection(geometry1, geometry2):
    """
    gets the intersect in lonlat space.
    geometry1 is split at the antimeridian
    (i.e. the 180 degree dateline)

    Parameters
    ----------
    geometry1 : OGRGeometry
        polygon geometry object in lonlat space
        is split by the antimeridian
    geometry2 : OGRGeometry
        geometry object
        should be the large one

    Returns
    -------
    boolean
        does geometry1 intersect with geometry2?
    """

    geometry1c=geometry1.Clone()
    geometry2c=geometry2.Clone()
    geometry1 = None
    geometry2 = None

    polygons = split_polygon_by_antimeridian(geometry1c)

    return polygons.Intersection(geometry2c)


def split_polygon_by_antimeridian(lonlat_polygon):
    """
    Function that splits a polygon at the antimeridian
    (i.e. the 180 degree dateline)

    Parameters
    ----------
    lonlat_polygon : OGRGeometry
        geometry object in lonlat space
        to be split by the antimeridian
    Returns
    -------
    wrapped_polygons : OGRGeometry
        MULTIPOLYGON comprising east and west parts of lonlat_polygon
        contains only one POLYGON if no intersect with antimeridian is given

    """

    # define the antimeridain
    antimeridian = LineString([(180, -90), (180, 90)])

    # use shapely for the splitting
    merged = linemerge([Polygon(lonlat_polygon.GetBoundary().GetPoints()).boundary, antimeridian])
    borders = unary_union(merged)
    polygons = polygonize(borders)

    # setup OGR multipolygon
    wrapped_polygons = ogr.Geometry(ogr.wkbMultiPolygon)
    geo_sr = osr.SpatialReference()
    geo_sr.SetWellKnownGeogCS("EPSG:4326")
    wrapped_polygons.AssignSpatialReference(geo_sr)

    # wrap the longitude coordinates
    # to get only longitudes out out [0, 180] or [-180, 0]
    for p in polygons:

        point_coords = p.exterior.coords[:]
        lons = [p[0] for p in point_coords]
        if np.mean(lons) > 180:
            wrapped_points = [(y[0] - 360, y[1], y[2]) if y[0] >= 180.0 else y for y in point_coords]
        else:
            wrapped_points = point_coords

        new_poly = Polygon(wrapped_points)
        wrapped_polygons.AddGeometry(ogr.CreateGeometryFromWkt(new_poly.wkt))

    return wrapped_polygons


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
    limits : list of numbers
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

    Returns
    -------
    geom_area : OGRGeometry
        a geometry representing the input extent as
        a) polygon-geometry when defined by a rectangle extent
        b) point-geometry when defined by extent through tuples of coordinates
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
            edge = segmentize_geometry(edge, segment=segment)

        geom_area = ogr.Geometry(ogr.wkbPolygon)
        geom_area.AddGeometry(edge)

        geom_area.AssignSpatialReference(osr_spref)

    return geom_area


def round_vertices_of_polygon(geometry, decimals=0):
    """
    'Cleans' the vertices of a polygon, so that it has rounded coordinates.

    Parameters
    ----------
    geometry : OGRGeometry
        a polygon geometry
    decimals : int
        optional. rounding precision. default is 0.

    Returns
    -------
    geometry_out : OGRGeometry
        a polygon geometry with rounded vertices
    """

    ring = geometry.GetGeometryRef(0)

    rounded_ring = ogr.Geometry(ogr.wkbLinearRing)

    n_points = ring.GetPointCount()

    for p in range(n_points):
        lon, lat, z = ring.GetPoint(p)
        rlon, rlat, rz = [np.round(lon, decimals=decimals),
                          np.round(lat, decimals=decimals),
                          np.round(z, decimals=decimals)]
        rounded_ring.AddPoint(rlon, rlat, rz)

    geometry_out = ogr.Geometry(ogr.wkbPolygon)
    geometry_out.AddGeometry(rounded_ring)

    return geometry_out


def _points2geometry(coords, osr_spref):
    """create a point geometry from a list of coordinate-tuples

    coords : list of tuples
        point-coords in terms of [(x1, y1), (x2, y2), ...]
    osr_spref : OGRSpatialReference
        spatial reference of the coordinates in extent

    Returns
    -------
    geom_area : OGRGeometry
        a point-geometry representing the input tuples of coordinates
    """

    u = []
    v = []
    for co in coords:
        if len(co) == 2:
            u.append(co[0])
            v.append(co[1])

    return create_multipoint_geom(u, v, osr_spref)


def setup_test_geom_spitzbergen():
    """
    Routine providing a geometry for testing

    Returns
    -------
    poly_spitzbergen : OGRGeometry
        4-corner polygon over high latitudes (is much curved on the globe)
    """

    ring_global = ogr.Geometry(ogr.wkbLinearRing)
    ring_global.AddPoint(8.391827331539572,77.35762113396143)
    ring_global.AddPoint(16.87007957357446,81.59290885863483)
    ring_global.AddPoint(40.5011949830408,79.73786853853339)
    ring_global.AddPoint(25.43098663332705,75.61353436967198)
    ring_global.AddPoint(8.391827331539572,77.35762113396143)

    poly_spitzbergen = ogr.Geometry(ogr.wkbPolygon)
    poly_spitzbergen.AddGeometry(ring_global)

    geom_global_wkt = '''GEOGCS[\"WGS 84\",
                           DATUM[\"WGS_1984\",
                                 SPHEROID[\"WGS 84\", 6378137, 298.257223563,
                                          AUTHORITY[\"EPSG\", \"7030\"]],
                                 AUTHORITY[\"EPSG\", \"6326\"]],
                           PRIMEM[\"Greenwich\", 0],
                           UNIT[\"degree\", 0.0174532925199433],
                           AUTHORITY[\"EPSG\", \"4326\"]]'''
    geom_global_sr = osr.SpatialReference()
    geom_global_sr.ImportFromWkt(geom_global_wkt)
    poly_spitzbergen.AssignSpatialReference(geom_global_sr)

    return poly_spitzbergen


def setup_geom_kamchatka():
    """
    Routine providing a geometry for testing

    Returns
    -------
    poly_kamchatka : OGRGeometry
        4-corner polygon close to the dateline at Kamchatka peninsula.
    """

    ring_global = ogr.Geometry(ogr.wkbLinearRing)
    ring_global.AddPoint(165.6045170932673, 59.05482187690058)
    ring_global.AddPoint(167.0124744118732, 55.02758744559601)
    ring_global.AddPoint(175.9512099050924, 54.36588084375806)
    ring_global.AddPoint(179.4591330039386, 56.57634572271662)
    ring_global.AddPoint(165.6045170932673, 59.05482187690058)

    poly_kamchatka = ogr.Geometry(ogr.wkbPolygon)
    poly_kamchatka.AddGeometry(ring_global)

    geom_global_wkt = '''GEOGCS[\"WGS 84\",
                           DATUM[\"WGS_1984\",
                                 SPHEROID[\"WGS 84\", 6378137, 298.257223563,
                                          AUTHORITY[\"EPSG\", \"7030\"]],
                                 AUTHORITY[\"EPSG\", \"6326\"]],
                           PRIMEM[\"Greenwich\", 0],
                           UNIT[\"degree\", 0.0174532925199433],
                           AUTHORITY[\"EPSG\", \"4326\"]]'''
    geom_global_sr = osr.SpatialReference()
    geom_global_sr.ImportFromWkt(geom_global_wkt)
    poly_kamchatka.AssignSpatialReference(geom_global_sr)

    return poly_kamchatka


def setup_test_geom_siberia_antimeridian_180plus():
    """
    Routine providing a geometry for testing

    Returns
    -------
    poly_siberia_antim_180plus : OGRGeometry
        4-corner polygon in Siberia, crossing the antimeridian
    """

    ring_global = ogr.Geometry(ogr.wkbLinearRing)
    ring_global.AddPoint(177.6584965942706,67.04864900747906)
    ring_global.AddPoint(179.0142461506587,65.34233852520839)
    ring_global.AddPoint(184.1800038679373,65.74423313395079)
    ring_global.AddPoint(183.1741580487398,67.46683765736415)
    ring_global.AddPoint(177.6584965942706,67.04864900747906)

    poly_siberia_antim_180plus = ogr.Geometry(ogr.wkbPolygon)
    poly_siberia_antim_180plus.AddGeometry(ring_global)

    geom_global_wkt = '''GEOGCS[\"WGS 84\",
                           DATUM[\"WGS_1984\",
                                 SPHEROID[\"WGS 84\", 6378137, 298.257223563,
                                          AUTHORITY[\"EPSG\", \"7030\"]],
                                 AUTHORITY[\"EPSG\", \"6326\"]],
                           PRIMEM[\"Greenwich\", 0],
                           UNIT[\"degree\", 0.0174532925199433],
                           AUTHORITY[\"EPSG\", \"4326\"]]'''
    geom_global_sr = osr.SpatialReference()
    geom_global_sr.ImportFromWkt(geom_global_wkt)
    poly_siberia_antim_180plus.AssignSpatialReference(geom_global_sr)

    return poly_siberia_antim_180plus


def setup_test_geom_siberia_alaska():
    """
    Routine providing a geometry for testing

    Returns
    -------
        poly_siberia_alaska : OGRGeometry
        4-corner polygon over Siberia and Alaska, crossing antimeridian
        and covering two Equi7Grid subgrids.
    """

    ring_global = ogr.Geometry(ogr.wkbLinearRing)
    ring_global.AddPoint(177.6545884597184,67.05574774066811)
    ring_global.AddPoint(179.0195867605756,65.33232820668778)
    ring_global.AddPoint(198.4723636216472,66.06909015550372)
    ring_global.AddPoint(198.7828129097253,68.14247939909886)
    ring_global.AddPoint(177.6545884597184,67.05574774066811)

    poly_siberia_alaska = ogr.Geometry(ogr.wkbPolygon)
    poly_siberia_alaska.AddGeometry(ring_global)

    geom_global_wkt = '''GEOGCS[\"WGS 84\",
                           DATUM[\"WGS_1984\",
                                 SPHEROID[\"WGS 84\", 6378137, 298.257223563,
                                          AUTHORITY[\"EPSG\", \"7030\"]],
                                 AUTHORITY[\"EPSG\", \"6326\"]],
                           PRIMEM[\"Greenwich\", 0],
                           UNIT[\"degree\", 0.0174532925199433],
                           AUTHORITY[\"EPSG\", \"4326\"]]'''
    geom_global_sr = osr.SpatialReference()
    geom_global_sr.ImportFromWkt(geom_global_wkt)
    poly_siberia_alaska.AssignSpatialReference(geom_global_sr)

    return poly_siberia_alaska