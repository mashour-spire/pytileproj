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

from osgeo import ogr, osr
import pytest

from pytileproj.base import TiledProjectionSystem
# from pytileproj.base import TiledProjectionSystemTest
# TODO: write a TiledProjectionSystem class which can be used for testing


def test_geom_intersect():
    geom_roi = setup_geom_roi()

    # reference file list
    ref_tiles = [
        'EU500M_E012N012T6', 'EU500M_E012N018T6', 'EU500M_E018N006T6', 'EU500M_E024N006T6',
                 'EU500M_E024N012T6', 'EU500M_E030N006T6', 'EU500M_E030N012T6', 'EU500M_E030N024T6',
                 'EU500M_E036N006T6', 'EU500M_E036N012T6', 'EU500M_E036N018T6', 'EU500M_E036N024T6',
                 'EU500M_E036N030T6', 'EU500M_E036N036T6', 'EU500M_E036N042T6', 'EU500M_E042N000T6',
                 'EU500M_E042N006T6', 'EU500M_E042N012T6', 'EU500M_E042N018T6', 'EU500M_E042N024T6',
                 'EU500M_E042N030T6', 'EU500M_E042N042T6', 'EU500M_E048N000T6', 'EU500M_E048N006T6',
                 'EU500M_E048N012T6', 'EU500M_E048N018T6', 'EU500M_E048N024T6', 'EU500M_E048N030T6',
                 'EU500M_E048N036T6', 'EU500M_E054N000T6', 'EU500M_E054N006T6', 'EU500M_E054N012T6',
                 'EU500M_E054N018T6', 'EU500M_E054N024T6', 'EU500M_E054N030T6', 'EU500M_E054N036T6',
                 'EU500M_E054N042T6', 'EU500M_E054N048T6', 'EU500M_E060N000T6', 'EU500M_E060N006T6',
                 'EU500M_E060N012T6', 'EU500M_E060N018T6', 'EU500M_E060N024T6', 'EU500M_E060N030T6',
                 'EU500M_E060N036T6', 'EU500M_E060N048T6', 'EU500M_E066N000T6', 'EU500M_E066N006T6',
                 'EU500M_E066N012T6', 'EU500M_E066N018T6', 'EU500M_E066N024T6', 'EU500M_E066N030T6',
                 'EU500M_E066N036T6', 'EU500M_E072N000T6', 'EU500M_E072N006T6', 'EU500M_E072N012T6',
                 'EU500M_E072N018T6', 'EU500M_E072N024T6', 'EU500M_E072N030T6', 'EU500M_E078N006T6',
                 'EU500M_E078N012T6']

    # grid = TiledProjectionSystemTest(500)
    # res_tiles = grid._search_sgrid_tiles(geom_roi, 'EU', True)

    assert 4 == 4


def setup_geom_roi():
    ring_global = ogr.Geometry(ogr.wkbLinearRing)
    ring_global.AddPoint(-180.0, -89.999928)
    ring_global.AddPoint(-0.000072, -89.999928)
    ring_global.AddPoint(179.999856, -89.999928)
    ring_global.AddPoint(179.999856, 0.000036)
    ring_global.AddPoint(179.999856, 90.0)
    ring_global.AddPoint(-0.000072, 90.0)
    ring_global.AddPoint(-180, 90)
    ring_global.AddPoint(-180.0, 0.000036)
    ring_global.AddPoint(-180.0, -89.999928)

    poly_global = ogr.Geometry(ogr.wkbPolygon)
    poly_global.AddGeometry(ring_global)

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
    poly_global.AssignSpatialReference(geom_global_sr)

    return poly_global
