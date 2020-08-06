# Copyright (c) 2016,Vienna University of Technology,
# Department of Geodesy and Geoinformation
# All rights reserved.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL VIENNA UNIVERSITY OF TECHNOLOGY, DEPARTMENT OF
# GEODESY AND GEOINFORMATION BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
# BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
# IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.


"""
Tests for the UTMGrid().
"""
import unittest
import numpy as np
import numpy.testing as nptest

from pytileproj.utmgrid import UTMGrid
from pytileproj.geometry import setup_test_geom_spitzbergen
from pytileproj.geometry import setup_geom_kamchatka


# ### for testing at BBM machine
# import os
# # gdal 2
# os.environ["GDAL_DATA"] = r"C:\ProgramData\OSGeoW\share\gdal"
# os.environ["GDAL_DRIVER_PAT"] = r"C:\ProgramData\OSGeoW\bin\gdalplugins"
#
# # gdal 3
# os.environ["GDAL_DATA"] = r"C:\Program Files\GDAL\gdal-data"
# os.environ["GDAL_DRIVER_PAT"] = r"C:\Program Files\GDAL\gdalplugins"

class TestBaseViaUTMGrid(unittest.TestCase):

    def test_lonlat2xy_doubles(self):
        """
        Tests lonlat to xy projection using double numbers.
        """
        utm = UTMGrid(500)
        x_should = 433124.249310
        y_should = 5338921.352324
        lon, lat = 14.1, 48.2
        sgrid_id, x, y = utm.lonlat2xy(lon, lat)
        assert sgrid_id == 'Z33N'
        nptest.assert_allclose(x_should, x)
        nptest.assert_allclose(y_should, y)


    def test_lonlat2xy_numpy_array(self):
        """
        Tests lonlat to xy projection using numpy arrays.
        """
        utm = UTMGrid(500)
        x_should = np.array([507840.292027,
                             210029.47])
        y_should = np.array([4983717.660043,
                             6820022.61])
        lon = np.array([15.1, 3.564943])
        lat = np.array([-45.3, 61.405307])
        sgrid_id, x, y = utm.lonlat2xy(lon, lat)
        nptest.assert_array_equal(sgrid_id, np.array(['Z33S', 'Z32N']))
        nptest.assert_allclose(x_should, x)
        nptest.assert_allclose(y_should, y)


    def ctest_lonlat2xy_MGRS_numpy_array(self):
        """
        Tests lonlat to xy projection using numpy arrays.
        """
        utm = UTMGrid(500)
        x_should = np.array([507840.292027,
                             210029.47,
                             2394255.583836])
        y_should = np.array([4983717.660043,
                             6820022.61,
                             2000000.0])
        lon = np.array([15.1, 3.564943, 90.0])
        lat = np.array([-45.3, 61.405307, -86.45])
        sgrid_id, x, y = utm.lonlat2xy_MGRS(lon, lat)
        nptest.assert_array_equal(sgrid_id, np.array(['Z33G', 'Z32V', 'Z00B']))
        nptest.assert_allclose(x_should, x)
        nptest.assert_allclose(y_should, y)


    def test_lonlat2xy_numpy_array_subgrid(self):
        """
        Tests lonlat to xy projection using numpy arrays.
        """
        utm = UTMGrid(500)
        x_should = np.array([492159.707973])
        y_should = np.array([5016282.339957])
        lon = np.array([-15.1])
        lat = np.array([45.3])
        sgrid_id, x, y = utm.lonlat2xy(lon, lat, subgrid='Z28N')
        nptest.assert_array_equal(sgrid_id, np.array(['Z28N']))
        nptest.assert_allclose(x_should, x)
        nptest.assert_allclose(y_should, y)


    def test_lonlat2xy_numpy_array_wrong_subgrid(self):
        """
        Tests lonlat to xy projection giving a wrong subgrid.
        """
        utm = UTMGrid(500)
        x_should = np.array([-5028208.49022517])
        y_should = np.array([20666035.74265257])
        lon = np.array([-15.1])
        lat = np.array([45.3])
        sgrid_id, x, y = utm.lonlat2xy(lon, lat, subgrid='Z44S')
        nptest.assert_array_equal(sgrid_id, np.array(['Z44S']))
        nptest.assert_allclose(x_should, x)
        nptest.assert_allclose(y_should, y)


    def test_xy2lonlat_doubles(self):
        """
        Tests xy to lonlat projection using double numbers.
        """
        utm = UTMGrid(500)
        x = 458119.890658
        y = 6312037.887621
        lon_should, lat_should = -105.45, -33.33
        lon, lat = utm.Z13S.xy2lonlat(x, y)
        nptest.assert_allclose(lon_should, lon)
        nptest.assert_allclose(lat_should, lat)


    def test_xy2lonlat_numpy_array(self):
        """
        Tests xy to lonlat projection using numpy arrays.
        """
        utm = UTMGrid(500)
        x = np.array([458119.890658])
        y = np.array([6312037.887621])
        lon_should, lat_should = -105.68849338, 56.95006105
        lon, lat = utm.Z13N.xy2lonlat(x, y)
        nptest.assert_allclose(lon_should, lon)
        nptest.assert_allclose(lat_should, lat)


    def test_ij2xy(self):
        """
        Tests tile indices to xy coordination in the subgrid projection.
        """
        utm = UTMGrid(500)
        x_should = 481500
        y_should = 9270500
        tile = utm.Z18N.tilesys.create_tile(x=481746, y=9270569)
        x, y = tile.ij2xy(963, 659)
        nptest.assert_allclose(x_should, x)
        nptest.assert_allclose(y_should, y)

        # pixel center position
        x_should = 750
        y_should = 9000750
        x, y = tile.ij2xy(1, 1, lowerleft=True, offset='center')
        nptest.assert_allclose(x_should, x)
        nptest.assert_allclose(y_should, y)

    def test_xy2ij(self):
        """
        Tests xy to tile array indices.
        """
        utm = UTMGrid(500)
        column_should = 963
        row_should = 659
        tile = utm.Z18N.tilesys.create_tile(x=481746, y=9270569)
        column, row = tile.xy2ij(481500, 9270500)
        nptest.assert_allclose(column_should, column)
        nptest.assert_allclose(row_should, row)


    def test_lonlat2ij_in_tile(self):
        """
        Tests xy to tile array indices.
        """
        utm = UTMGrid(500)
        column_should = 604
        row_should = 276
        tile_should = 'Z34N500M_E000N048T6'
        tilename, i, j = utm.lonlat2ij_in_tile(18.507, 44.571, lowerleft=True)
        nptest.assert_allclose(i, column_should)
        nptest.assert_allclose(j, row_should)
        nptest.assert_equal(tilename, tile_should)


    def test_xy2ij_in_tile(self):
        """
        Tests xy to tile array indices.
        """
        utm = UTMGrid(500)
        column_should = 604
        row_should = 276
        tile_should = 'Z34N500M_E000N048T6'
        tilename, i, j = utm.Z34N.tilesys.xy2ij_in_tile(x=302051, y=4938318, lowerleft=True)
        nptest.assert_allclose(i, column_should)
        nptest.assert_allclose(j, row_should)
        nptest.assert_equal(tilename, tile_should)


    def test_decode_tilename(self):
        """
        Tests the decoding of tilenames.
        """
        utm_500 = UTMGrid(500)
        utm_10 = UTMGrid(10)

        assert utm_500.Z24N.tilesys.decode_tilename('Z24N500M_E000N006T6') == \
               ('Z24N', 500, 600000, 0, 600000, 'T6')
        assert utm_10.Z01N.tilesys.decode_tilename('Z01N010M_E085N091T1') == \
               ('Z01N', 10, 100000, 8500000, 9100000, 'T1')

        assert utm_500.Z24N.tilesys.decode_tilename('Z24N500M_E000N006T6') == \
               ('Z24N', 500, 600000, 0, 600000, 'T6')
        assert utm_10.Z01N.tilesys.decode_tilename('Z01N010M_E085N091T1') == \
               ('Z01N', 10, 100000, 8500000, 9100000, 'T1')

        with nptest.assert_raises(ValueError) as excinfo:
            utm_10.Z01N.tilesys.decode_tilename('E000N006T6')
        assert str(excinfo.exception).startswith(
            '"tilename" is not properly defined!')


    def test_find_overlapping_tilenames(self):
        """
        Tests search for tiles which share the same extent_m but
        with different resolution and tilecode.
        """
        utm_500 = UTMGrid(500)
        utm_10 = UTMGrid(10)

        tiles1_should = ['Z33N025M_E000N006T3', 'Z33N025M_E000N009T3',
                         'Z33N025M_E003N006T3', 'Z33N025M_E003N009T3']
        tiles1 = utm_500.Z33N.tilesys.get_congruent_tiles_from_tilename(
            'Z33N500M_E000N006T6',
            target_sampling=25)
        assert sorted(tiles1) == sorted(tiles1_should)

        tiles2_should = ['E000N006T3', 'E000N009T3', 'E003N006T3',
                         'E003N009T3']
        tiles2 = utm_500.Z33N.tilesys.get_congruent_tiles_from_tilename('E000N006T6',
                                                                        target_tiletype='T3')
        assert sorted(tiles2) == sorted(tiles2_should)

        tiles3_should = ['Z33N500M_E000N012T6']
        tiles3 = utm_10.Z33N.tilesys.get_congruent_tiles_from_tilename('E004N015T1',
                                                                       target_sampling=500)
        assert sorted(tiles3) == sorted(tiles3_should)

        tiles4_should = ['E003N009T3']
        tiles4 = utm_10.Z33N.tilesys.get_congruent_tiles_from_tilename('E004N011T1',
                                                                       target_tiletype='T3')
        assert sorted(tiles4) == sorted(tiles4_should)


    def test_search_tiles_lon_lat_extent(self):
        """
        Tests searching for tiles with input of lon lat extent
        """
        # TODO: STILL NEEDS TO BE CAREFULLY CHECKED! OVER THE ARCTICS! Also need to be
        # adapted for coverland!
        # two tiles are different between GDAL3 and GDAL2!

        utm = UTMGrid(500)

        tiles = utm.search_tiles_in_roi(bbox=[(-10, 80), (5, 85)],
                                        coverland=True)

        desired_tiles = ['Z31N500M_E000N084T6', 'Z31N500M_E000N090T6',
                         'Z00Y500M_E018N012T6', 'Z00Z500M_E018N012T6',
                         'Z29N500M_E000N084T6', 'Z29N500M_E000N090T6',
                         'Z30N500M_E000N084T6', 'Z30N500M_E000N090T6']
        assert sorted(tiles) == sorted(desired_tiles)

        tiles_all = utm.search_tiles_in_roi(bbox=[(-179.9, -89.9), (179.9, 89.9)],
                                            coverland=True)

        assert len(tiles_all) == 3638



    def test_search_tiles_lon_lat_extent_by_points(self):
        """
        Tests searching for tiles with input of lon lat points
        """
        # TODO: STILL NEEDS TO BE CAREFULLY CHECKED! Also need to be
        # adapted for coverland!

        utm = UTMGrid(500)

        tiles = utm.search_tiles_in_roi(points=[(10, 40), (5, 50), (-90.9, -1.2), (-175.2, 66)],
                                        coverland=True)

        desired_tiles = ['Z31N500M_E006N054T6', 'Z01N500M_E000N072T6',
                         'Z32N500M_E000N042T6', 'Z15S500M_E006N096T6']
        assert sorted(tiles) == sorted(desired_tiles)


    def test_search_tiles_spitzbergen(self):
        """
        Tests the tile searching over Spitzbergen in the polar zone; ROI #defined
        by a 4-corner polygon over high latitudes (is much curved on the globe).
        """
        # TODO: STILL NEEDS TO BE CAREFULLY CHECKED! Also need to be
        # adapted for coverland!

        grid = UTMGrid(500)

        spitzbergen_geom = setup_test_geom_spitzbergen()
        spitzbergen_geom_tiles = sorted(
            ['Z31N500M_E006N084T6', 'Z33N500M_E000N084T6', 'Z33N500M_E000N090T6',
             'Z33N500M_E006N084T6', 'Z33N500M_E006N090T6', 'Z35N500M_E000N078T6',
             'Z35N500M_E000N084T6', 'Z35N500M_E000N090T6', 'Z35N500M_E006N084T6',
             'Z37N500M_E000N084T6'])
        tiles = sorted(grid.search_tiles_in_roi(spitzbergen_geom,
                                                coverland=False))
        assert sorted(tiles) == sorted(spitzbergen_geom_tiles)

        spitzbergen_geom_tiles = sorted(
            ['Z31N500M_E006N084T6', 'Z33N500M_E000N084T6', 'Z33N500M_E000N090T6',
             'Z33N500M_E006N084T6', 'Z33N500M_E006N090T6', 'Z35N500M_E000N078T6',
             'Z35N500M_E000N084T6', 'Z35N500M_E000N090T6', 'Z35N500M_E006N084T6',
             'Z37N500M_E000N084T6'])
        tiles = sorted(grid.search_tiles_in_roi(spitzbergen_geom,
                                                coverland=True))
        assert sorted(tiles) == sorted(spitzbergen_geom_tiles)


    def test_search_tiles_kamchatka(self):
        """
        Tests the tile searching over Kamchatka in far east Sibiria;

        This test is especially nice, as it contains also a tile that covers both,
        the ROI and the continental zone, but the intersection of the tile and
        the ROI is outside of the zone.

        Furthermore, it also covers zones that consist of a multipolygon, as it
        is located at the 180deg/dateline.
        """

        grid = UTMGrid(500)

        kamchatka_geom = setup_geom_kamchatka()
        kamchatka_geom_tiles = sorted(
            ['Z58N500M_E000N060T6', 'Z58N500M_E006N060T6',
             'Z59N500M_E000N060T6', 'Z59N500M_E006N060T6',
             'Z60N500M_E000N060T6', 'Z60N500M_E006N060T6'])
        tiles = sorted(
            grid.search_tiles_in_roi(kamchatka_geom, coverland=False))

        assert sorted(tiles) == sorted(kamchatka_geom_tiles)


    def test_identify_tiles_overlapping_xybbox(self):
        """
        Tests identification of tiles covering a bounding box
        given in UTM coordinats
        """
        # TODO: implemented check for out of bounds/zone tiles
        utm_500 = UTMGrid(500)
        utm_10 = UTMGrid(10)

        tiles1_should = ['Z33N500M_E000N054T6', 'Z33N500M_E006N054T6']

        tiles2_should = ['Z33N010M_E005N058T1', 'Z33N010M_E005N059T1',
                         'Z33N010M_E005N060T1', 'Z33N010M_E005N061T1']

        tiles1 = utm_500.Z33N.tilesys.identify_tiles_overlapping_xybbox(
            [559745, 5852882, 611111, 5952882])

        tiles2 = utm_10.Z33N.tilesys.identify_tiles_overlapping_xybbox(
            [559745, 5852882, 571111, 6102882])

        assert sorted(tiles1) == sorted(tiles1_should)
        assert sorted(tiles2) == sorted(tiles2_should)


    def test_get_covering_tiles(self):
        """
        Tests the search for co-locating tiles of other type.

        """
        utm_500 = UTMGrid(500)
        utm_10 = UTMGrid(10)

        fine_tiles = ['Z33N010M_E005N058T1', 'Z33N010M_E005N059T1',
                      'Z33N010M_E005N060T1', 'Z33N010M_E005N061T1']

        target_tiletype = utm_500.get_tiletype()
        target_sampling = utm_500.core.sampling

        # invoke the results as tile name in short form
        coarse_tiles_shortform = utm_10.Z33N.tilesys.collect_congruent_tiles(fine_tiles,
                                                                             target_tiletype=target_tiletype)

        # invoke the results as tile name in long form
        coarse_tiles_longform = utm_10.Z33N.tilesys.collect_congruent_tiles(fine_tiles,
                                                                            target_sampling=target_sampling)

        assert sorted(coarse_tiles_shortform) == ['E000N054T6', 'E000N060T6']
        assert sorted(coarse_tiles_longform) == ['Z33N500M_E000N054T6', 'Z33N500M_E000N060T6']


if __name__ == '__main__':
    unittest.main()