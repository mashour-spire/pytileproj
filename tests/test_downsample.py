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
Tests for the PixelDownsampler class.
"""


import numpy.testing as nptest
import numpy as np
from pytileproj.downsample import calc_pixel_index_pattern
from pytileproj.downsample import translate_pixelmaps
from pytileproj.downsample import fast_mask_counting
from pytileproj.downsample import PixelDownsampler
from pytileproj.downsample import downsampling_gauss_filter


def test_calc_pixel_index_pattern():

    a_out = calc_pixel_index_pattern(75, 500)
    a_should = [7, 6, 7]
    nptest.assert_array_equal(a_should, a_out)

    b_out = calc_pixel_index_pattern(13, 4)
    b_should = [0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0]
    nptest.assert_array_equal(b_should, b_out)


def test_translate_pixelmaps():

    spacing_fine = 10
    spacing_coarse = 75
    target_bbox = [4800000, 1200000, 4800150, 1200300]
    _, _, pix_count, cum_pixels_x, cum_pixels_y, needed_bbox, needed_shape = \
        translate_pixelmaps(
            spacing_fine, spacing_coarse, target_bbox, correct_boundary=True)

    pix_count_should = (
        [[49, 56, 49, 56, 49, 56],
         [56, 64, 56, 64, 56, 64],
         [49, 56, 49, 56, 49, 56],
         [56, 64, 56, 64, 56, 64],
         [49, 56, 49, 56, 49, 56],
         [56, 64, 56, 64, 56, 64],
         [49, 56, 49, 56, 49, 56],
         [56, 64, 56, 64, 56, 64]])
    nptest.assert_array_equal(pix_count_should, pix_count)

    cum_pixels_x_should = np.array([7, 15, 22, 30, 37, 45])
    nptest.assert_array_equal(cum_pixels_x_should, cum_pixels_x)

    cum_pixels_y_should = np.array([7, 15, 22, 30, 37, 45, 52, 60])
    nptest.assert_array_equal(cum_pixels_y_should, cum_pixels_y)

    needed_bbox_should = [4799850, 1199850, 4800300, 1200450]
    nptest.assert_array_equal(needed_bbox_should, needed_bbox)

    needed_shape_should = (60, 45)
    nptest.assert_array_equal(needed_shape_should, needed_shape)


def test_fast_mask_counting():

    mask = np.ones((9, 6), dtype='int8')
    mask[1:4, [1, 4]] = 0
    mask[3, 5] = 0
    coarse_shape = (6, 4)
    pattern_x = np.array((1, 3, 4, 6), dtype='uint16')
    pattern_y = np.array((1, 3, 4, 6, 7, 9), dtype='uint16')

    a_out = fast_mask_counting(mask, coarse_shape, pattern_x, pattern_y)

    a_should = (
        [[0, 0, 0, 0],
         [0, 2, 0, 2],
         [0, 1, 0, 2],
         [0, 0, 0, 0],
         [0, 0, 0, 0],
         [0, 0, 0, 0]])

    nptest.assert_array_equal(a_should, a_out)


def test_downsampling_gauss_filter():

    a = np.arange(6 * 6, dtype=np.float32).reshape((6, 6)) * 77 - 566
    a[1, 3] = np.nan

    a_out = downsampling_gauss_filter(
        a, no_data_value=None, use_hard_coded_kernel=False)

    a_should = np.array(
        [[-445.795013, -385.967133, -331.785095, -282.566467, -169.724945, -95.139290],
         [-86.827858, -27.000000, 39.146229,
             127.000000, 214.853775, 263.827850],
         [375.172150, 435.000000, 532.150391,
             654.122620, 694.210571, 725.827881],
         [837.172119, 897.000000, 974.000000,
             1051.000000, 1128.000000, 1187.827881],
         [1299.172119, 1359.000000, 1436.000000,
             1513.000000, 1590.000000, 1649.827881],
         [1658.139282, 1717.967163, 1794.967163, 1871.967163, 1948.967163, 2008.795044]])

    nptest.assert_allclose(a_should, a_out, atol=1e-7)


def test_downsample_via_pixel_indices():

    a = np.arange(9 * 6, dtype=np.float32).reshape((9, 6)) - 56
    a[1, 3] = np.nan

    ds = PixelDownsampler(200, 300, [0, 0, 1200, 1800])
    a_out = ds.downsample_via_pixel_indices(a)

    a_should = np.array([
        [-53.65834427, -52.34365845, -50.62223434, -49.67818069],
        [-46.66547775, -45.12936783, -43.07934952, -42.46389008],
        [-37.66547775, -36.35079193, -34.62936783, -33.68531418],
        [-28.66547775, -27.5, -26., -24.83452225],
        [-19.66547775, -18.5, -17., -15.83452225],
        [-12.67261124, -11.50713348, -10.00713348, -8.84165573]])

    nptest.assert_allclose(a_should, a_out)
