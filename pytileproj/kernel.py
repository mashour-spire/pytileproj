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

import numpy as np
from scipy.stats import norm

from astropy.convolution import kernels


class Kernel(object):

    def __init__(self, size=3):
        self.size = size
        self.array = np.array((size, size)).astype(np.float32)

    def hamming_window(self):
        '''
        creates a 2D hamming window

        '''

        window = np.hamming(self.size)
        kernel = np.outer(window, window)
        kernel /= np.sum(kernel)
        self.array = kernel

    def mean_kernel(self):
        '''
        creates a 2D mean average window

        '''

        kernel = np.full(
            (self.size, self.size), 1. / self.size**2, dtype=np.float32)
        self.array = kernel

    def gauss_fwhm_kernel(self, fwhm=2.0):
        '''
        Creates a 2D gauss window defined by the desired full width at half maximum (FWHM).
        A window is created that features the centres' half weight at a distance
        of FWHM/2 from the center.

        '''

        # factor between sigma of filter and FWHM
        kappa = 2.0 * np.sqrt(2 * np.log(2))

        # "sigma" auf gaussian kernel
        sdev = fwhm / kappa

        # recommended size of array carrying the kernel
        rec_size = kernels._round_up_to_odd_integer(4.0 * sdev)
        if self.size != rec_size:
            print(
                "Kernel size is not properly chosen, to reflect gaussian function!")

        # normalized gaussian kernel
        GaussKernel = kernels.Gaussian2DKernel(
            stddev=sdev, x_size=self.size, y_size=self.size)
        kernel = GaussKernel._array / np.sum(GaussKernel._array)

        self.array = kernel

    def gauss_weightratio_kernel(self, dist=1, row=0.5):
        '''
        Creates a 2D gauss window defined by the desired sum of weighting
        (ratio of weigthing, roe) between (centre - dist) and (centre + dist).
        This means, the inter-quartile range (IQR) has the
        extent of 2*dist.
        '''

        # interquartile range (IQR) of normal distribution
        iqr = 1.0 / norm.ppf(0.5 + row / 2)

        # standard deviaton of the gaussian kernel
        sdev = dist * iqr / 2

        # recommended size of array carrying the kernel
        rec_size = kernels._round_up_to_odd_integer(4.0 * sdev)
        if self.size != rec_size:
            print(
                "Kernel size is not properly chosen, to reflect gaussian function!")

        # normalized gaussian kernel
        GaussKernel = kernels.Gaussian2DKernel(
            stddev=sdev, x_size=self.size, y_size=self.size)
        kernel = GaussKernel._array / np.sum(GaussKernel._array)

        self.array = kernel
