# Code reproduced from buildH by Patrick Fuchs under BSD3 LICENSCE
# see below for original copyright notice
#
# BSD 3-Clause License
#
# Copyright (c) 2019, Patrick FUCHS
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
#* Redistributions of source code must retain the above copyright notice, this
#  list of conditions and the following disclaimer.
#
#* Redistributions in binary form must reproduce the above copyright notice,
#  this list of conditions and the following disclaimer in the documentation
#  and/or other materials provided with the distribution.
#
#* Neither the name of the copyright holder nor the names of its
#  contributors may be used to endorse or promote products derived from
#  this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

"""
Unit tests for buildH_calcOP.
Test functions from module hydrogen
"""

import collections
import numpy as np
from numpy.testing import assert_almost_equal
import pytest

from fast_forward import hydrogen


class TestComputeHydrogen:
    """Test class for the functions which reconstruct hydrogen from atom and helpers."""

    # Data tuple for get_CH() test
    Data = collections.namedtuple('Data', ['atom', 'helper1', 'helper2', 'helper3', 'H_coord'])

    @pytest.mark.parametrize('data',[
        (Data(np.array([26.87, 45.09, 26.03], dtype=np.float32),
              np.array([27.85, 46.07, 25.38], dtype=np.float32),
              np.array([27.65, 44.79, 27.32], dtype=np.float32),
              np.array([25.53, 45.59, 26.01], dtype=np.float32),
              np.array([26.6186898, 44.1429609, 25.552444], dtype=np.float32))),
        (Data(np.array([09.37, 46.61, 21.07], dtype=np.float32),
              np.array([08.62, 47.02, 19.80], dtype=np.float32),
              np.array([09.71, 47.93, 21.77], dtype=np.float32),
              np.array([08.64, 45.68, 21.88], dtype=np.float32),
              np.array([10.2773800, 46.049355, 20.8454200], dtype=np.float32))),
        (Data(np.array([30.83, 37.09, 27.29], dtype=np.float32),
              np.array([30.72, 36.37, 25.95], dtype=np.float32),
              np.array([31.63, 38.33, 26.88], dtype=np.float32),
              np.array([29.69, 37.27, 28.14], dtype=np.float32),
              np.array([31.3003330, 36.458855, 28.0440180], dtype=np.float32))),
        (Data(np.array([40.79, 32.48, 25.71], dtype=np.float32),
              np.array([41.58, 32.30, 24.41], dtype=np.float32),
              np.array([39.48, 31.71, 25.54], dtype=np.float32),
              np.array([41.50, 31.96, 26.85], dtype=np.float32),
              np.array([40.6256800, 33.5420680, 25.8919559], dtype=np.float32))),
    ])
    def test_get_CH(self, data):
        """Test for get_CH().
        3 helpers are needed and 1 hydrogen is returned.
        Parameters
        ----------
        data: namedtuple
                structure holding input and reference data.
        """
        assert_almost_equal(hydrogen.get_CH(data.atom, data.helper1, data.helper2, data.helper3),
                            data.H_coord)


    # Data tuple for get_CH2() test
    Data = collections.namedtuple('Data', ['atom', 'helper1', 'helper2', 'H1_coord', 'H2_coord'])

    @pytest.mark.parametrize('data',[
        (Data(np.array([32.42, 45.49, 26.87], dtype=np.float32),
              np.array([33.40, 46.51, 27.28], dtype=np.float32),
              np.array([31.74, 45.57, 25.50], dtype=np.float32),
              np.array([31.6241872, 45.5047426, 27.6146938], dtype=np.float32),
              np.array([32.9347447, 44.5299254, 26.9072779], dtype=np.float32))),
        (Data(np.array([28.81, 46.72, 30.22], dtype=np.float32),
              np.array([28.65, 45.64, 29.15], dtype=np.float32),
              np.array([28.1 , 46.33, 31.52], dtype=np.float32),
              np.array([28.3837052, 47.6521510, 29.8492298], dtype=np.float32),
              np.array([29.8712052, 46.8618091, 30.4245299], dtype=np.float32))),
        (Data(np.array([78.72, 55.93, 30.95], dtype=np.float32),
              np.array([77.34, 55.74, 31.57], dtype=np.float32),
              np.array([78.85, 55.49, 29.49], dtype=np.float32),
              np.array([79.4352465, 55.3560309, 31.5391384], dtype=np.float32),
              np.array([78.9706437, 56.9893365, 31.0055434], dtype=np.float32))),
        (Data(np.array([15.51, 14.20, 40.00], dtype=np.float32),
              np.array([14.35, 14.44, 39.03], dtype=np.float32),
              np.array([15.67, 15.44, 40.89], dtype=np.float32),
              np.array([16.4287332, 14.0301596, 39.4385875], dtype=np.float32),
              np.array([15.2963201, 13.3288386, 40.6192891], dtype=np.float32))),
    ])
    def test_get_CH2(self, data):
        """Test for get_CH2().
        2 helpers are needed and 2 hydrogen are returned.
        Parameters
        ----------
        data: namedtuple
                structure holding input and reference data.
        """
        assert_almost_equal(hydrogen.get_CH2(data.atom, data.helper1, data.helper2),
                            (data.H1_coord, data.H2_coord))


    # Data tuple for get_CH3() test
    Data = collections.namedtuple('Data', ['atom', 'helper1', 'helper2', 'H1_coord',
                                           'H2_coord', 'H3_coord'])

    @pytest.mark.parametrize('data',[
        (Data(np.array([34.42, 46.94, 26.31], dtype=np.float32),
              np.array([33.40, 46.51, 27.28], dtype=np.float32),
              np.array([32.42, 45.49, 26.87], dtype=np.float32),
              np.array([35.0616142, 47.6932027, 26.7672876], dtype=np.float32),
              np.array([33.9312885, 47.3632873, 25.4324530], dtype=np.float32),
              np.array([35.0224909, 46.0819597, 26.0118850], dtype=np.float32))),
        (Data(np.array([31.05, 43.87, 46.61], dtype=np.float32),
              np.array([30.47, 43.67, 45.21], dtype=np.float32),
              np.array([30.13, 45.04, 44.62], dtype=np.float32),
              np.array([31.2970770, 42.9006050, 47.0428147], dtype=np.float32),
              np.array([30.3152883, 44.3722516, 47.2393194], dtype=np.float32),
              np.array([31.9512322, 44.4797600, 46.5462135], dtype=np.float32))),
        (Data(np.array([63.69, 53.48, 16.38], dtype=np.float32),
              np.array([64.97, 52.74, 16.35], dtype=np.float32),
              np.array([65.30, 52.00, 17.57], dtype=np.float32),
              np.array([63.5437287, 53.9916433, 15.4287233], dtype=np.float32),
              np.array([63.7126355, 54.2125359, 17.1868312], dtype=np.float32),
              np.array([62.8701735, 52.7812533, 16.5465560], dtype=np.float32))),
        (Data(np.array([05.02, 23.08, 20.84], dtype=np.float32),
              np.array([04.46, 24.37, 21.26], dtype=np.float32),
              np.array([03.10, 24.66, 20.78], dtype=np.float32),
              np.array([06.0242550, 22.9693904, 21.2490700], dtype=np.float32),
              np.array([05.0646060, 23.0403350, 19.7516357], dtype=np.float32),
              np.array([04.3870316, 22.2722336, 21.2073760], dtype=np.float32))),
    ])
    def test_get_CH3(self, data):
        """Test for get_CH3().
        2 helpers are needed and 3 hydrogen are returned.
        Parameters
        ----------
        data: namedtuple
                structure holding input and reference data.
        """
        assert_almost_equal(hydrogen.get_CH3(data.atom, data.helper1, data.helper2),
                            (data.H1_coord, data.H2_coord, data.H3_coord))


    # Data tuple for get_CH_double_bond() test
    Data = collections.namedtuple('Data', ['atom', 'helper1', 'helper2', 'H1_coord'])

    @pytest.mark.parametrize('data',[
        (Data(np.array([20.64, 39.80, 31.91], dtype=np.float32),
              np.array([21.13, 41.14, 31.36], dtype=np.float32),
              np.array([20.40, 39.52, 33.25], dtype=np.float32),
              np.array([20.4645443, 38.9986618, 31.1922436], dtype=np.float32))),
        (Data(np.array([15.41, 13.27, 34.47], dtype=np.float32),
              np.array([15.05, 14.68, 34.00], dtype=np.float32),
              np.array([15.71, 12.86, 35.77], dtype=np.float32),
              np.array([15.4351896, 12.4968746, 33.7020546], dtype=np.float32))),
        (Data(np.array([07.29, 24.63, 35.47], dtype=np.float32),
              np.array([07.65, 24.12, 34.08], dtype=np.float32),
              np.array([07.81, 25.78, 36.05], dtype=np.float32),
              np.array([06.5726283, 24.0505400, 36.0511270], dtype=np.float32))),
    ])
    def test_get_CH_double_bond(self, data):
        """Test for get_CH_double_bond().
        2 helpers are needed and 1 hydrogen is returned.
        Parameters
        ----------
        data: namedtuple
                structure holding input and reference data.
        """
        assert_almost_equal(hydrogen.get_CH_double_bond(data.atom, data.helper1, data.helper2),
                            data.H1_coord)
