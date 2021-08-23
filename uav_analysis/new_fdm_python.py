#!/usr/bin/env python3
# Copyright (C) 2021, Miklos Maroti
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

from typing import List, Tuple

import numpy
import os

PROP_DIR = os.path.join(os.path.dirname(__file__), 'data', 'propeller')


class Propeller(object):
    def __init__(self, cname: str):
        self.cname = cname

    @staticmethod
    def propread(filename: str) -> Tuple[List[float], List[float], numpy.ndarray]:
        """
        This method reads the propeller file and returns two lists (the speed
        of the propeller in radians per second, and the advance ratio), and
        a rectangular table of shape (len(vs), len(js), 2) containing the Ct and
        Cp values for each combination of speeds and advance ratios. This 
        implementation is based on `subroutine propread` in `propread.f`.
        """

        # we first read the data into a Dict[float, List[Tuple[float, float, float]]]
        data = dict()
        with open(filename, 'r') as file:
            rpm = None
            for line in file.readlines():
                line = line.strip()
                if line.startswith("PROP RPM = "):
                    line = line.split()
                    rpm = float(line[3])
                    assert rpm not in data
                    data[rpm] = []
                elif rpm and line and line[:1].isdigit():
                    line = line.split()
                    j = float(line[1])
                    ct = float(line[3])
                    cp = float(line[4])
                    assert not data[rpm] or data[rpm][-1][0] < j
                    data[rpm].append((j, ct, cp))

        # convert the data into a Dict[float, numpy.ndarray]
        data = {key: numpy.array(val, dtype=numpy.float64)
                for key, val in data.items()}

        # the original code interpolates this to a common shape
        rpms = numpy.array(sorted(data.keys()), dtype=numpy.float64)
        js = data[rpms[0]][:, 0]
        assert list(js) == sorted(js)

        table = numpy.empty((len(rpms), len(js), 2), dtype=numpy.float64)
        table[0] = data[rpms[0]][:, 1:]

        for idx, rpm in enumerate(rpms):
            if idx == 0:
                continue
            table[idx, :, 0] = numpy.interp(js, data[rpm][:, 0], data[rpm][:, 1])
            table[idx, :, 1] = numpy.interp(js, data[rpm][:, 0], data[rpm][:, 2])

        omegas = rpms * (2.0 * numpy.pi / 60.0)
        return omegas, js, table


if __name__ == '__main__':
    omegas, js, table = Propeller.propread(
        os.path.join(PROP_DIR, 'PER3_875x50.dat'))
