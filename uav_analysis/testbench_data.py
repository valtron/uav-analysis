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

from typing import List, Dict

import csv
import os
import sys
import sympy
import numpy
import zipfile

from uav_analysis.fdm_input import parse_fdm_input
from uav_analysis.func_approx import approximate


class TestbenchData():
    def __init__(self):
        self.output_csv = []     # List[Dict[str, str]]
        self.flightdyn_inp = {}  # Dict[str, Dict]

    def load(self, filename: str):
        with zipfile.ZipFile(filename) as file:
            for name in file.namelist():
                if os.path.basename(name) == 'output.csv':
                    with file.open(name) as content:
                        lines = content.readlines()
                        lines = [line.decode('ascii') for line in lines]
                        reader = csv.DictReader(lines)
                        self.output_csv.extend([dict(line) for line in reader])
                elif os.path.basename(name) == 'FlightDyn.inp':
                    guid = os.path.basename(os.path.dirname(name))
                    assert guid not in self.flightdyn_inp
                    with file.open(name) as content:
                        lines = content.readlines()
                        lines = [line.decode('ascii') for line in lines]
                        data = parse_fdm_input(lines)
                        self.flightdyn_inp[guid] = data

    def get_tables(self, fields: List[str]) -> Dict[str, numpy.ndarray]:
        result = {field: [] for field in fields}

        for entry in self.output_csv:
            if entry['AnalysisError'] != 'False':
                continue
            if int(entry['Interferences']) != 0:
                continue

            entry2 = self.flightdyn_inp[entry['GUID']]

            for field in fields:
                if field in entry:
                    value = float(entry[field])
                elif field in entry2:
                    value = float(entry2[field])
                else:
                    raise ValueError("Unknown field " + field)

                result[field].append(value)

        return {key: numpy.array(val) for key, val in result.items()}

    def get_table(self, field: str) -> numpy.ndarray:
        result = self.get_tables([field])
        return result[field]

    def plot2d(self, field1: str, field2: str):
        from matplotlib import pyplot

        tables = self.get_tables([field1, field2])
        fig, ax1 = pyplot.subplots()
        ax1.scatter(
            tables[field1],
            tables[field2],
            s=5.0)
        ax1.set_xlabel(field1)
        ax1.set_ylabel(field2)
        pyplot.show()


if __name__ == '__main__':
    assert len(sys.argv) >= 2
    data = TestbenchData()
    data.load(sys.argv[1])

    data.plot2d('Length_0', 'Length_1')
    data.plot2d('Length_0', 'aircraft.X_fuseuu')
    data.plot2d('Length_1', 'aircraft.X_fuseuu')

    length_0 = sympy.Symbol('Length_0')
    length_1 = sympy.Symbol('Length_1')

    param_a0 = sympy.Symbol('param_a')

    param_b0 = sympy.Symbol('param_b0')
    param_b1 = sympy.Symbol('param_b1')

    param_c0 = sympy.Symbol('param_c0')
    param_c1 = sympy.Symbol('param_c1')
    param_c2 = sympy.Symbol('param_c2')

    param_d0 = sympy.Symbol('param_d0')
    param_d1 = sympy.Symbol('param_d1')
    param_d2 = sympy.Symbol('param_d2')
    param_d3 = sympy.Symbol('param_d3')

    func = param_a0 + param_b0 * length_0 + param_b1 * length_1 \
        + param_c0 * length_0 * length_0 + param_c1 * length_0 * length_1 + param_c2 * length_1 * length_1 \
        + param_d0 * length_0 ** 3 + param_d1 * length_0 ** 2 * length_1 \
        + param_d2 * length_0 * length_1 ** 2 + param_d3 * length_1 ** 3

    input_data = data.get_tables(['Length_0', 'Length_1'])
    print(input_data)

    for name in ['aircraft.mass',
                 'aircraft.x_cm', 'aircraft.y_cm', 'aircraft.z_cm',
                 'aircraft.Ixx', 'aircraft.Iyy', 'aircraft.Izz',
                 'aircraft.X_fuseuu', 'aircraft.Y_fusevv', 'aircraft.Z_fuseww']:
        print("Approximating", name)
        output_data = data.get_table(name)
        approximate(func, input_data, output_data)
