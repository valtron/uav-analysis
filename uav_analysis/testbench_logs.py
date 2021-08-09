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
import torch
import zipfile

from uav_analysis.fdm_input import parse_fdm_input
from constraint_prog.point_cloud import PointCloud
from constraint_prog.func_approx import approximate


def read_testbench_zip(filename: str) -> Dict:
    testbench_data = {
        'output.csv': None,   # List[Dict[str, str]]
        'flightdyn.inp': {},  # Dict[str, Dict]
    }

    with zipfile.ZipFile(filename) as file:
        for name in file.namelist():
            if os.path.basename(name) == 'output.csv':
                assert testbench_data['output.csv'] is None
                with file.open(name) as content:
                    lines = content.readlines()
                    lines = [line.decode('ascii') for line in lines]
                    reader = csv.DictReader(lines)
                    testbench_data['output.csv'] = [dict(line) for line in reader]
            elif os.path.basename(name) == 'FlightDyn.inp':
                guid = os.path.basename(os.path.dirname(name))
                assert guid not in testbench_data['flightdyn.inp']
                with file.open(name) as content:
                    lines = content.readlines()
                    lines = [line.decode('ascii') for line in lines]
                    data = parse_fdm_input(lines)
                    testbench_data['flightdyn.inp'][guid] = data

    return testbench_data


def extract_float_table(testbench_data: Dict, fields: List[str]) -> numpy.ndarray:
    table = []
    for entry in testbench_data['output.csv']:
        if entry['AnalysisError'] != 'False':
            continue
        if int(entry['Interferences']) != 0:
            continue

        entry2 = testbench_data['flightdyn.inp'][entry['GUID']]

        row = []
        for field in fields:
            if field in entry:
                value = float(entry[field])
            elif field.startswith('aircraft/'):
                field = field[field.find('/') + 1:]
                try:
                    value = entry2['aircraft'][field]
                except KeyError:
                    print("WARNING: unknown field", field)
            else:
                print("WARNING: unknown field", field)

            row.append(value)
        table.append(row)

    return numpy.array(table, dtype=numpy.float64)


if __name__ == '__main__':
    testbench_data = read_testbench_zip(sys.argv[1])
    # print(testbench_data['output.csv'][0])
    # print(testbench_data['flightdyn.inp'][testbench_data['output.csv'][0]['GUID']])

    names = ['Length_0', 'Length_1', 'aircraft/mass',
             'aircraft/x_cm', 'aircraft/y_cm', 'aircraft/z_cm',
             'aircraft/x_fuse', 'aircraft/y_fuse', 'aircraft/z_fuse',
             'aircraft/X_fuseuu', 'aircraft/Y_fusevv', 'aircraft/Z_fuseww',
             'aircraft/Ixx', 'aircraft/Iyy', 'aircraft/Izz']
    table = extract_float_table(testbench_data, names)
    table = torch.tensor(table, dtype=torch.float32)
    points = PointCloud(names, table)

    length_0 = sympy.Symbol('Length_0')
    length_1 = sympy.Symbol('Length_1')

    param_a = sympy.Symbol('param_a')

    param_b0 = sympy.Symbol('param_b0')
    param_b1 = sympy.Symbol('param_b1')

    param_c0 = sympy.Symbol('param_c0')
    param_c1 = sympy.Symbol('param_c1')
    param_c2 = sympy.Symbol('param_c2')

    param_d0 = sympy.Symbol('param_d0')
    param_d1 = sympy.Symbol('param_d1')
    param_d2 = sympy.Symbol('param_d2')
    param_d3 = sympy.Symbol('param_d3')

    func = param_a + param_b0 * length_0 + param_b1 * length_1 \
        + param_c0 * length_0 * length_0 + param_c1 * length_0 * length_1 + param_c2 * length_1 * length_1 \
        + param_d0 * length_0 ** 3 + param_d1 * length_0 ** 2 * length_1 \
        + param_d2 * length_0 * length_1 ** 2 + param_d3 * length_1 ** 3

    approximate(func, points, 'aircraft/Ixx')
