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
