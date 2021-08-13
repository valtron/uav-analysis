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

from typing import List, Dict, Any, Union

import csv
import os
import numpy
import zipfile


def parse_fortran_value(value: str) -> Any:
    value = value.strip()
    if value.lower() == '.true.':
        return True
    elif value.lower() == '.false.':
        return False
    elif value.startswith('\''):
        return value.strip('\'')
    elif ',' in value:
        value = value.split(',')
        return [parse_fortran_value(v) for v in value]
    elif 'd' in value or '.' in value:
        return float(value.replace('d', 'e'))
    else:
        return int(value)


def parse_fdm_input(lines: Union[str, List[str]]) -> Dict[str, Any]:
    """
    Reads a new_fdm input file and returns the parameters as a dictionary.
    """
    if isinstance(lines, str):
        with open(lines, 'r') as f:
            lines = f.readlines()

    design = dict()
    state = 0
    for line in lines:
        line = line.strip()
        if line == "&aircraft_data":
            assert state == 0
            state = 1
            continue
        elif line == "/":
            assert state == 1
            state = 2
            continue
        elif state != 1:
            continue

        pos = line.find('!')
        if pos >= 0:
            line = line[:pos].strip()
        if not line:
            continue

        pos = line.find('=')
        assert pos >= 0
        attr = line[:pos].strip().replace('%', '.')
        value = line[pos + 1:].strip()

        design[attr] = parse_fortran_value(value)

    if state != 2:
        print("ERROR: Could not load aircraft design")

    return design


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

    def has_field(self, field: str) -> bool:
        for entry in self.output_csv:
            if entry['AnalysisError'] != 'False':
                continue
            if int(entry['Interferences']) != 0:
                continue

            entry2 = self.flightdyn_inp[entry['GUID']]
            if field not in entry and field not in entry2:
                return False

        return True

    def get_fields(self) -> List[str]:
        fields = set()
        for entry in self.output_csv:
            fields = fields.union(entry.keys())

            entry2 = self.flightdyn_inp[entry['GUID']]
            fields = fields.union(entry2.keys())

        return sorted(list(fields))

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


def run(args=None):
    import argparse

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('file', type=str,  nargs="+", metavar='FILE',
                        help='a zip log files to read')
    parser.add_argument('--plot', type=str, nargs=2, metavar='VAR',
                        help="plots the given pair of values")
    args = parser.parse_args(args)

    data = TestbenchData()
    for file in args.file:
        print("Reading", file)
        data.load(file)

    print("fields:", ','.join(data.get_fields()))

    if args.plot:
        data.plot2d(args.plot[0], args.plot[1])
