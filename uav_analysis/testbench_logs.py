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
import zipfile

from uav_analysis.fdm_input import parse_fdm_input


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


def extract_float_table(testbench_data: Dict, fields: List[str]):
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
                value = entry2['aircraft'][field]
            else:
                print("WARNING: unknown field", field)

            row.append(value)

        print(row)
        table.append(row)

    return table


if __name__ == '__main__':
    testbench_data = read_testbench_zip(sys.argv[1])
    # print(testbench_data['output.csv'][0])
    # print(testbench_data['flightdyn.inp'][testbench_data['output.csv'][0]['GUID']])
    extract_float_table(testbench_data, ['Length_0', 'Length_1',
                                         'aircraft/mass', 'aircraft/x_cm', 'aircraft/y_cm', 'aircraft/z_cm'])
