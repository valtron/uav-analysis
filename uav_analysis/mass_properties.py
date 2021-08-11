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

import sympy

from uav_analysis.testbench_data import TestbenchData
from uav_analysis.func_approx import approximate


def quad_copter(data: 'TestbenchData') -> Dict[str, sympy.Expr]:
    length_0 = sympy.Symbol('Length_0')
    length_1 = sympy.Symbol('Length_1')

    input_data = data.get_tables(['Length_0', 'Length_1'])

    param_idx = 0

    def param():
        nonlocal param_idx
        param_idx += 1
        return sympy.Symbol("param_" + str(param_idx))

    mass = param() + param() * length_0 + param() * length_1
    subs, _ = approximate(mass, input_data, data.get_table('aircraft.mass'))
    print(subs)
    mass = mass.subs(subs)

    return {
        'mass': mass
    }


if __name__ == '__main__':
    import sys

    data = TestbenchData()
    data.load(sys.argv[1])

    print(quad_copter(data))
