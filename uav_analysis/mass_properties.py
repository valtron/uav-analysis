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


def quad_copter_props(data: 'TestbenchData') -> Dict[str, sympy.Expr]:
    L0 = sympy.Symbol('Length_0')  # arm
    L1 = sympy.Symbol('Length_1')  # support

    input_data = data.get_tables(['Length_0', 'Length_1'])

    param = 0

    def P():
        nonlocal param
        param += 1
        return sympy.Symbol("param_" + str(param))

    result = dict()

    def fit(name, expr):
        subs, error = approximate(expr, input_data, data.get_table(name))
        print("INFO:", name, "approx error", error)
        result[name] = expr.subs(subs)

    fit('aircraft.mass', P() + P() * L0 + P() * L1)

    fit('aircraft.x_cm',
        (P() + P() * L0 + P() * L1 + P() * L1 ** 2) / result['aircraft.mass'])
    fit('aircraft.y_cm',
        (P() + P() * L0 + P() * L1 + P() * L1 ** 2) / result['aircraft.mass'])
    fit('aircraft.z_cm',
        (P() + P() * L0 + P() * L1 + P() * L1 ** 2) / result['aircraft.mass'])

    fit('aircraft.X_fuseuu', P() + P() * L0 + P() * L1)
    fit('aircraft.Y_fusevv', P() + P() * L0 + P() * L1)
    fit('aircraft.Z_fuseww', P() + P() * L0 + P() * L1)

    fit('aircraft.Ixx', P() + P() * L0 + P() * L1
        + P() * L0 ** 2 + P() * L0 * L1 + P() * L1 ** 2
        + P() * L0 ** 3 + P() * L0 ** 2 * L1 + P() * L0 * L1 ** 2 + P() * L1 ** 3)
    fit('aircraft.Iyy', P() + P() * L0 + P() * L1
        + P() * L0 ** 2 + P() * L0 * L1 + P() * L1 ** 2
        + P() * L0 ** 3 + P() * L0 ** 2 * L1 + P() * L0 * L1 ** 2 + P() * L1 ** 3)
    fit('aircraft.Izz', P() + P() * L0 + P() * L1
        + P() * L0 ** 2 + P() * L0 * L1 + P() * L1 ** 2
        + P() * L0 ** 3 + P() * L0 ** 2 * L1 + P() * L0 * L1 ** 2 + P() * L1 ** 3)

    fit('propeller(1).x', P() * L0)
    fit('propeller(1).y', P() * L0)
    fit('propeller(1).z', P())
    fit('propeller(2).x', P() * L0)
    fit('propeller(2).y', P() * L0)
    fit('propeller(2).z', P())
    fit('propeller(3).x', P() * L0)
    fit('propeller(3).y', P() * L0)
    fit('propeller(3).z', P())
    fit('propeller(4).x', P() * L0)
    fit('propeller(4).y', P() * L0)
    fit('propeller(4).z', P())

    return result


if __name__ == '__main__':
    import sys

    data = TestbenchData()
    data.load(sys.argv[1])

    formulas = quad_copter_props(data)
    print(formulas)
