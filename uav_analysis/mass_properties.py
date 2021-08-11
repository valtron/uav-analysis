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

    def C():
        nonlocal param
        param += 1
        return sympy.Symbol("param_" + str(param))

    result = dict()

    def fit(name, expr):
        subs, error = approximate(expr, input_data, data.get_table(name))
        print("INFO:", name, "approx error", error)
        result[name] = expr.subs(subs)

    fit('aircraft.mass', C() + C() * L0 + C() * L1)

    fit('aircraft.x_cm', C() / result['aircraft.mass'])
    fit('aircraft.y_cm', C() / result['aircraft.mass'])
    fit('aircraft.z_cm', (C() + C() * L0 + C() * L1 + C() * L1 ** 2) / result['aircraft.mass'])

    fit('aircraft.X_fuseuu', C() + C() * L0 + C() * L1)
    fit('aircraft.Y_fusevv', C() + C() * L0 + C() * L1)
    fit('aircraft.Z_fuseww', C() + C() * L0 + C() * L1)

    fit('aircraft.Ixx', C() + C() * L0 + C() * L1
        + C() * L0 ** 2 + C() * L0 * L1 + C() * L1 ** 2
        + C() * L0 ** 3 + C() * L0 ** 2 * L1 + C() * L0 * L1 ** 2 + C() * L1 ** 3)
    fit('aircraft.Iyy', C() + C() * L0 + C() * L1
        + C() * L0 ** 2 + C() * L0 * L1 + C() * L1 ** 2
        + C() * L0 ** 3 + C() * L0 ** 2 * L1 + C() * L0 * L1 ** 2 + C() * L1 ** 3)
    fit('aircraft.Izz', C() + C() * L0 + C() * L1
        + C() * L0 ** 2 + C() * L0 * L1 + C() * L1 ** 2
        + C() * L0 ** 3 + C() * L0 ** 2 * L1 + C() * L0 * L1 ** 2 + C() * L1 ** 3)

    fit('propeller(1).x', C() * L0)
    fit('propeller(1).y', C() * L0)
    fit('propeller(1).z', C())
    fit('propeller(2).x', C() * L0)
    fit('propeller(2).y', C() * L0)
    fit('propeller(2).z', C())
    fit('propeller(3).x', C() * L0)
    fit('propeller(3).y', C() * L0)
    fit('propeller(3).z', C())
    fit('propeller(4).x', C() * L0)
    fit('propeller(4).y', C() * L0)
    fit('propeller(4).z', C())

    return result


if __name__ == '__main__':
    import sys

    data = TestbenchData()
    data.load(sys.argv[1])

    formulas = quad_copter_props(data)
    print(formulas)
