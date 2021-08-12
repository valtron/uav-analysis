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
import numpy

from uav_analysis.testbench_data import TestbenchData
from uav_analysis.func_approx import approximate


def inertia_matrix(ixx, ixy, ixz, iyy, iyz, izz) -> numpy.ndarray:
    """
    Takes 6 values or arrays of shape [*] and returns the inertia tensor 
    matrix of shape [*, 3, 3].
    """
    return numpy.array([[ixx, ixy, ixz], [ixy, iyy, iyz], [ixz, iyz, izz]])


def translate_inertia(
        inertia: numpy.ndarray,
        mass: numpy.ndarray,
        offset: numpy.ndarray) -> numpy.ndarray:
    """
    Returns the moment of inertia matrix with respect to a different point from the
    moment of inertia through the center of mass point. The inertia matrix must be 
    of shape [*, 3, 3], the mass must be of shape [*], and the offset be of shape
    [*, 3], the returned matrix is of shape [*, 3, 3]. If you would like to transfer
    the moment of inertia back to the center of mass, then use negative mass value.
    """
    mat = numpy.array([
        [offset[1] ** 2 + offset[2] ** 2, -offset[0] * offset[1], -offset[0] * offset[2]],
        [-offset[0] * offset[1], offset[0] ** 2 + offset[2] ** 2, -offset[1] * offset[2]],
        [-offset[0] * offset[2], -offset[1] * offset[2], offset[0] ** 2 + offset[1] ** 2]
    ])
    return inertia + mass * mat


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

    if False:
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
        fit('aircraft.Ixy', C() + C() * L0 + C() * L1
            + C() * L0 ** 2 + C() * L0 * L1 + C() * L1 ** 2
            + C() * L0 ** 3 + C() * L0 ** 2 * L1 + C() * L0 * L1 ** 2 + C() * L1 ** 3)
        fit('aircraft.Ixz', C() + C() * L0 + C() * L1
            + C() * L0 ** 2 + C() * L0 * L1 + C() * L1 ** 2
            + C() * L0 ** 3 + C() * L0 ** 2 * L1 + C() * L0 * L1 ** 2 + C() * L1 ** 3)
        fit('aircraft.Iyz', C() + C() * L0 + C() * L1
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

    cm_inertia = inertia_matrix(
        data.get_table('aircraft.Ixx'),
        data.get_table('aircraft.Ixy'),
        data.get_table('aircraft.Ixz'),
        data.get_table('aircraft.Iyy'),
        data.get_table('aircraft.Iyz'),
        data.get_table('aircraft.Izz')
    )
    mass = data.get_table('aircraft.mass')
    offset = numpy.array([
        data.get_table('aircraft.x_cm'),
        data.get_table('aircraft.y_cm'),
        data.get_table('aircraft.z_cm'),
    ])
    rf_inertia = translate_inertia(cm_inertia, mass, offset)

    expr = C() + C() * L0 + C() * L1 + C() * L0 ** 2 + C() * L0 * L1 + C() * L1 ** 2 + \
        C() * L0 ** 3 + C() * L0 ** 2 * L1 + C() * L0 * L1 ** 2 + C() * L1 ** 3
    subs, error = approximate(expr, input_data, rf_inertia[0, 0])
    print("INFO: reference frame Ixx approx error", error)

    return result


if __name__ == '__main__':
    import sys

    data = TestbenchData()
    data.load(sys.argv[1])

    formulas = quad_copter_props(data)
    print(formulas)
