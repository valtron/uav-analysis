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

from typing import Dict

import sympy

from uav_analysis.testbench_data import TestbenchData
from uav_analysis.mass_properties import quad_copter_fixed_bemp2


def get_formulas1(filename: str) -> Dict[str, sympy.Expr]:
    data = TestbenchData()
    data.load('/home/mmaroti/Downloads/MiklosConf2.zip')
    return quad_copter_fixed_bemp2(data)


def get_formulas2() -> Dict[str, sympy.Expr]:
    Length_0 = sympy.Symbol('Length_0')  # arm length
    Length_1 = sympy.Symbol('Length_1')  # support leg
    Length_8 = sympy.Symbol('Length_8')  # battery offset
    Length_9 = sympy.Symbol('Length_9')  # battery offset

    return {'aircraft.mass': 0.000227122507201286*Length_0 + 0.000227122507201246*Length_1 + 2.22035363865589, 'aircraft.x_cm': (0.806101730552653*Length_8 + 0.806101730552652*Length_9 + 0.00618926278194512)/(0.000227122507201286*Length_0 + 0.000227122507201246*Length_1 + 2.22035363865589), 'aircraft.y_cm': (0.806101730552652*Length_8 - 0.806101730552651*Length_9 + 0.00618853054686912)/(0.000227122507201286*Length_0 + 0.000227122507201246*Length_1 + 2.22035363865589), 'aircraft.z_cm': (0.00230756467316657*Length_0 + 0.000113561253600658*Length_1**2 + 0.00346134700975123*Length_1 - 63.5384903704215)/(0.000227122507201286*Length_0 + 0.000227122507201246*Length_1 + 2.22035363865589), 'aircraft.Ixx': 3.78537512001683e-5*Length_0**3 + 0.000113561253600669*Length_0**2*Length_1 + 0.524693277856849*Length_0**2 - 6.43892518904808*Length_0 + 7.57075024002329e-5*Length_1**3 + 0.00346134700984928*Length_1**2 + 0.0549322159353924*Length_1 + 0.569999999999982*Length_8**2 - 1.13999999999999*Length_8*Length_9 + 0.00875189910994965*Length_8 + 0.570000000000005*Length_9**2 - 0.00875189911019538*Length_9 - (0.131806068394536*(Length_8 - 0.999999999999999*Length_9 + 0.00767710862328301)**2/(0.000102291141035883*Length_0 + 0.000102291141035865*Length_1 + 1)**2 + 818.897382477748*(3.63175873350744e-5*Length_0 + 1.78728284129211e-6*Length_1**2 + 5.44763810026333e-5*Length_1 - 1)**2/(0.000102291141035883*Length_0 + 0.000102291141035865*Length_1 + 1)**2)*(0.000227122507201286*Length_0 + 0.000227122507201246*Length_1 + 2.22035363865589) + 6850.25986589179, 'aircraft.Iyy': 3.78537512002041e-5*Length_0**3 + 0.000113561253600664*Length_0**2*Length_1 + 0.524693277856803*Length_0**2 - 6.43892518902789*Length_0 + 7.57075024003283e-5*Length_1**3 + 0.0034613470098049*Length_1**2 + 0.0549322159428197*Length_1 + 0.569999999999982*Length_8**2 + 1.14*Length_8*Length_9 + 0.00875292998125496*Length_8 + 0.57*Length_9**2 + 0.0087529299814994*Length_9 - (0.131806068394536*(Length_8 + 1.0*Length_9 + 0.00767801698887539)**2/(0.000102291141035883*Length_0 + 0.000102291141035865*Length_1 + 1)**2 + 818.897382477748*(3.63175873350744e-5*Length_0 + 1.78728284129211e-6*Length_1**2 + 5.44763810026333e-5*Length_1 - 1)**2/(0.000102291141035883*Length_0 + 0.000102291141035865*Length_1 + 1)**2)*(0.000227122507201286*Length_0 + 0.000227122507201246*Length_1 + 2.22035363865589) + 6850.2599397045, 'aircraft.Izz': 7.57075024003728e-5*Length_0**3 + 0.000227122507201336*Length_0**2*Length_1 + 1.04938655571365*Length_0**2 - 12.9291022697404*Length_0 - 3.27264989944383e-16*Length_1**3 + 1.5470692105743e-13*Length_1**2 + 0.00436257502198411*Length_1 + 1.13999999999997*Length_8**2 - 9.66601953571952e-15*Length_8*Length_9 + 0.0175048290907453*Length_8 + 1.13999999999999*Length_9**2 + 1.030871682827e-6*Length_9 - (0.131806068394536*(Length_8 - 0.999999999999999*Length_9 + 0.00767710862328301)**2/(0.000102291141035883*Length_0 + 0.000102291141035865*Length_1 + 1)**2 + 0.131806068394536*(Length_8 + 1.0*Length_9 + 0.00767801698887539)**2/(0.000102291141035883*Length_0 + 0.000102291141035865*Length_1 + 1)**2)*(0.000227122507201286*Length_0 + 0.000227122507201246*Length_1 + 2.22035363865589) + 8148.94905326617, 'aircraft.Ixy': -9.85631692595163e-18*Length_0**3 + 1.3029845074765e-19*Length_0**2*Length_1 + 1.34995350435647e-14*Length_0**2 - 6.08612325450176e-12*Length_0 - 4.53764077099293e-18*Length_1**3 + 2.43930311835976e-15*Length_1**2 - 4.35460430237572e-13*Length_1 - 0.569999999999995*Length_8**2 + 1.69889567496315e-15*Length_8*Length_9 - 0.00875241454584585*Length_8 + 0.569999999999996*Length_9**2 + 5.15435785055202e-7*Length_9 - 1264.81806045157 + (0.806101730552652*Length_8 - 0.806101730552651*Length_9 + 0.00618853054686912)*(0.806101730552653*Length_8 + 0.806101730552652*Length_9 + 0.00618926278194512)/(0.000227122507201286*Length_0 + 0.000227122507201246*Length_1 + 2.22035363865589), 'aircraft.Ixz': -8.44694021480435e-19*Length_0**3 - 3.6801770862108e-20*Length_0**2*Length_1 + 1.11668973057976e-15*Length_0**2 - 4.64107113368816e-13*Length_0 + 4.23216264270344e-19*Length_1**3 - 1.17167790260786e-16*Length_1**2 + 1.76803213825051e-14*Length_1 + 2.72289014600044e-16*Length_8**2 + 5.99821390014281e-16*Length_8*Length_9 + 30.3677962727528*Length_8 + 5.67114146985124e-17*Length_9**2 + 30.3677962727528*Length_9 + 0.0506073323905293 + (0.806101730552653*Length_8 + 0.806101730552652*Length_9 + 0.00618926278194512)*(0.00230756467316657*Length_0 + 0.000113561253600658*Length_1**2 + 0.00346134700975123*Length_1 - 63.5384903704215)/(0.000227122507201286*Length_0 + 0.000227122507201246*Length_1 + 2.22035363865589), 'aircraft.Iyz': 2.15664972175755e-18*Length_0**3 - 8.93168512177589e-20*Length_0**2*Length_1 - 2.85858219450388e-15*Length_0**2 + 1.23829231924835e-12*Length_0 + 1.40768808860066e-18*Length_1**3 - 8.29361988670351e-16*Length_1**2 + 1.65109327869914e-13*Length_1 + 6.10613902735457e-16*Length_8**2 - 4.58450634456839e-16*Length_8*Length_9 + 30.3677962727528*Length_8 + 1.06762992767642e-16*Length_9**2 - 30.3677962727528*Length_9 + 0.0505971853390777 + (0.806101730552652*Length_8 - 0.806101730552651*Length_9 + 0.00618853054686912)*(0.00230756467316657*Length_0 + 0.000113561253600658*Length_1**2 + 0.00346134700975123*Length_1 - 63.5384903704215)/(0.000227122507201286*Length_0 + 0.000227122507201246*Length_1 + 2.22035363865589), 'aircraft.X_fuseuu': 0.000234581361221017*Length_0**3 - 0.369566254724661*Length_0**2 + 206.215520990621*Length_0 - 0.000204946520509182*Length_1**3 - 0.0173268238789224*Length_1**2 + 32.2428262451149*Length_1 - 24217.1455819846, 'aircraft.Y_fusevv': 0.000173457176201281*Length_0**3 - 0.286003650043083*Length_0**2 + 168.247280840272*Length_0 - 0.000433334445471943*Length_1**3 + 0.106381792222529*Length_1**2 + 13.0019500935199*Length_1 - 17659.9994636514, 'aircraft.Z_fuseww': -0.000592348683105564*Length_0**3 + 0.839357563849155*Length_0**2 - 351.00958903508*Length_0 - 0.00154063159082293*Length_1**3 + 0.73678633744287*Length_1**2 - 109.096063879292*Length_1 + 60452.0986740916, 'Prop_0_x': -0.707106781186547*Length_0, 'Prop_0_y': -0.707106781186547*Length_0, 'Prop_0_z': -35.5616666666667, 'Prop_1_x': 0.707106781186547*Length_0, 'Prop_1_y': -0.707106781186547*Length_0, 'Prop_1_z': -35.5616666666667, 'Prop_2_x': 0.707106781186547*Length_0, 'Prop_2_y': 0.707106781186547*Length_0, 'Prop_2_z': -35.5616666666667, 'Prop_3_x': -0.707106781186547*Length_0, 'Prop_3_y': 0.707106781186547*Length_0, 'Prop_3_z': -35.5616666666666}


def get_params(formulas: Dict[str, sympy.Expr], subs: Dict[str, float]) -> Dict[str, float]:
    result = dict()
    for key, val in formulas.items():
        if isinstance(val, float):
            result[key] = val
        else:
            result[key] = float(val.subs(subs))
    return result


if __name__ == '__main__':
    formulas = get_formulas2()
    subs = {'Length_0': 100, 'Length_1': 50, 'Length_8': 0, 'Length_9': 0}
    print(get_params(formulas, subs))
