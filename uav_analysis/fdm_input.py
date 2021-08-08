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

import sys


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


def parse_fdm_input(lines: Union[str, List[str]]) -> Dict:
    if isinstance(lines, str):
        with open(lines, 'r') as f:
            lines = f.readlines()

    design = {
        'aircraft': {
            'num_propellers': 0,
            'num_wings': 0,
            'num_batteries': 0
        },
        'control': {},
        'propellers': [],
        'wings': [],
        'batteries': []
    }

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

        pos = line.find('%')
        assert pos >= 0
        part = line[:pos].strip()
        rest = line[pos + 1:]

        pos = rest.find('=')
        assert pos >= 0
        attr = rest[:pos].strip()
        value = rest[pos + 1:].strip()

        if part.endswith(')'):
            pos = part.find('(')
            assert pos >= 0
            index = int(part[pos + 1:-1])
            part = part[:pos]
        else:
            index = -1

        value = parse_fortran_value(value)
        if part == 'aircraft':
            design['aircraft'][attr] = value
        elif part == 'propeller':
            assert index >= 1
            while len(design['propellers']) < index:
                design['propellers'].append(dict())
            design['propellers'][index - 1][attr] = value
        elif part == 'wing':
            assert index >= 1
            while len(design['wings']) < index:
                design['wings'].append(dict())
            design['wings'][index - 1][attr] = value
        elif part == 'battery':
            assert index >= 1
            while len(design['batteries']) < index:
                design['batteries'].append(dict())
            design['batteries'][index - 1][attr] = value
        elif part == 'control':
            design['control'][attr] = value
        else:
            print("WARNING: unexpected property", line)

    if state != 2:
        print("ERROR: Could not load aircraft design")

    for name in ['propellers', 'wings', 'batteries']:
        if design['aircraft']['num_' + name] != len(design[name]):
            print("WARNING: incorrect number of " + name)
        else:
            del design['aircraft']['num_' + name]

    return design


if __name__ == '__main__':
    import json
    design = parse_fdm_input(sys.argv[1])
    print(json.dumps(design, indent=2))
