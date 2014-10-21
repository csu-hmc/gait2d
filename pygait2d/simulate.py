#!/usr/bin/env python
# -*- coding: utf-8 -*-

from collections import OrderedDict

import yaml


def load_constants(constants, path):
    """Parses a yaml file and builds an ordered dictionary that maps SymPy
    symbols to floats."""

    with open(path, 'r') as f:
        constant_values_dict = yaml.load(f)

    res = OrderedDict()

    for c in constants:
        res[c] = constant_values_dict[c.name]

    return res


def map_values_to_autolev_symbols(constants):
    """Returns a dictionary mapping the autoleve symbol names to the ones
    used in Python model.

    Parameters
    ==========
    constants : dictionary
        Maps python symbol names to floats.

    Returns
    =======
    d : dictionary
        Maps autolev symbol names to floats.

    """

    d = {}
    d['TrunkMass'] = constants['ma']
    d['TrunkInertia'] = constants['ia']
    d['TrunkCMy'] = constants['ya']
    d['ThighMass'] = constants['mb']
    d['ThighInertia'] = constants['ib']
    d['ThighCMy'] = constants['yb']
    d['ThighLen'] = constants['lb']
    d['ShankMass'] = constants['mc']
    d['ShankInertia'] = constants['ic']
    d['ShankCMy'] = constants['yc']
    d['ShankLen'] = constants['lc']
    d['FootMass'] = constants['md']
    d['FootInertia'] = constants['id']
    d['FootCMx'] = constants['xd']
    d['FootCMy'] = constants['yd']
    d['ContactY'] = constants['fyd']
    d['ContactHeelX'] = constants['hxd']
    d['ContactToeX'] = constants['txd']
    d['ContactStiff'] = constants['kc']
    d['ContactDamp'] = constants['cc']
    d['ContactV0'] = constants['vs']
    d['ContactFric'] = constants['mu']

    return d
