# -*- coding: utf-8 -*-
"""
Created on Fri Mar  8 09:31:31 2024

@author: fentonlab
"""

import sys
import os 

dir_current = os.path.dirname(__file__)
print("CURRENT DIR", dir_current)
dir_utils = os.path.abspath(os.path.join(dir_current, '..', 'utils-tools'))

print("UTILS DIR", dir_utils)
sys.path.append(dir_utils)

from utils_phase_frequency import *
print(utils_phase_frequency.__file__)