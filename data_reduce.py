from __future__ import division, print_function, absolute_import

import os
import sys
import pickle
import traceback
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import utils # make sure this is from the current directory!

#import matplotlib
#matplotlib.use('TkAgg')


root_dir = '/Users/mikemon/data/Pro_EM/20141121/'

obj = 'SDSSJ0651+2844/SDSSJ0651+2844.spe'
bias = 'bias/bias.spe'
dark = 'dark/dark_10s.spe'
flat = 'dome_flat/dome_flat_BG40_3s.spe'

fjson='reduce_config.json'

com_str = utils.create_reduce_config(rootdir=root_dir,obj=obj,fjson=fjson, bias=bias, dark=dark, flat=flat)

print(com_str)

