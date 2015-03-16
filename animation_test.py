#!/usr/bin/env python

import matplotlib
matplotlib.use('TkAgg')
from pylab import *
import numpy as np
#from matplotlib.widgets import Slider, Button, RadioButtons

x = np.arange(0,10,0.01)

for f in np.arange(1,2,0.01):
    y=np.sin(x*f+5*f)
    #print x,f,y
    if f==1:
        l, = plot(x,y)
        xlabel('x')
        ylabel('y')
        show(block=False)
    else:
        l.set_ydata(y)
    pause(0.03)



