from __future__ import print_function
import matplotlib
matplotlib.use('TkAgg')
import os
import warnings
import read_spe
import matplotlib.pyplot as plt

file_name = 'SDSSJ0651+2844.spe'
spe = read_spe.File(file_name)

(frame, metadata) = spe.get_frame(-1)
print(metadata)
plt.imshow(frame)
plt.show()
print(metadata)
