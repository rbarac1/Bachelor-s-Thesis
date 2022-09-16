import numpy as np

import linac_classes as lc

import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D
from IPython import display

def update_lines(num, walks, lines):
    for line, walk in zip(lines, walks):
        # NOTE: there is no .set_data() for 3 dim data...
        line.set_data(walk[:num, :2].T)
        line.set_3d_properties(walk[:num, 2])
    return lines
    
# Attaching 3D axis to the figure
fig = plt.figure()
ax = fig.add_subplot(111,projection="3d")


###########################################################################################################
