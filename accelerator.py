from turtle import color
import particle_gun as pg
import numpy as np

import linac_classes as lc
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D

from IPython import display

#print(particle_gun.particles)
def update_lines(num, walks, lines):
    for line, walk in zip(lines, walks):
        # NOTE: there is no .set_data() for 3 dim data...
        line.set_data(walk[:num, :2].T)
        line.set_3d_properties(walk[:num, 2])
    return lines
    
# Attaching 3D axis to the figure
fig = plt.figure()
ax = fig.add_subplot(111,projection="3d")
#print("AAAAAAAAAAAAAAAAAAAAAAAAA", type(ax))





pg.linac.set_electrode_distance(0.5)
pg.linac.add_segment(lc.collimator(0.5, 0.015))
pg.linac.add_segment(lc.electrode(1.0, 0.1))





Nsteps = pg.linac.evolve()
#print(Nsteps)
# Create lines initially without data
particles = [pg.linac.particles[i].r[:, [2, 0, 1]] for i in range(len(pg.linac.particles))]
#particles = [particle_gun.particles[i] for i in range(len(particle_gun.linac.particles))]
#print(particles[0][10])
lines = [ax.plot([], [], [], color="blue")[0] for _ in particles]
#print(particle_gun.particles[0][10])
# Setting the axes properties
ax.set(xlim3d=(0, pg.linac.seg_positions2[-1]+(pg.linac.seg_positions2[-1]-pg.linac.seg_positions1[-1])*0.5), xlabel='Z')
xy_max = pg.linac.segments[-1].R*1.35
ax.set(ylim3d=(-xy_max, xy_max), ylabel='X')
ax.set(zlim3d=(-xy_max, xy_max), zlabel='Y')


# ax.set(xlim3d=(-0.02, 0.02), xlabel='X')
# ax.set(ylim3d=(-0.02, 0.02), ylabel='Y')
# ax.set(zlim3d=(0, 1), zlabel='Z')

# Creating the Animation object
ani = animation.FuncAnimation(fig, update_lines, Nsteps, fargs=(particles, lines), interval=100)


#cylinder plotting
def data_for_cylinder_along_z(center_x,center_y,radius,height_z, z0):
    z = np.linspace(z0, z0+height_z, 50)
    theta = np.linspace(0, 2*np.pi, 50)
    theta_grid, z_grid=np.meshgrid(theta, z)
    x_grid = radius*np.cos(theta_grid) + center_x
    y_grid = radius*np.sin(theta_grid) + center_y
    return x_grid,y_grid,z_grid



for i in range(len(pg.linac.segments)):
    #print("A")
    Xc,Yc,Zc = data_for_cylinder_along_z(0.,0.,pg.linac.segments[i].R,pg.linac.segments[i].l, pg.linac.seg_positions1[i])
    ax.plot_surface(Zc, Xc, Yc, color="red", alpha=0.5)

# Make data.
# X = np.linspace(-0.1, 0.1)
# Y = np.linspace(-0.1, 0.1)
# X, Y = np.meshgrid(X, Y)
# R = np.sqrt(X**2 + Y**2)
# Z = np.sin(R)

# ax.plot_surface(X, Y, Z)
ax.view_init(20,-70)
ax.set_title("Linac")

#plt.show()

# saving to m4 using ffmpeg writer
writervideo = animation.FFMpegWriter(fps=60)
ani.save(r'/Users/roccobarac/Documents/Bachelors_Thesis_files/simulation.mp4', writer=writervideo)
plt.close()

#print(pg.linac.beam[0].r)