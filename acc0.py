from turtle import color
import particle_gun
import numpy as np

import linac_classes as lc
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D

from IPython import display

# particle_gun.linac.evolve()
# particle_gun.linac.add_segment(lc.collimator(0.5, 0.015))

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


# Create lines initially without data
particles = [particle_gun.particles[i][:, [2, 0, 1]] for i in range(len(particle_gun.linac.particles))]
#particles = [particle_gun.particles[i] for i in range(len(particle_gun.linac.particles))]
#print(particles[0][10])
lines = [ax.plot([], [], [], color="blue")[0] for _ in particles]
#print(particle_gun.particles[0][10])
# Setting the axes properties
ax.set(xlim3d=(0, 1), xlabel='Z')
ax.set(ylim3d=(-0.02, 0.02), ylabel='X')
ax.set(zlim3d=(-0.02, 0.02), zlabel='Y')


# ax.set(xlim3d=(-0.02, 0.02), xlabel='X')
# ax.set(ylim3d=(-0.02, 0.02), ylabel='Y')
# ax.set(zlim3d=(0, 1), zlabel='Z')

# Creating the Animation object
ani = animation.FuncAnimation(fig, update_lines, 300, fargs=(particles, lines), interval=10)

# Make data.
# X = np.linspace(-0.1, 0.1)
# Y = np.linspace(-0.1, 0.1)
# X, Y = np.meshgrid(X, Y)
# R = np.sqrt(X**2 + Y**2)
# Z = np.sin(R)

# ax.plot_surface(X, Y, Z)
#ax.view_init(180,45)


#cylinder plotting
def data_for_cylinder_along_z(center_x,center_y,radius,height_z):
    z = np.linspace(0, height_z, 50)
    theta = np.linspace(0, 2*np.pi, 50)
    theta_grid, z_grid=np.meshgrid(theta, z)
    x_grid = radius*np.cos(theta_grid) + center_x
    y_grid = radius*np.sin(theta_grid) + center_y
    return x_grid,y_grid,z_grid




Xc,Yc,Zc = data_for_cylinder_along_z(0.,0.,particle_gun.linac.segments[0].R,particle_gun.linac.segments[0].l)
ax.plot_surface(Zc, Xc, Yc, alpha=0.5)




plt.show()

# saving to m4 using ffmpeg writer
# writervideo = animation.FFMpegWriter(fps=60)
# ani.save(r'/Users/roccobarac/Downloads/simulation.mp4', writer=writervideo)
# plt.close()