import accelerator_construction_protons as ac
import numpy as np

import linac_classes as lc
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Attaching 3D axis to the figure
fig = plt.figure(figsize=(12,7))
ax = fig.add_subplot(111,projection="3d")


#letting the particles through the accelerator
N_part = 1e4
ac.linac.change_to_collimator(ac.linac.segments[0].l, ac.linac.segments[0].R)
ac.create_particles(N_part, ac.vt, ac.vt/100, part="proton", acc=ac.linac)
print("mass and charge:", ac.linac.beam[0].m, ac.linac.beam[0].q)
ac.linac.evolve()
print(len(ac.linac.beam))


#cylinder plotting
def data_for_cylinder_along_z(center_x,center_y,radius,height_z, z0):
    z = np.linspace(z0, z0+height_z, 50)
    theta = np.linspace(0, 2*np.pi, 50)
    theta_grid, z_grid=np.meshgrid(theta, z)
    x_grid = radius*np.cos(theta_grid) + center_x
    y_grid = radius*np.sin(theta_grid) + center_y
    return x_grid,y_grid,z_grid


for i in range(len(ac.linac.segments)):
    #print("A")
    Xc,Yc,Zc = data_for_cylinder_along_z(0.,0.,ac.linac.segments[i].R,ac.linac.segments[i].l, ac.linac.seg_positions1[i])
    ax.plot_surface(Zc, Xc, Yc, color="blue", alpha=0.6)



#particle plotting
for i in range(len(ac.linac.particles)):
    x = ac.linac.particles[i].r[:,0]
    y = ac.linac.particles[i].r[:,1]
    z = ac.linac.particles[i].r[:,2]

    ax.plot(z,x,y, c="r", lw=0.38 if i==0 else 0.065, label="proton" if i==0 else "")



#plot details
ax.set(xlim3d=(0.2, ac.linac.seg_positions2[-1]+(ac.linac.seg_positions2[-1]-ac.linac.seg_positions1[-1])*0), xlabel='Z[m]')
ax.set(ylim3d=(-ac.linac.segments[1].R*1.3, ac.linac.segments[1].R*1.3), ylabel='X[m]')
ax.set(zlim3d=(-ac.linac.segments[1].R*1.3, ac.linac.segments[1].R*1.3), zlabel='Y[m]')

#rotating the view
ax.view_init(12,-35)

#axis label distances from the numbers
ax.yaxis.labelpad = 20
ax.zaxis.labelpad = 10
ax.xaxis.labelpad = 10

#axis label sizes
ax.xaxis.label.set_size(16)
ax.yaxis.label.set_size(16)
ax.zaxis.label.set_size(16)

plt.legend()
#plt.tight_layout()
plt.title("Simulacija s protonima", fontsize=20)
plt.savefig(r'/Users/roccobarac/Documents/Bachelors_Thesis_files/linac_protoni.png')
plt.title("Simulation with protons", fontsize=20)
plt.savefig(r'/Users/roccobarac/Documents/Bachelors_Thesis_files/linac_protons.png')
#plt.close()
#plt.show()
plt.close()