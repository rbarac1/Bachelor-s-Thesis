import numpy as np

import linac_classes as lc
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Attaching 3D axis to the figure
fig = plt.figure(figsize=(12,7))
ax = fig.add_subplot(111,projection="3d")



linac = lc.accelerator()

def create_particles(N, v_ave, v_dev, part="electron", acc=linac):
    for i in range(int(N)):
        v = np.random.normal(v_ave, v_dev)
        theta = np.random.random()*20*np.pi/180
        fi = np.random.random()*360*np.pi/180

        #print(v, theta, fi, v*np.cos(theta))
        #print(v*np.sin(theta)*np.cos(fi), v*np.sin(theta)*np.sin(fi), v*np.cos(theta))
        if part=="electron":
            acc.add_particle(lc.electron([v*np.sin(theta)*np.cos(fi), v*np.sin(theta)*np.sin(fi), v*np.cos(theta)]))
        else:
            acc.add_particle(lc.proton([v*np.sin(theta)*np.cos(fi), v*np.sin(theta)*np.sin(fi), v*np.cos(theta)]))

#letting the particles through the accelerator
N_part = 30
linac.add_segment(lc.collimator(0.5, 0.05))
create_particles(N_part, 3e4, 3e1)

linac.d.append(2)
linac.t0 = 100
for part in linac.particles:
    part.starts = 4

#linac.evolve()
#print(linac.beam[0].r[-1,2])

#custom parts
linac.add_segment(lc.electrode(0.2, 0.3))
linac.add_segment(lc.electrode(0.5, 0.3))
linac.add_segment(lc.electrode(0.7, 0.3))
linac.add_segment(lc.collimator(0.5, 0.05))

linac.d = [0.05,0.3,0.6,0.8]


#particle plotting
for i in range(len(linac.particles)):
    x = linac.particles[i].r[:,0]
    y = linac.particles[i].r[:,1]
    z = linac.particles[i].r[:,2]

    #ax.plot(z,x,y, c="b", label="particle" if i==0 else "")


#cylinder plotting
def data_for_cylinder_along_z(center_x,center_y,radius,height_z, z0):
    z = np.linspace(z0, z0+height_z, 50)
    theta = np.linspace(0, 2*np.pi, 50)
    theta_grid, z_grid=np.meshgrid(theta, z)
    x_grid = radius*np.cos(theta_grid) + center_x
    y_grid = radius*np.sin(theta_grid) + center_y
    return x_grid,y_grid,z_grid


for i in range(len(linac.segments)):
    #print("A")
    Xc,Yc,Zc = data_for_cylinder_along_z(0.,0.,linac.segments[i].R,linac.segments[i].l, linac.seg_positions1[i])
    ax.plot_surface(Zc, Xc, Yc, color="red", alpha=0.7)




#plot details
# ax.set(xlim3d=(0, linac.seg_positions2[-1]+(linac.seg_positions2[-1]-linac.seg_positions1[-1])*0.5), xlabel='Z')
# ax.set(ylim3d=(-linac.segments[1].R*1.3, linac.segments[1].R*1.3), ylabel='X')
# ax.set(zlim3d=(-linac.segments[1].R*1.3, linac.segments[1].R*1.3), zlabel='Y')

ax.set(xlim3d=(0, linac.seg_positions2[-1]+(linac.seg_positions2[-1]-linac.seg_positions1[-1])*0.5), xlabel='Z[m]')
ax.set(ylim3d=(-linac.segments[1].R*1.3, linac.segments[1].R*1.3), ylabel='X[m]')
ax.set(zlim3d=(-linac.segments[1].R*1.3, linac.segments[1].R*1.3), zlabel='Y[m]')

#plt.yticks(rotation=20)
ax.yaxis.labelpad = 20
ax.zaxis.labelpad = 10
ax.xaxis.labelpad = 8


ax.xaxis.label.set_size(16)
ax.yaxis.label.set_size(16)
ax.zaxis.label.set_size(16)
#plt.legend()
#plt.title("Collimator test", fontsize=20)
plt.tight_layout()
plt.show()