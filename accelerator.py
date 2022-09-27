from turtle import color
import accelerator_construction_protons as ac
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
fig = plt.figure(figsize=(12,7))
ax = fig.add_subplot(111,projection="3d")
#print("AAAAAAAAAAAAAAAAAAAAAAAAA", type(ax))





# ac.linac.set_electrode_distance(0.7)
# ac.linac.add_segment(lc.collimator(0.3, 0.01))
# ac.linac.add_segment(lc.electrode(0.2, 0.1))



N_part = 1e3

N_part = int(N_part)

ac.linac.change_to_collimator(ac.linac.segments[0].l, ac.linac.segments[0].R)
ac.create_particles(N_part, ac.vt, ac.vt/(1e2), acc=ac.linac)



Nsteps = ac.linac.evolve()
#print(Nsteps)
# Create lines initially without data
particles = [ac.linac.particles[i].r[:, [2, 0, 1]] for i in range(len(ac.linac.particles))]
#particles = []
#particles = [particle_gun.particles[i] for i in range(len(particle_gun.linac.particles))]
#print(particles[0][10])

#print("\n\n\n\n\n\n\n\n\nPARTICLES[0], len PRINT:", particles[0], len(particles[0]))

for i in range(N_part):
    k = 0
    for j in range(len(ac.linac.particles[i].t)-1):
        #print("t[j], t, t[j+1]",ac.linac.particles[i].t[j], t, ac.linac.particles[i].t[j+1])
        if (j%100==0)==False:
            print("i, len(linac.particles)", i, N_part)
            print("j=", j)

            #some weird bugfix
            # if i==1 or bf==1:
            #     bf = 1
            #     if i==0:
            #         bf2 = 1
            #         break

            if len(particles[i])>(j-k):
                print("r = ", particles[i][j-k])
                particles[i] = np.delete(particles[i], j-k, 0)
                k += 1
            else:
                break
        
        #else:
            #print("\n\n\n\n\nt = ", j*dt)
            #print("\n\n\n\n\n")

    # if bf2==1:
    #     break




lines = [ax.plot([], [], [], color="red", linewidth=0.46)[0] for _ in particles]
#print(particle_gun.particles[0][10])



# ax.set(xlim3d=(-0.02, 0.02), xlabel='X')
# ax.set(ylim3d=(-0.02, 0.02), ylabel='Y')
# ax.set(zlim3d=(0, 1), zlabel='Z')

Nsteps = len(particles[0])

for part in particles:
    if len(part)>Nsteps:
        Nsteps = len(part)

Nsteps += -1


# Creating the Animation object
ani = animation.FuncAnimation(fig, update_lines, Nsteps, fargs=(particles, lines), interval=17)


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

# Make data.
# X = np.linspace(-0.1, 0.1)
# Y = np.linspace(-0.1, 0.1)
# X, Y = np.meshgrid(X, Y)
# R = np.sqrt(X**2 + Y**2)
# Z = np.sin(R)

# Setting the axes properties
#plot details
ax.set(xlim3d=(0.2, ac.linac.seg_positions2[-1]+(ac.linac.seg_positions2[-1]-ac.linac.seg_positions1[-1])*0), xlabel='Z[m]')
ax.set(ylim3d=(-ac.linac.segments[1].R*1.3, ac.linac.segments[1].R*1.3), ylabel='X[m]')
ax.set(zlim3d=(-ac.linac.segments[1].R*1.3, ac.linac.segments[1].R*1.3), zlabel='Y[m]')


# ax.plot_surface(X, Y, Z)
# ax.view_init(20,-70)
ax.view_init(12,-35)
ax.set_title("Simulacija sa protonima")


#axis label distances from the numbers
ax.yaxis.labelpad = 20
ax.zaxis.labelpad = 10
ax.xaxis.labelpad = 10

#axis label sizes
ax.xaxis.label.set_size(16)
ax.yaxis.label.set_size(16)
ax.zaxis.label.set_size(16)


#plt.show()

# saving to m4 using ffmpeg writer
writervideo = animation.FFMpegWriter(fps=60)
ani.save(r'/Users/roccobarac/Documents/Bachelors_Thesis_files/new_test.mp4', writer=writervideo)
plt.close()

#print(pg.linac.beam[0].r)