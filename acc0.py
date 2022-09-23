from turtle import color
import accelerator_construction as ac
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


N_part = 30
ac.linac.segments[-1].R = 0.05
ac.create_particles(N_part, 3e4, 3e1, acc=ac.linac)
ac.linac.evolve()
print(ac.linac.beam[0].r[-1,2])

tmax = ac.linac.beam[0].t[-1]
dt = tmax/100
t = 0
# pos = []

# for i in range(N_part):
#     pos.append([])

# while t<tmax:
#     for i in range(N_part):
#         for j in range(len(pg.linac.particles[i].t)):
#             if t>0 and pg.linac.particles[i].t[j-1]<t and pg.linac.particles[i].t[j]>t:
#                 print("i, len(pos), len(linac.particles)", i, len(pos), len(pg.linac.particles))
#                 print("j=", j)
#                 pos[i].append(pg.linac.particles[i].r[j, :])

#     t += dt


# Create lines initially without data
#print("Before x-z", pos)
#particles = [pos[i][:, [2, 0, 1]] for i in range(len(pg.linac.particles))]



particles = [ac.linac.particles[i].r[:, [2, 0, 1]] for i in range(len(ac.linac.particles))]

print("Before", particles)
k = 0       #when an element is deleted, the elements after it are all move by 1
bf = 0  #bugfix variable
bf2 = 0

# for i in range(N_part):
#     k = 0
#     for j in range(len(pg.linac.particles[i].t)-1):
#         t = (j-k)*dt
#         print("t[j], t, t[j+1]",pg.linac.particles[i].t[j], t, pg.linac.particles[i].t[j+1])
#         if (pg.linac.particles[i].t[j]<=t and pg.linac.particles[i].t[j+1]>=t)==False:
#             print("i, len(linac.particles)", i, N_part)
#             print("j=", j)

#             #some weird bugfix
#             if i==1 or bf==1:
#                 bf = 1
#                 if i==0:
#                     bf2 = 1
#                     break

#             print("r = ", particles[i][j-k])
#             particles[i] = np.delete(particles[i], j-k, 0)
#             k += 1
#             #particles[i][j].pop()
        
#         else:
#             print("\n\n\n\n\nt = ", j*dt)
#             print("\n\n\n\n\n")

#     if bf2==1:
#         break

for i in range(N_part):
    k = 0
    for j in range(len(ac.linac.particles[i].t)-1):
        print("t[j], t, t[j+1]",ac.linac.particles[i].t[j], t, ac.linac.particles[i].t[j+1])
        if (j%100==0)==False:
            print("i, len(linac.particles)", i, N_part)
            print("j=", j)

            #some weird bugfix
            if i==1 or bf==1:
                bf = 1
                if i==0:
                    bf2 = 1
                    break

            print("r = ", particles[i][j-k])
            particles[i] = np.delete(particles[i], j-k, 0)
            k += 1
            #particles[i][j].pop()
        
        else:
            print("\n\n\n\n\nt = ", j*dt)
            print("\n\n\n\n\n")

    if bf2==1:
        break


print("\n\n\n\n\n\nAfter", particles)
#particles = [particle_gun.particles[i] for i in range(len(particle_gun.linac.particles))]
#print(particles[0][10])
lines = [ax.plot([], [], [], color="blue")[0] for _ in particles]
#print(particle_gun.particles[0][10])
# Setting the axes properties
ax.set(xlim3d=(0, ac.linac.seg_positions2[-1]+(ac.linac.seg_positions2[-1]-ac.linac.seg_positions1[-1])*0.5), xlabel='Z')
ax.set(ylim3d=(-ac.linac.segments[1].R*1.3, ac.linac.segments[1].R*1.3), ylabel='X')
ax.set(zlim3d=(-ac.linac.segments[1].R*1.3, ac.linac.segments[1].R*1.3), zlabel='Y')


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
    ax.plot_surface(Zc, Xc, Yc, color="pink", alpha=0.5)




#Xc,Yc,Zc = data_for_cylinder_along_z(0.,0.,pg.linac.segments[0].R,pg.linac.segments[0].l)
#ax.plot_surface(Zc, Xc, Yc, alpha=0.5)

ax.set_title("Linac")


#plt.show()

# saving to m4 using ffmpeg writer
writervideo = animation.FFMpegWriter(fps=60)
ani.save(r'/Users/roccobarac/Downloads/first_simulation.mp4', writer=writervideo)
plt.close()