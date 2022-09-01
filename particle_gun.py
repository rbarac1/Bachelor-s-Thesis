import numpy as np
import linac_classes as lc
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

#generate N electrons in random directions from the source with an average velocity v_ave and deviation v_dev

N = 1
v_ave = 3e1
v_dev = v_ave/10

#particles = []

linac = lc.accelerator()

#electrons
for i in range(int(N)):
    v = np.random.normal(v_ave, v_dev)
    theta = np.random.random()*10*np.pi/180*0
    fi = np.random.random()*360*np.pi/180*0

    #print(v, theta, fi, v*np.cos(theta))
    #print(v*np.sin(theta)*np.cos(fi), v*np.sin(theta)*np.sin(fi), v*np.cos(theta))
    linac.add_particle(lc.electron([v*np.sin(theta)*np.cos(fi), v*np.sin(theta)*np.sin(fi), v*np.cos(theta)], dt=1e-4))


#collimator


#print(linac.particles[0].r)
#print(linac.seg_positions)
#linac.add_segment(lc.collimator(0.5, 0.015))
#linac.evolve()
# for i in range(int(N)):
#     print(linac.particles[i].v[0])
# print("positions")
# for i in range(int(N)):
#     print(linac.particles[i].r[10])
#linac.evolve()
#print(type(particles[0])==type(part.electron(10)))
#print(linac.particles[0].r[-1])
#print(linac.particles[0].r[-1,2])
#print(len(linac.particles), len(linac.beam))

#for i in range(len(linac.particles)):
    #plt.plot(linac.particles[i].t, np.sqrt(linac.particles[i].r[:,0]**2+linac.particles[i].r[:,1]**2))
    #plt.plot(linac.particles[i].r[:,2], np.sqrt(linac.particles[i].r[:,0]**2+linac.particles[i].r[:,1]**2))
#plt.show()

#particles = [linac.particles[i].r for i in range(len(linac.particles))]

# print("list positions")
# for i in range(int(N)):
#     print(particles[i][10])
