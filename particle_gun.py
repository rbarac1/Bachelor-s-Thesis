import numpy as np
import linac_classes as lc
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

#generate N electrons in random directions from the source with an average velocity v_ave and deviation v_dev

# N = 1
# v_ave = 3e4
# v_dev = v_ave/10

#particles = []

linac = lc.accelerator()
 

#create particles with random directions
def create_particles(N, v_ave, v_dev, part="electron", acc=linac):
    for i in range(int(N)):
        v = np.random.normal(v_ave, v_dev)
        theta = np.random.random()*10*np.pi/180
        fi = np.random.random()*360*np.pi/180

        #print(v, theta, fi, v*np.cos(theta))
        #print(v*np.sin(theta)*np.cos(fi), v*np.sin(theta)*np.sin(fi), v*np.cos(theta))
        if part=="electron":
            acc.add_particle(lc.electron([v*np.sin(theta)*np.cos(fi), v*np.sin(theta)*np.sin(fi), v*np.cos(theta)]))
        else:
            acc.add_particle(lc.proton([v*np.sin(theta)*np.cos(fi), v*np.sin(theta)*np.sin(fi), v*np.cos(theta)]))


#find the median velocity of a given collection of particles
def vtest(v_ave, v_dev, N=1e3, part="electron", acc=linac):
    acc.reset()
    create_particles(N, v_ave, v_dev, part, acc=acc)
    return acc.vtest_evolve()


#construct an accelerator by accelerating a single particle
def testrun(vtest, N_electrodes, part="electron", acc=linac):
    if part=="electron":
        acc.add_particle(lc.electron([0, 0, vtest]))
    else:
        acc.add_particle(lc.proton([0, 0, vtest]))
    acc.testrun_evolve(N_electrodes)

#linac.set_electrode_distance(0.7)
linac.add_segment(lc.collimator(0.3, 0.01))
#v = vtest(3e4, 3e1)
#linac.remove_particles()
#print("v=", v)
# linac.remove_particles()
# linac.d.append(0.3)
# create_particles(1, v, v/10)
# linac.evolve()
testrun(3e4, 4)
print("segments=",len(linac.segments))
print("all d", linac.d)
print("beta_f=",linac.beam[0].v[-1,2]/lc.c)

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
