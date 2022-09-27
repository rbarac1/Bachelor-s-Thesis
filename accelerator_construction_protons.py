import numpy as np
import linac_classes as lc
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


################################################################
#This file is used to construct an Alvarez type linac given the type of particles and it's initial velocities
################################################################


# N = 1
# v_ave = 3e4
# v_dev = v_ave/10

#particles = []

linac = lc.accelerator()
 

#create N particles with random directions
def create_particles(N, v_ave, v_dev, part="proton", acc=linac, C=0, A=0):
    for i in range(int(N)):
        v = np.random.normal(v_ave, v_dev)     #speed
        theta = np.random.random()*10*np.pi/180     #polar angle
        fi = np.random.random()*360*np.pi/180       #azimuthal angle

        vx = v*np.sin(theta)*np.cos(fi)
        vy = v*np.sin(theta)*np.sin(fi)
        vz = v*np.cos(theta)

        if part=="proton":
            acc.add_particle(lc.proton([vx, vy, vz]))
        elif part=="ion":
            acc.add_particle(lc.ion([vx, vy, vz],C,A))
        else:
            acc.add_particle(lc.electron([vx, vy, vz], C, A))



#find the median velocity of a given collection of particles
def vtest(v_ave, v_dev, N=1e3, part="proton", acc=linac, C=0, A=0):
    acc.reset()
    create_particles(N, v_ave, v_dev, part, acc=acc, C=C, A=A)
    return acc.vtest_evolve()


#construct an accelerator by accelerating a single particle
def testrun(vtest, N_electrodes, part="proton", acc=linac, C=0, A=0):
    if part=="proton":
        acc.add_particle(lc.proton([0, 0, vtest]))
    elif part=="ion":
        acc.add_particle(lc.ion([0, 0, vtest],C,A))
    else:
        acc.add_particle(lc.electron([0, 0, vtest]))
    acc.testrun_evolve(N_electrodes)

#just creates a test particle
def create_test_particle(vtest, part="proton", acc=linac, C=0, A=0):
    if part=="proton":
        acc.add_particle(lc.proton([0, 0, vtest]))
    elif part=="ion":
        acc.add_particle(lc.ion([0, 0, vtest],C,A))
    else:
        acc.add_particle(lc.electron([0, 0, vtest]))

#linac.set_electrode_distance(0.7)
linac.add_segment(lc.collimator(0.5, 0.01))
vt = 0.08*lc.c
v = vt
v = vtest(vt, vt/(1e2))
#print(len(linac.beam))
linac.remove_particles()
#print(len(linac.beam))
#print("v=", v/lc.c)
# linac.d.append(0.3)
#create_particles(1, v, 0)
#print(len(linac.beam))
# linac.evolve()
linac.set_electrode_radius(0.2)
testrun(v, 10)
print("segments=",len(linac.segments))
print("all d", linac.d)
print("beta_0=", v/lc.c)
print("beta_f=",linac.beam[0].v[-1,2]/lc.c)
print("first collimator escape time: ", linac.t0)
print("total time: ", linac.beam[0].t[-1])
print("Energies at the end of each segment in MeV: ", linac.beam[0].E_MeV)
linac.remove_particles()
linac.change_to_collimator(linac.segments[0].l, linac.segments[0].R)

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
