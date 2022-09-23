import accelerator_construction as ac
import numpy as np

import linac_classes as lc
import matplotlib.pyplot as plt

print("\n\nfirst gap distance: ", 1000*((ac.linac.seg_positions2[1]+ac.linac.d[1]/2)-(ac.linac.seg_positions2[0]+ac.linac.d[0]/2)))
print("last gap distance: ", 1000*((ac.linac.seg_positions2[-2]+ac.linac.d[-1]/2)-(ac.linac.seg_positions2[-3]+ac.linac.d[-2]/2)))



N_part = 100

#ac.linac.change_to_collimator(ac.linac.segments[0].l, ac.linac.segments[0].R)

#l_beam = 0.05

#print(int(N_part))

for i in range(int(N_part)):
    ac.create_test_particle(ac.v, part="ion", C=-1, A=1)
    #ac.linac.beam[-1].r = np.array([[0,0,l_beam*(i+1)/2/N_part]])
    #ac.linac.beam[-1].t = [2*l_beam*(i+1)/2/N_part/ac.v]
    #ac.linac.beam[-1].r = np.append(ac.linac.beam[-1].r, np.array([0,0,l_beam*(i+1)/2/N_part]),0)
    #print("TEST")
    #print(l_beam*(i+1)/2/N_part, ac.linac.beam[-1].r)
    ac.linac.beam[-1].v[0,2] = np.random.normal(ac.v, ac.v/(1e4))

#ac.linac.t0 += (l_beam/2)/ac.v
#print(ac.linac.t0)
ac.linac.evolve()

for part in ac.linac.beam:
    print(part.r[-1, 2], part.v[-1, 2])


print("\n\nFIRST COLLIMATOR:")
print("N = ", len(ac.linac.v1))
v1 = np.mean(ac.linac.v1)
v1std = np.std(ac.linac.v1)
E1 = np.mean(ac.linac.E1)
E1std = np.std(ac.linac.E1)
print("beta, beta_std: ", np.mean(ac.linac.v1), np.std(ac.linac.v1))
print("E, E_std: ", np.mean(ac.linac.E1), np.std(ac.linac.E1))

print("\n\nSECOND COLLIMATOR:")
print("N = ", len(ac.linac.v2))
print("beta, beta_std: ", np.mean(ac.linac.v2), np.std(ac.linac.v2))
print("E, E_std: ", np.mean(ac.linac.E2), np.std(ac.linac.E2))

print(ac.linac.beam[0].r[:,2])

v2 = np.mean(ac.linac.v2)
v2std = np.std(ac.linac.v2)

E2 = np.mean(ac.linac.E2)
E2std = np.std(ac.linac.E2)

for i in range(len(ac.linac.v1)):
    if np.abs(ac.linac.v1[i]-v1)>5*v1std:
        del ac.linac.v1[i]
        del ac.linac.E1[i]

for i in range(len(ac.linac.v2)):
    if np.abs(ac.linac.v2[i]-v2)>5*v2std:
        del ac.linac.v2[i]
        del ac.linac.E2[i]


#1st collimator plot
plt.figure(figsize=(12,7))

#1st collimator velocity
plt.subplot(1,2,1)
plt.hist(ac.linac.v1, 50, facecolor='g', alpha=0.85)
plt.xlabel(r'$\beta$', fontsize=14)
plt.ylabel('Number of particles', fontsize=14)
plt.title('Speed histogram', fontsize=16)
plt.xlim(v1-4*v1std, v1+4*v1std)
#plt.ylim(0, 0.03)
plt.grid(True)


#1st collimator energy
plt.subplot(1,2,2)
plt.hist(ac.linac.E1, 50, facecolor='g', alpha=0.85)
plt.xlabel("E[MeV]", fontsize=14)
plt.ylabel('Number of particles', fontsize=14)
plt.title('Energy histogram', fontsize=16)
plt.xlim(E1-4*E1std, E1+4*E1std)
#plt.ylim(0, 0.03)
plt.grid(True)

plt.tight_layout()
plt.suptitle("1st collimator", fontsize=20)
plt.subplots_adjust(top=0.88)
plt.savefig(r'/Users/roccobarac/Documents/Bachelors_Thesis_files/linac4-1mod.png')
plt.close()




#2nd collimator plot
plt.figure(figsize=(12,7))

#2nd collimator velocity
plt.subplot(1,2,1)
plt.hist(ac.linac.v2, 50, facecolor='g', alpha=0.85)
plt.xlabel(r'$\beta$', fontsize=14)
plt.ylabel('Number of particles', fontsize=14)
plt.title('Speed histogram', fontsize=16)
plt.xlim(v2-4*v2std, v2+4*v2std)
#plt.ylim(0, 0.03)
plt.grid(True)


#2nd collimator energy
plt.subplot(1,2,2)
plt.hist(ac.linac.E2, 50, facecolor='g', alpha=0.85)
plt.xlabel("E[MeV]", fontsize=14)
plt.ylabel('Number of particles', fontsize=14)
plt.title('Energy histogram', fontsize=16)
plt.xlim(E2-4*E2std, E2+4*E2std)
#plt.ylim(0, 0.03)
plt.grid(True)

plt.tight_layout()
plt.suptitle("2nd collimator", fontsize=20)
plt.subplots_adjust(top=0.88)
plt.savefig(r'/Users/roccobarac/Documents/Bachelors_Thesis_files/linac4-2mod.png')
plt.close()







#######################################################
#croatian plots
#1st collimator plot
plt.figure(figsize=(12,7))

#1st collimator velocity
plt.subplot(1,2,1)
plt.hist(ac.linac.v1, 50, facecolor='g', alpha=0.85)
plt.xlabel(r'$\beta$', fontsize=14)
plt.ylabel('Broj čestica', fontsize=14)
plt.title('Histogram brzina', fontsize=16)
plt.xlim(v1-4*v1std, v1+4*v1std)
#plt.ylim(0, 0.03)
plt.grid(True)


#1st collimator energy
plt.subplot(1,2,2)
plt.hist(ac.linac.E1, 50, facecolor='g', alpha=0.85)
plt.xlabel("E[MeV]", fontsize=14)
plt.ylabel('Broj čestica', fontsize=14)
plt.title('Histogram energija', fontsize=16)
plt.xlim(E1-4*E1std, E1+4*E1std)
#plt.ylim(0, 0.03)
plt.grid(True)

plt.tight_layout()
plt.suptitle("1. kolimator", fontsize=20)
plt.subplots_adjust(top=0.88)
plt.savefig(r'/Users/roccobarac/Documents/Bachelors_Thesis_files/linac4-1_cromod.png')
plt.close()




#2nd collimator plot
plt.figure(figsize=(12,7))

#2nd collimator velocity
plt.subplot(1,2,1)
plt.hist(ac.linac.v2, 50, facecolor='g', alpha=0.85)
plt.xlabel(r'$\beta$', fontsize=14)
plt.ylabel('Broj čestica', fontsize=14)
plt.title('Histogram brzina', fontsize=16)
plt.xlim(v2-4*v2std, v2+4*v2std)
#plt.ylim(0, 0.03)
plt.grid(True)


#2nd collimator energy
plt.subplot(1,2,2)
plt.hist(ac.linac.E2, 50, facecolor='g', alpha=0.85)
plt.xlabel("E[MeV]", fontsize=14)
plt.ylabel('Broj čestica', fontsize=14)
plt.title('Histogram energija', fontsize=16)
plt.xlim(E2-4*E2std, E2+4*E2std)
#plt.ylim(0, 0.03)
plt.grid(True)

plt.tight_layout()
plt.suptitle("2. kolimator", fontsize=20)
plt.subplots_adjust(top=0.88)
plt.savefig(r'/Users/roccobarac/Documents/Bachelors_Thesis_files/linac4-2_cromod.png')
plt.close()