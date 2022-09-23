import numpy as np
import matplotlib.pyplot as plt


plt.figure(figsize=(12,7))

plt.subplot(1,2,1)
plt.plot([0,1], [1,2])
plt.xlabel(r'$\beta$', fontsize=14)
plt.ylabel('Number of particles', fontsize=14)
plt.title('Speed distribution histogram', fontsize=16)

plt.subplot(1,2,2)
plt.plot([0,1], [0,1])
plt.xlabel(r'$\beta$', fontsize=14)
plt.ylabel('Number of particles', fontsize=14)
plt.title('Speed distribution histogram', fontsize=16)
plt.tight_layout()
plt.suptitle("TEST", fontsize=20)
plt.subplots_adjust(top=0.88)
plt.savefig(r'/Users/roccobarac/Documents/Bachelors_Thesis_files/test.png')
plt.close()


plt.figure(figsize=(12,7))
plt.subplot(1,2,1)
plt.plot([0,1], [1,3])
plt.xlabel(r'$\beta$', fontsize=14)
plt.ylabel('Number of particles', fontsize=14)
plt.title('Speed distribution histogram', fontsize=16)

plt.subplot(1,2,2)
plt.plot([0,1], [0,120])
plt.xlabel(r'$\beta$', fontsize=14)
plt.ylabel('Number of particles', fontsize=14)
plt.title('Speed distribution histogram', fontsize=16)
plt.tight_layout()
plt.suptitle("TEST", fontsize=20)
plt.subplots_adjust(top=0.88)
plt.savefig(r'/Users/roccobarac/Documents/Bachelors_Thesis_files/test2.png')
plt.close()