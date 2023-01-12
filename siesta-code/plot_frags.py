#!/usr/bin/env python3

import pickle
import numpy as np
import mass_spec as ms
import matplotlib.pyplot as plt

job = 'MetINim-z4'
parent_mass = 254
nr_runs = 20

with open('fragments.pickle','rb') as f:
    fragments = pickle.load(f) #fragments consist of list with [geo][ion][frag]
    print(fragments)

frag_data = []
for geo in range(len(fragments)):
    for ion in range(len(fragments[geo])):
        for frag in range(len(fragments[geo][ion])):
            mol = [] #Changed from str to array to fit into enumerate, must change after in this script
            for atom in range(len(fragments[geo][ion][frag])):
                element = ''.join([i for i in fragments[geo][ion][frag][atom] if not i.isdigit()])
                mol.append(element) #get atom to molecule
            frag_mass = ms.get_mass(mol)

            frag_data.append(frag_mass)
arr_data = np.array(frag_data)
hist = np.histogram(arr_data,bins=int(np.max(arr_data))-1)
#To start masspec at 0
mass = np.insert(hist[1],0,0)
count = np.append(hist[0],0)
count = np.insert(count,0,0)
percent = count/nr_runs*100

#print(hist[0])
#print(count)
fig,ax = plt.subplots()

with open(job+'-fragment_data.txt','w') as f:
    f.write('Fragments obtained from Siesta molecular dynamics for system '+job+'\n')
    f.write('Mass (m/z)\tIntensity\tPercent\n')
    for i in range(mass.size):
        if count[i] > 0 :
            f.write(f'{mass[i]:>5}\t{count[i]:>10}\t{percent[i]:>10}\n')
ax.plot(mass,count)
ax.set_xlabel('m/z')
ax.set_ylabel('intensity (arb.u.)')
plt.savefig(job+'-mass_spec.png')
plt.show()


