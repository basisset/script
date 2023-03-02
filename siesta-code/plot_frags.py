#!/usr/bin/env python3

import pickle
import numpy as np
import mass_spec as ms
import matplotlib.pyplot as plt

job = 'Nim-z3'
parent_mass = 128  #metINim 254 Nim 128
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
            if frag_mass == 99:
                print(mol)
            frag_data.append(frag_mass)

arr_data = np.array(frag_data)
count = np.zeros(parent_mass+1)
for i in range(parent_mass+1):
    count[i]=frag_data.count(i)

#h,b,p = plt.hist(arr_data, bins=255)
#bins = np.linspace(0, parent_mass+1, 1000)
#count, bins = np.histogram(arr_data, bins=bins)
#print(bins[0], bins[1])
#mass = 0.5*(bins[1:] + bins[:-1])
#To start masspec at 0
percent = count/nr_runs*100

#fig,ax = plt.subplots()
with open(job+'-fragment_data.txt','w') as f:
    f.write('Fragments obtained from Siesta molecular dynamics for system '+job+'\n')
    f.write('Mass (m/z)\tIntensity\tPercent\n')
    for i in range(count.size):
        if count[i] > 0 :
            f.write(f'{i:>5}\t{count[i]:>10}\t{percent[i]:>10}\n')
    print('File saved!') 
#ax.plot(mass,count)
#ax.set_xlabel('m/z')
#ax.set_ylabel('intensity (arb.u.)')
#plt.savefig(job+'-mass_spec.png')
#plt.show()


