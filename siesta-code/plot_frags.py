#!/usr/bin/env python3

import pickle
import numpy as np
import mass_spec as ms
import matplotlib.pyplot as plt

with open('fragments.pickle','rb') as f:
    fragments = pickle.load(f) #fragments consist of list with [geo][ion][frag]
    print(fragments)

frag_data = []
for geo in range(len(fragments)):
    for ion in range(len(fragments[geo])):
        for frag in range(len(fragments[geo][ion])):
            mol = []#Changed from str to array to fit into enumerate, must change after in this script
            for atom in range(len(fragments[geo][ion][frag])):
                element = ''.join([i for i in fragments[geo][ion][frag][atom] if not i.isdigit()])
                mol.append(element) #get atom to molecule
            print(mol)
            frag_mass = ms.get_mass(mol)
            frag_data.append(frag_mass)
arr_data = np.array(frag_data)
print(arr_data)
hist = np.histogram(arr_data,bins=int(np.max(arr_data))-1)
print(hist)
count = np.append(hist[0],0)
plt.plot(hist[1],count)

plt.show()


