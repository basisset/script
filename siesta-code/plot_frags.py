#!/usr/bin/env python3

import pickle
import numpy as np
import mass_spec as ms

with open('fragments.pickle','rb') as f:
    fragments = pickle.load(f) #fragments consist of list with [geo][ion][frag]

frag_data = []
for geo in range(len(fragments)):
    for ion in range(len(fragments[geo])):
        for frag in range(len(fragments[geo][ion])):
            mol = ''
            for atom in range(len(fragments[geo][ion][frag])):
                mol += str(fragments[geo][ion][frag][atom]).split()[0][0]
            frag_mass = ms.get_mass(mol)
            frag_data.append(f'{mol} {frag_mass}')
            print(frag_data)
            


