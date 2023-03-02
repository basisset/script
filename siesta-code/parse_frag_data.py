#!/usr/bin/env python3

import numpy as np
import scipy as sp
import sys
import pickle
from statistics import mean, stdev
from analyze_trajectories import *
import matplotlib.pyplot as plt
from elementdata import *

plt.rcParams["font.family"] = "times new roman"    # Changing font style and size for plots
plt.rc('font', size=12) 


#Reads data files for verlet runs
mol = 'Nim'
z = 'z3'
verlet=['MD.out','MD2.out']
nr_runs = 20


thermalization_list=[]
for run in verlet:
    time_pos, timeserie, orblegend, specieslegend, numberlegend = parse_timestep(run)
    thermalization_list.append(time_pos)# Creates list of time positions 
index_to_atom, atom_to_index=make_atom_dictionary_from_timeserie(time_pos)
print(f'{len(thermalization_list)} input files found! ')

neighbors_list = get_neighborlist(time_pos[0],1.8)# choose a specific time for which we can identify neighbours

#BrINim
#neighbors_list[1].extend([7])# The iodine-carbon bond is longer, and is therefore manually added
if mol == 'metINim':
    print(f'> Analyzing {mol}')
    neighbors_list[0].extend([7])
    #Iodine at place 0 and closest carbon is nr 7 in metINim
if mol == 'BrINim': 
    print(f'> Analyzing {mol}')
    neighbors_list[0].extend([8])
    #Iodine at place 0 and carbon at 8 in BrINim
else:
    print(f'> Analyzing {mol}')
print(f'> Neighbor list[atom][neighbor]: {neighbors_list}')                               
print(f'> Number of neighbor lists is {len(neighbors_list)}')

print('> Calculating mean and standard deviation for atomic distances...')
mean_distances_dict, distance_list = mean_distance_dict(thermalization_list, index_to_atom, neighbors_list)
print(mean_distances_dict)

#Change to save file for mean and std from user input
#inp = input('> Do you want to print mean and standard deviation? (Y/n) \n')
#if inp in ['Y', 'y','']:
#    for k in range(len(neighbors_list)):
#        for j in neighbors_list[k]:
#            print(f"Mean distance between {index_to_atom[str(k)]} and {index_to_atom[str(j)]}: \t"
#              f"{mean(mean_distances_dict[str((index_to_atom[str(k)],index_to_atom[str(j)]))])}")
#            print(f'Standard deviation: \t\t\t{stdev(mean_distances_dict[str((index_to_atom[str(k)],index_to_atom[str(j)]))])}',end='\n\n')




all_keys = list(mean_distances_dict.keys())
sorted_keys = [sorted(e) for e in list(mean_distances_dict.keys())]
print(f'> Length of neighbor keys is {len(all_keys)}')

index=[]
for i, key in enumerate(sorted_keys):
    index.append([i for i, keyi in enumerate(sorted_keys) if key==keyi])

for indexpair in index:
    if len(indexpair) > 1:
        if all_keys[indexpair[1]] in mean_distances_dict:
            del mean_distances_dict[all_keys[indexpair[1]]]

#Removed incorrect bonds (H is only bound to one)
mean_distances_dict = rm_incorrect_bonds(index_to_atom,mean_distances_dict)

#runs2=[['output-13.out'],['output-20.out'],['output-1.out'],['output-2.out'],['output-3.out'],['output-4.out'],['output-5.out'],['output-6.out'],['output-7.out'],['output-8.out'],['output-9.out'],['output-10.out'],['output-11.out'],['output-12.out'],['output-14.out'],['output-15.out'],['output-16.out'],['output-17.out'],['output-18.out'],['output-19.out']]
#ionization_list = []
#for geo in runs2:
#    for run in geo:
#        time_pos, timeserie, orblegend, specieslegend, numberlegend = parse_timestep(run)
#        ionization_list.append(time_pos)


runs=['output-1.out','output-2.out','output-3.out','output-4.out','output-5.out','output-6.out','output-7.out','output-8.out','output-9.out','output-10.out','output-11.out','output-12.out','output-13.out','output-14.out','output-15.out','output-16.out','output-17.out','output-18.out','output-19.out','output-20.out']
#runs=['output-1.out','output-2.out','output-3.out']

ionization_list=[]
for run in runs:
    time_pos, timeserie, orblegend, specieslegend, numberlegend = parse_timestep(run)
    ionization_list.append(time_pos)# Creates list of time positions 


n_geo = len(runs)   # Number of different starting geometries
n_ion = 1    # Number of different ionization numbers (# electrons removed)
ion_dict = {}
for key in list(mean_distances_dict.keys()):
    index_i = int(atom_to_index[key.split("'")[1]])
    index_j = int(atom_to_index[key.split("'")[3]])
    for geo in range(n_geo):
        for ion in range(n_ion):
            current_run=(n_ion*geo)+ion
            ion_dist_list = [dist_timestep(ionization_list[geo][t],index_i,index_j) for t in range(len(ionization_list[geo]))]

            if key not in ion_dict:
                ion_dict[key]=[None]*n_geo
                ion_dict[key]=[[None]*n_ion for x in ion_dict[key]]
            else:
                pass
            ion_dict[key][geo][ion]=ion_dist_list

#Bond integrity over time
#atom_pairs = ["('N3', 'C2')", "('I1', 'C3')"]#Key C2,N3 doesnt work
#atom_pairs = ["('I1', 'C2')"]
print(ion_dict.keys())
atom_pairs = ["('N3', 'C2')","('C2', 'C3')","('N2', 'C2')","('N2', 'C1')","('N1', 'C1')","('N1', 'C3')","('C1', 'C4')"]
nr_broken = 0

ion_run = 1 #number of ionizations
ion = ion_run - 1
i = None
j = None
time = [t for t in range(len(thermalization_list[0]))]
bond_intr_mean = np.zeros(len(time))
bond_intr_stat = atom_pairs.copy() #Bond breakage statistics
for bond,atom_pair in enumerate(atom_pairs):
    fig, ax = plt.subplots()
    for geo in range(n_geo):
        bond_intr = bond_broken_2(ion_dict[atom_pair][geo][ion],len(ion_dict[atom_pair][geo][ion]),mean(mean_distances_dict[atom_pair]), stdev(mean_distances_dict[atom_pair]),10)
        bond_intr_mean += bond_intr
        nr_broken += (bond_intr[-1]<0.5).sum()
        ax.plot(time, bond_intr,zorder=1000, color='gray')
        i = atom_pair.split("'")[1]
        j = atom_pair.split("'")[3]
        ax.set(xlabel='Time [fs]',ylabel='Bond integrity',title=f'Bond integrity for atom pair {i} {j} for {mol} {z}') 
    bond_intr_mean = bond_intr_mean/n_geo 
    bond_intr_stat[bond] += f': {(nr_broken/nr_runs)*100}%'
    nr_broken = 0
    ax.fill_between(time,bond_intr_mean,1,zorder=1, color='orange')
    plt.savefig(f'{i}_{j}-bondIntr-{mol}-{z}.png')
    print(bond_intr_stat)
    plt.show()
    bond_intr_mean = np.zeros(len(time))



total_fragments = frags_from_dists(mean_distances_dict, atom_to_index, ion_dict, lamda=10, cutoff_BI=0.5)
print(f'> Inputfiles:\n {runs}\nFragments: \n{total_fragments}')
print(f'> Bond integrity in final time step per bond is:\n {bond_intr_stat}')
with open('fragments.pickle','wb') as f:
    pickle.dump(total_fragments, f)
print('> Done')

#Write code which also shows mass of each fragment






