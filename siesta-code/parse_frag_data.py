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
runs = ['output-10.out','output-7.out']
#runs=['output-1.out','output-3.out','output-4.out','output-5.out','output-6.out','output-7.out','output-8.out','output-9.out','output-10.out']
#runs=['output-1.out','output-2.out','output-3.out','output-4.out','output-5.out','output-6.out','output-7.out','output-8.out','output-9.out','output-10.out','output-11.out','output-12.out','output-14.out','output-15.out','output-16.out','output-17.out','output-18.out','output-19.out']
#runs=['output-1.out','output-2.out','output-3.out','output-4.out','output-5.out','output-6.out','output-7.out','output-8.out','output-9.out','output-10.out','output-11.out','output-12.out','output-13.out','output-14.out','output-15.out','output-16.out','output-17.out','output-18.out','output-19.out','output-20.out']


thermalization_list=[]
for run in runs:
    time_pos, timeserie, orblegend, specieslegend, numberlegend = parse_timestep(run)
    thermalization_list.append(time_pos)# Creates list of time positions 
index_to_atom, atom_to_index=make_atom_dictionary_from_timeserie(time_pos)
print(f'{len(thermalization_list)} input files found! ')

neighbors_list = get_neighborlist(time_pos[0],1.8)# choose a specific time for which we can identify neighbours
#BrINim
#neighbors_list[0].extend([8])
#neighbors_list[1].extend([7])# The iodine-carbon bond is longer, and is therefore manually added
#neighbors_list[0].extend([7])#Iodine at place 0 and closest carbon is nr 7 in metINim
#Iodine at place 0 and carbon at 8 in BrINim
print(f'> Neighbor list[atom][neighbor]: {neighbors_list}')                               
print(f'> Number of neighbor lists is {len(neighbors_list)}')

print('> Calculating mean and standard deviation for atomic distances...')
mean_distances_dict, distance_list = mean_distance_dict(thermalization_list, index_to_atom, neighbors_list)


#Change to save file for mean and std from user input
inp = input('> Do you want to print mean and standard deviation? (Y/n) \n')
if inp in ['Y', 'y','']:
    for k in range(len(neighbors_list)):
        for j in neighbors_list[k]:
            print(f"Mean distance between {index_to_atom[str(k)]} and {index_to_atom[str(j)]}: \t"
              f"{mean(mean_distances_dict[str((index_to_atom[str(k)],index_to_atom[str(j)]))])}")
            print(f'Standard deviation: \t\t\t{stdev(mean_distances_dict[str((index_to_atom[str(k)],index_to_atom[str(j)]))])}',end='\n\n')



# Analyzing bond distance
#atom = 4# Only one neigbor-list is chosen
#time = [x for x in range(len(thermalization_list[0]))]
#for granne in neighbors_list[atom]:# For each atom in the list
    #fig, ax = plt.subplots()
#    print(f'####Bonds to atom {index_to_atom[str(granne)]} analyzed')
  #  print(f'Lengths of files: {len(time)} and {len(distance_list[i][str(granne)])}')#in case output list is shorter than time (crashed calculation)
    #ax.plot(time, distance_list[atom][str(granne)])
    #ax.set(xlabel='Time [fs]', ylabel='Distance [AA]',title=f'Distance between atom pair {index_to_atom[str(atom)]} {index_to_atom[str(granne)]}')
    #fig.savefig(f'dist{index_to_atom[str(i)]} {index_to_atom[str(j)]}')
#    plt.show()

#Analyzing bond integrity for one pair
#atom_pair = "('N3', 'C2')"
#atom_i = atom_pair.split("'")[1]
#i = atom_to_index[atom_i]
#atom_j = atom_pair.split("'")[3]
#j = atom_to_index[atom_j]
#lamda = 10

#bond_dists = distance_list[int(i)][str(j)]
#fig, ax = plt.subplots()
#ax.plot(bond_dists, bond_broken_2(bond_dists,len(bond_dists), mean(mean_distances_dict[atom_pair]),
#                                  stdev(mean_distances_dict[atom_pair]), lamda), c='orange')
#ax.set(xlabel='Bond distance [AAngstroem]', ylabel='Bond integrity', title=f'Bond integrity for atom pair {atom_i} {atom_j}')
#plt.show()

#Comparing bond intergity for atom pair with general form of bond integrity using bond_broken_2
#d = np.arange(min(bond_dists),3,0.01)
#fig, ax = plt.subplots()
#I = int(i)
#J = str(j)
#for lamda in [1,10,100]:
#    ax.plot(d,bond_broken_2(d,len(d),mean(mean_distances_dict[atom_pair]),stdev(mean_distances_dict[atom_pair]), lamda),label=f'BI function l={lamda}')
#    ax.plot(distance_list[I][J], bond_broken_2(distance_list[I][J],len(distance_list[I][J]), mean(mean_distances_dict[atom_pair]), stdev(mean_distances_dict[atom_pair]), lamda),label=f'Thermalization run l={lamda}')
#
#ax.set_ylim([0,1.1])
#ax.legend()
##ax.set(xlabel='Bond distance [Å]', ylabel='Bond integrity', title=f'Bond integrity for {atom_i} {atom_j} with lambda={lamda}')
#ax.set(xlabel='Bond distance [AA]', ylabel='Bond integrity')
##fig.savefig(f'lamda{lamda} {atom_i} {atom_j}')
#plt.show()


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





n_geo = len(runs)   # Number of different starting geometries # går programmet in o kollar alla körningar?
n_ion = 1    # Number of different ionization numbers (# electrons removed)
ion_dict = {}
for key in list(mean_distances_dict.keys()):
    index_i = int(atom_to_index[key.split("'")[1]])
    index_j = int(atom_to_index[key.split("'")[3]])
    for geo in range(n_geo):
        for ion in range(n_ion):
            current_run=(n_ion*geo)+ion
            ion_dist_list = [dist_timestep(thermalization_list[geo][t],index_i,index_j) for t in range(len(thermalization_list[geo]))]

            if key not in ion_dict:
                ion_dict[key]=[None]*n_geo
                ion_dict[key]=[[None]*n_ion for x in ion_dict[key]]
            else:
                pass
            ion_dict[key][geo][ion]=ion_dist_list

#Bond integrity over time
#atom_pairs = ["('N3', 'C2')", "('I1', 'C3')"]#Key C2,N3 doesnt work
#atom_pairs = ["('I1', 'C2')"]#Key C2,N3 doesnt work
atom_pairs = ["('N3', 'C2')"]

ion_run = 1 #number of ionizations
ion = ion_run - 1
i = None
j = None
time = [t for t in range(len(thermalization_list[0]))]
bond_intr_mean = np.zeros(len(time))
for atom_pair in atom_pairs:
    fig, ax = plt.subplots()
    for geo in range(n_geo):
        bond_intr = bond_broken_2(ion_dict[atom_pair][geo][ion],len(ion_dict[atom_pair][geo][ion]),mean(mean_distances_dict[atom_pair]), stdev(mean_distances_dict[atom_pair]),10)
        bond_intr_mean += bond_intr
        ax.plot(time, bond_intr,zorder=1000, color='gray')
        i = atom_pair.split("'")[1]
        j = atom_pair.split("'")[3]
        ax.set(xlabel='Time [fs]',ylabel='Bond integrity',title=f'Bond integrity for atom pair {i} {j} for geometry nr {min(range(geo+1))+1}-{max(range(geo+1))+1}') 
    bond_intr_mean = bond_intr_mean/n_geo 
    ax.fill_between(time,bond_intr_mean,1,zorder=1, color='orange')
#    plt.savefig(f'{i}_{j}-bondIntr.png')
    #plt.show()
    bond_intr_mean = np.zeros(len(time))


### To produce heat maps ###
#atom_pair = "('I1', 'C8')"
#i = atom_pair.split("'")[1]
#j = atom_pair.split("'")[3]
#
#mean_g_dist = [[] for _ in range(n_geo)]
#all_g = [[] for _ in range(n_geo)]
#for ion in range(n_ion):
#    all_g[ion] = [ion_dict[atom_pair][g][ion] for g in range(n_geo) if len(ion_dict[atom_pair][g][ion]) == len(time)]
#
#    for t in range(len(time)):
#        mean_g_dist[ion].append(mean([all_g[ion][g][t] for g in range(len(all_g[ion]))]))
#
#
#z_mesh = np.divide(np.linspace(0,n_ion,10),np.float(30))
#time_mesh = [t for t in range(len(ionization_list[0]))]
#all_mean_integrity = np.transpose([bond_broken_2(mean_g_dist[current_i], len(mean_g_dist[current_i]),
#                    mean(mean_distances_dict[atom_pair]), stdev(mean_distances_dict[atom_pair]),10) for current_i in range(10)])
#
#
#fig, ax = plt.subplots(figsize=(5,4))
##fig, ax = plt.subplots()
#p = plt.contourf(z_mesh, time_mesh, all_mean_integrity, levels=100, vmin=0., vmax=1.0,
#             alpha=1, cmap='plasma')
#
#x_labels = [round(x/30,2) for x in range(1,11)]   #Specifying the ticks on x-axis, to show simulated values exactly, 1/30-10/30
#plt.xticks([x/30 for x in range(1,11)], x_labels)
#
##fig.colorbar(p, ticks=[0,0.2,0.4,0.6,0.8,1], label='Mean bond integrity')
##plt.clim(0,1)
##plt.rc('font', size=12)
#ax.set(xlabel='$\overline{z}$ [e/N]', ylabel='Time [fs]')
#fig.savefig(f'BI {i} {j}')
##ax.set(xlabel='$\overline{z}$ [e/N]', ylabel='Time [fs]', title=f'Mean bond integrity for atom pair {i} {j}')
#plt.show()



total_fragments = frags_from_dists(mean_distances_dict, atom_to_index, ion_dict, lamda=10, cutoff_BI=0.5)
print(total_fragments)
with open('fragments.pickle','wb') as f:
    pickle.dump(total_fragments, f)
print('> Done')

#Write code which also shows mass of each fragment






