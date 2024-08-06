#!/usr/bin/env Python3
import numpy as np
import pandas as pd
import csv

nr_sims=50
#Updating names to fit a xyz-file
def rename(atom):
    if atom in ['OW', 'O1', 'O2']:
        return 'O'
    if atom in ['HW1', 'HW2', 'H1', 'H2']:
        return 'H'
    if atom=='I1':
        return 'I'
    if atom in ['N1', 'N2', 'N3']:
        return 'N'
    if atom in ['C1', 'C2', 'C3']:
        return 'C'
    if atom=='CL':
        return 'Cl'
    if atom=='BR1':
        return 'Br'

    raise RuntimeError(f'unknown atom: {atom}')
#To determine distances between molecules (sometimes the dists are too close between atoms after shifting)
def atomic_distance(moldata, totdata, threshold=1.5):
    totdata_arr = totdata.to_numpy()
    moldata_arr = moldata.to_numpy()
    to_rm = np.arange(len(totdata_arr))
    rm_count = 0
    closest_atom = np.zeros(3)
    closest_h2o = np.zeros(3)
    atoms_in_mol = len(moldata_arr)
    for i,molatom in enumerate(totdata_arr[:atoms_in_mol][:]):
        for j,solatom in enumerate(totdata_arr[atoms_in_mol:][:]): 
            #if j+atoms_in_mol in to_rm:
            #    continue
            distance = np.sqrt((molatom[5]-solatom[5])**2 + (molatom[6]-solatom[6])**2 + (molatom[7]-solatom[7])**2)
            #finding closest atom
            if (j==0 and i==0) and (distance > threshold):
                closest_atom = [j+atoms_in_mol, distance, totdata_arr[j+atoms_in_mol][2]]
            if (closest_atom[1]>distance) and (distance > threshold):
                closest_atom = [j+atoms_in_mol, distance, totdata_arr[j+atoms_in_mol][2]]
            # Distance condition
            #if (distance < threshold):
            #    if j+atoms_in_mol not in to_rm:
            #        rm_count += 1
            #    # Storing indeces of atoms to remove
            #    to_rm[rm_count] = j+atoms_in_mol
            #    print(f'Distance is: {distance}')
    #to_rm = to_rm[to_rm != 0]
    if closest_atom[2]=='OW':
        closest_h2o[0] = closest_atom[0] #index of closest atom
        closest_h2o[1] = closest_atom[0]+1
        closest_h2o[2] = closest_atom[0]+2
    if closest_atom[2]=='HW1':
        closest_h2o[1] = closest_atom[0] #index of closest atom
        closest_h2o[0] = closest_atom[0]-1
        closest_h2o[2] = closest_atom[0]+1
    if closest_atom[2]=='HW2':
        closest_h2o[2] = closest_atom[0] #index of closest atom
        closest_h2o[1] = closest_atom[0]-1
        closest_h2o[0] = closest_atom[0]-2
    
    to_rm[:atoms_in_mol] = 0
    for i,index in enumerate(to_rm):
        if index in closest_h2o:
            to_rm[i] = 0
    
    to_rm = to_rm[to_rm != 0]
    
    #Checks so that an atom in the closest molecule doesnt get removed. 
    #for atom in closest_h2o:
    #    if atom in to_rm:
    #        to_rm = np.delete(to_rm, np.where(to_rm==atom))
    print(f'> Atoms in original file: {len(totdata_arr)}')
    print(f'> Atoms to remove: {len(to_rm)}')
    print(f'> In positions: {to_rm}')
    return to_rm.astype(int), totdata_arr, closest_h2o

def delete_atoms(array_of_atoms, index_to_rm):
    data_del = np.delete(array_of_atoms, index_to_rm, axis=0)
    return data_del


#Reading and aranging data
for run in np.arange(nr_sims):
    print(f'> Reading file {run}...')
    fname = f'2-bromo-5-iodo-4-nitroimidazole_{run}.pdb'
    preamble = pd.read_csv(fname, nrows=5)
    end = pd.read_csv(fname, skiprows=29283)

    dtype = np.dtype([('ATOM', '|S4'),('index', np.int64),('atom', '|S3'),('molecule', '|S3'),('1', np.int64), ('x', np.float64), ('y', np.float64), ('z', np.float64), ('1.0', np.float64), ('0.0', np.float64), ('element', '|S2')])

    soldata = pd.read_csv(fname, skiprows=17, nrows=29266, sep='\s+', names=['ATOM','index','atom','molecule','1','x','y','z','1.0','0.0','element'])
    moldata = pd.read_csv(fname, skiprows=5, nrows=12, sep='\s+', names=['ATOM','index','atom','molecule','1','x','y','z','1.0','0.0','element'])
    new_moldata = moldata.copy()
    new_soldata = soldata.copy()

#Shifting data
    box_size = [30, 30, 330]
    center_index = 0#index of the atom which we want to put in the middle
    center_coord = moldata['z'].iloc[center_index]
    shift = (box_size[2]/2)-center_coord
    new_moldata['z'] = new_moldata['z']+shift

#if coordinate of atom is outside box then it is moved to the other side
    new_moldata.loc[new_moldata['z'] > box_size[2],'z'] = new_moldata.loc[new_moldata['z'] > box_size[2],'z'] - box_size[2]
    new_soldata.loc[new_soldata['z'] > box_size[2],'z'] = new_soldata.loc[new_soldata['z'] > box_size[2],'z'] - box_size[2]

    tot_data = pd.concat([new_moldata,new_soldata])
#tot_data.to_csv('temp.csv', columns=['atom', 'x', 'y', 'z'],index=False, sep='\t')

    to_rm, totdata_arr, closest_h2o = atomic_distance(new_moldata, tot_data)
    data_del = delete_atoms(totdata_arr, to_rm)
    updated_atoms = len(data_del)

    new_totdata = pd.DataFrame(data_del, columns=['ATOM','index','atom','molecule','1','x','y','z','1.0','0.0','element'])
#tot_data['atom'] = tot_data['atom'].apply(rename)
    new_totdata['atom'] = new_totdata['atom'].apply(rename)


    new_totdata.to_csv(f'sim_{run}.xyz', columns=['atom', 'x', 'y', 'z'],index=False, sep='\t')
    with open(f'sim_{run}.xyz', 'r+') as file:
        file_data = file.read()
        file.seek(0,0)
        file.write(f'{updated_atoms}' + '\n' + file_data)

