#!/usr/bin/env python   
import os, sys
import numpy as np
import shutil
import matplotlib.pyplot as plt
from statistics import mean, stdev
from numpy import linalg as LA
from scipy import interpolate
from itertools import combinations

# To analyze preparsed with "the naming convention", run e.g. ./analyze_preparsed.py ALA 1 10 1 10 C4 "H10 H11 H12" "Alanine C-H methyl bonds"
class atom:
    def __init__(self):
        self.name=""
        self.rvec=np.zeros(3)
        self.dvec = np.zeros(3) # direct
        self.pdos = np.zeros(1)
        self.sumdos = np.zeros(1)
        self.color = 0 # color index is used to find color from list in plotter, otherwise it's messier to change.
        self.phonons = []
        self.speciesName = ""
        self.speciesNumber = 0
        self.specnum=0
        self.speciesZNumber = 0
        self.mass = 0.0
        self.hirshfeldcharge=0.0
        self.mulliken_legend=[]
        self.mulliken_charges=[]

    def distance(self,center=np.asarray([0.0,0.0,0.0])):
        return np.linalg.norm(np.subtract(self.rvec,center))
        # return float(np.sqrt(self.rvec[0]**2+self.rvec[1]**2+self.rvec[2]**2))

    def in_cluster(self,maxrad,center=np.asarray([0.0,0.0,0.0]),minrad=0.0):
        return (self.distance(center) <= float(maxrad) and self.distance(center) >= float(minrad))


class lattice:
     def __init__(self):
         self.bravais=np.zeros((3,3))
         self.reciprocal=np.zeros((3,3))
         self.atoms=[]
         self.lattparam=0.0
         self.indSpecies=[]   # Number of atoms for one individual specie
         self.numSpecies=0    # Number of species
         self.indSpeciesNames=[] 
         self.coordtype=""


def deleteContent(fName): #Clearing the previous input file
            with open(fName, "w"):
                pass

def Date_and_Time():
    from time import gmtime, strftime #Current date and time
    t = strftime("%Y-%m-%d %H:%M:%S", gmtime())
    return t

def parse_text_bond_data(filename):
    bond_integrity=[]
    f=open(filename,'r')
    for i, line in enumerate(f.readlines()):
        if line.split()[-1][-1] != "]":
            full_line = line.split()
        else:
            full_line = full_line + line.split()
            bond_integrity.append(np.asarray(filter(None,[element.strip('[]') for element in full_line[1:]])).astype(np.float))
    return np.asarray(bond_integrity)

def get_neighborlist(timestep,rmax):
    neighborlist=[]
    rmin = 0.1 # Do not include self
    for i, atm in enumerate(timestep):
        neighborlist.append(find_atoms_within_radius(timestep,atm.rvec,rmax,rmin))
    return neighborlist


def find_atoms_within_cartesian(cluster,xlim,ylim,zlim):
    indices=[]
    for i, atm in enumerate(cluster):
        within = ((float(xlim[0])<= float(atm.rvec[0]) <= float(xlim[1])) and
                  (float(ylim[0])<= float(atm.rvec[1]) <=float(ylim[1])) and
                  (float(zlim[0])<= float(atm.rvec[2]) <=float(zlim[1])))
        if within:
            indices.append(i)
    return indices

def find_atoms_within_radius(cluster,center,rmax,rmin=0.0): 
    indices=[]
    for i, atm in enumerate(cluster):
        if (atm.in_cluster(rmax,center,rmin)):
            indices.append(i)
    return indices

def get_neighborlist(timestep,rmax):
    neighborlist=[]
    rmin = 0.1 # Do not include self
    for i, atm in enumerate(timestep):
        neighborlist.append(find_atoms_within_radius(timestep,atm.rvec,rmax,rmin))
    return neighborlist



#checking if element is int
def Is_Int(s):
    try:
        int(s)
        return True
    except ValueError:
        return False

#Parsing .ANI file
def parse_ANI(filename):
    f = open(filename, 'r')
    contents = f.readlines()
    f.close()
    atoms=[]
    time_serie=[]
    for i in range(len(contents)):
        if (Is_Int(contents[i])):
            atoms_in_timestep=int(contents[i].split()[0])	   
            for j in range(i+2,i+2+int(atoms_in_timestep)):
                atoms.append(atom())
                atoms[-1].rvec=[float(contents[j].split()[k]) for k in range(1,4)]
                atoms[-1].name=contents[j].split()[0]
            time_serie.append(atoms)
            atoms=[]
    return atoms_in_timestep, time_serie


def distR(D):
    N = np.loadtxt(D, dtype=np.float, delimiter=',')    
    Q = [np.linalg.norm(a-b) for a, b in combinations(N, 2)]
    return Q


def dist_timestep(timestep,atom1,atom2):
    return np.linalg.norm(np.subtract(timestep[atom2].rvec,timestep[atom1].rvec))


def bond_broken(dist,mean,T=150):
    B=[]
    for num in range(0,T):
        try:
            a = np.sqrt((np.sum(dist[num]-mean))**2)-0.5
        except:
                a = a # will this work for simulations that broke before T=150?
        b = 0.03*a
        c = np.exp(b)
        d = 1+c
        e = 1/d 
        B.append(e)
    
    return np.asarray(B)

def bond_broken_2(dist, T, mean, sigma, lamda):
    B=[]
    for num in range(0,T):
        e = (1 + np.exp(lamda*(dist[num]-mean-sigma-0.5)))**(-1)
        B.append(e)
    
    return np.asarray(B)

def parse_hirsh(filename):
    f = open(filename, 'r')
    contents = f.readlines()
    f.close()
    timesteps=[]
    charges=[]
    numatm=0
    for i in range(len(contents)):
        if ("NumberOfAtoms" in contents[i]):
            numatm=contents[i].split()[1]
        if ("Atom #    Qatom  Species" in contents[i]):
            for j in range(i+1,i+1+int(numatm)):
                charges.append(contents[j].split()[1])
            timesteps.append(np.asarray(charges,dtype=float))
            charges=[]
    return timesteps


def mean_distance_dict(thermalization_list, index_to_atom, neighbors_list):
    """Returns a dictionary with keys in the form '(Atom_A_index, Atom_B_index)' with values in the form of lists, where the elements are mean values for the distance between atom A and B over time for each thermalization run in thermalization list (which contains information from parse_timestep function). index_to_atom is a dictionary made from make_atom_dictionary_from_timeserie() and neighbors_list from get_neighborlist()"""
    mean_distance_dict = {} # Dictionary with every atom pair_kj, in tuple-form, as keys and their values are mean distances over
                        # all time positions for every thermalization run

    for run_index in range(len(thermalization_list)):
        distance_list = []
        for k in range(len(neighbors_list)): # atom_k, atom_j = atom pair_kj
            distance_lexi = {}
            for j in neighbors_list[k]:      # distance_list has dicts for every atom_k with neighbor atom_j as key and with 
                                        # distance_atom pair_kj(t) as values
                distance_lexi[str(j)] = [dist_timestep(thermalization_list[run_index][t_i],k,j) for t_i in
                                         range(len(thermalization_list[run_index]))]# Replaces previous list !!!!!
            distance_list.append(distance_lexi)    
            for j in neighbors_list[k]:
                KJ = str((index_to_atom[str(k)],index_to_atom[str(j)]))
                 
                if mean_distance_dict.get(KJ) is not None:
                    mean_distance_dict[KJ].extend([mean(distance_list[k][str(j)])])#adds mean distance values
                else:
                    mean_distance_dict[KJ] = [mean(distance_list[k][str(j)])]
        if run_index == 0:
             print(distance_list[2])
    return mean_distance_dict, distance_list


def frags_from_dists(mean_distances_dict, atom_to_index, ion_dict, lamda, cutoff_BI):
    """Returns a list of lists for which fragments there are at the last timestep, given a lambda value and a cutoff value for the bond integrity. The dimensions of the list are: total_fragments[geo][ion][fragment][atom]."""
    l=lamda #Lambda value
    #n_geo = 11
    #n_ion = 10  #Number of ionizations?
    #n_geo = 1
    #n_ion = 1
    
    broken_bonds_dict={}
    for bond in list(mean_distances_dict.keys()):
        n_geo = len(ion_dict[bond])
        broken_bonds_dict[bond]=[None]*n_geo
        for geo in range(n_geo):
            n_ion = len(ion_dict[bond][geo])
            broken_bonds_dict[bond][geo]=[None]*n_ion
            for ion in range(n_ion):
                BI = bond_broken_2(ion_dict[bond][geo][ion],len(ion_dict[bond][geo][ion]),
                                mean(mean_distances_dict[bond]), stdev(mean_distances_dict[bond]),l)
                if BI[-1] <= cutoff_BI:
                    broken_bonds_dict[bond][geo][ion]="broken"
                else:
                    broken_bonds_dict[bond][geo][ion]="intact"

    total_fragments=[None]*n_geo
    total_fragments=[[None]*n_ion for x in total_fragments]

    print(f'####################{n_geo}####################')
    for geo in range(n_geo):
        for ion in range(n_ion):
            polyatomic=[]
            monoatomic=[]
            for bond in broken_bonds_dict.keys():
                atoms= [x for x in atom_to_index.keys() if bond.split("'")[1]==x or bond.split("'")[3]==x]
                if broken_bonds_dict[bond][geo][ion]=="intact":
                    found=False
                    merged=False
                    for j in range (len(polyatomic)):
                        if (atoms[0] in polyatomic[j]) and (atoms[1] not in polyatomic[j]):
                            for k in range(len(polyatomic)):
                                if atoms[1] in polyatomic[k] and atoms[0] not in polyatomic[k]:
                                    polyatomic[j].extend(polyatomic[k])
                                    polyatomic[k]=[]
                                    merged=True
                                    break
                            if not merged:
                                polyatomic[j].append(atoms[1])
                            found=True
                        elif (atoms[1] in polyatomic[j]) and (atoms[0] not in polyatomic[j]):
                            for k in range(len(polyatomic)):
                                if atoms[0] in polyatomic[k] and atoms[1] not in polyatomic[k]:
                                    polyatomic[j].extend(polyatomic[k])
                                    polyatomic[k]=[]
                                    merged=True
                                    break
                            if not merged:
                                polyatomic[j].append(atoms[0])
                            found=True
                        elif (atoms[0] in polyatomic[j]) and (atoms[1] in polyatomic[j]):
                            found=True
                        else:
                            pass
                    if found==False: #This is the case where the atoms do not occur anywhere in the current version of "polyatomic"
                        polyatomic.append(atoms)
                elif broken_bonds_dict[bond][geo][ion]=="broken":
                    atom0_in_somefrag=False
                    atom1_in_somefrag=False
                    for other_bond in broken_bonds_dict.keys():
                        if broken_bonds_dict[other_bond][geo][ion]=="intact":
                            if (atoms[0]==other_bond.split("'")[1] or atoms[0]==other_bond.split("'")[3]):
                                atom0_in_somefrag=True
                            if (atoms[1]==other_bond.split("'")[1] or atoms[1]==other_bond.split("'")[3]):
                                atom1_in_somefrag=True
                    if not atom0_in_somefrag and atoms[0] not in monoatomic:
                        monoatomic.append(atoms[0])
                    if not atom1_in_somefrag and atoms[1] not in monoatomic:
                        monoatomic.append(atoms[1])
                    else:
                        pass
                else:
                    print("broken_bonds_dict["+bond+"]["+geo+"]["+ion+"] was not assigned a value")
       
            for i in range (len(monoatomic)):
                monoatomic[i]=[monoatomic[i]]
            empty_indices=[]
            for j in range(len(polyatomic)):
                if not polyatomic[j]:
                    empty_indices.append(j)
            for empty_index in reversed(empty_indices):
                del polyatomic[empty_index]
            fragments=[]
            fragments.extend(polyatomic)
            fragments.extend(monoatomic)
            total_fragments[geo][ion]=fragments
    return total_fragments


def write_xyz_anim(filename,timesteps,skipstep=1):
    f = open(filename,'w')
    for i, step in enumerate(timesteps):
        if (np.mod(i,skipstep)<0.5):
            f.write(str(len(step))+"\n")
            f.write('Timestep: '+str(i*skipstep)+"\n")
            for atm in step:
                f.write(str(atm.name)+" "+str(atm.rvec[0])+" "+str(atm.rvec[1])+" "+str(atm.rvec[2])+"\n")


def parse_hirsh_from_file(ion,lastion,acid):
    all_mean_hirsh=[]
    all_std_hirsh=[]
    for ionstage in range(ion,lastion):
        hirsh=[]
        print("acid, ionstage: ", str(acid), str(ionstage))
        for geostage in range(geometry,lastgeometry):
            Sim = './startgeo{0}_ionization{1}'.format(geostage,ionstage)
            os.chdir(Sim)
            try:
                hirsh.append(parse_hirsh("./stdout"))
            except:
                print("Failed to parse Hirshfeld for: {0}/startgeo{1}_ionization{2}".format(acid, geostage, ionstage))
            os.chdir("..")
        #print np.asarray(hirsh).mean(0).shape
        mean_data_name='{0}_hirshfeld_charge_{1}_hirshrun.dat'.format(acid,ionstage)
        np.savetxt(mean_data_name,np.asarray(hirsh).mean(0))
        mean_data_name='{0}_stdev_hirshfeld_charge_{1}_hirshrun.dat'.format(acid,ionstage)
        np.savetxt(mean_data_name,np.asarray(hirsh).std(0))
        all_mean_hirsh.append(np.asarray(hirsh).mean(0))
        all_std_hirsh.append(np.asarray(hirsh).std(0))
    return all_mean_hirsh, all_std_hirsh

def parse_eigenvalues(filename):
    f = open(filename, 'r')
    contents = f.readlines()
    f.close()
    timeserie_eig=[]
    timeserie_occ=[]

    for i, line in enumerate(contents):
        if ("Timestep" in line):
            current_step=int(line.split()[1])
            print(current_step)
            num_eigens=int(line.split()[3])
            print(num_eigens)
            eigenvalues=[]
            occupations=[]
            for j in range(i+1,i+num_eigens+1):
                eigenvalues.append(np.asarray(contents[j].split()[0:2], dtype=float))
                occupations.append(np.asarray(contents[j].split()[3:5], dtype=float))
            # Transpose to get spin-channels as timeserie[itime][ispin][:]
            timeserie_eig.append(np.transpose(np.asarray(eigenvalues)))
            timeserie_occ.append(np.transpose(np.asarray(occupations)))

    return np.asarray(timeserie_eig), np.asarray(timeserie_occ)

def read_preparsed_hirsh(acid,ion):
    mean_data_name='{0}_hirshfeld_charge_{1}_hirshrun.dat'.format(acid,ion)
    mean_hirsh=np.loadtxt(mean_data_name, dtype=np.float)
    mean_data_name='{0}_stdev_hirshfeld_charge_{1}_hirshrun.dat'.format(acid,ion)
    std_hirsh=np.loadtxt(mean_data_name, dtype=np.float)
    return mean_hirsh, std_hirsh

def make_atom_dictionary(filename):
    natoms, md_verlet = parse_ANI(filename)
    atomdict={}
    name_list=[atm.name for atm in md_verlet[0]]
    for i, atm in enumerate(md_verlet[0]):
        new_atom_number=name_list[0:i].count(atm.name)+1
        key=str(i)
        value=atm.name+str(new_atom_number)
        atomdict[key]=value
    inverted_dict = dict(map(reversed, atomdict.items()))
    return atomdict, inverted_dict

def make_atom_dictionary_from_timeserie(timeserie):
    atomdict={}
    name_list=[atm.name for atm in timeserie[0]]
    #print(name_list), print(type(name_list[0]))
    for i, atm in enumerate(timeserie[0]):
        print(str(i), atm.name)
        new_atom_number=name_list[0:i].count(atm.name)+1
        key=str(i)
        value=atm.name+str(new_atom_number)
        atomdict[key]=value
    inverted_dict = dict(map(reversed, atomdict.items()))
    return atomdict, inverted_dict


def parse_xyz(filename):
    xyz=[]
    f = open(filename, 'r')
    contents = f.readlines()
    f.close()
    for line in contents:
        xyz.append(np.asarray(line.split()[1:4], dtype=float))
    return np.transpose(np.asarray(xyz))

def parse_timestep(filename, outfile=None):
    with open(filename, 'r') as f:
        contents = f.readlines()
        #här print("filename: "+str(filename))
        print("length of file: "+str(len(contents)))
    numatm=0
    basissize='SZP'
    time_pos=[]
    time_mulliken=[]
    timesteps=[]
    specieslegend={}
    numberlegend={}
    mulls=[]
    orblegend=[]
    for i in range(len(contents)):
        if ("NumberOfAtoms" in contents[i]):
            numatm=int(contents[i].split()[1])
            print("Number of Atoms: "+str(numatm))
            break

    for i in range(len(contents)):
        if ("PAO.BasisSize" in contents[i]):
            basissize=str(contents[i].split()[1])
            print("Basis Size: "+str(basissize))
            break


    for i in range(len(contents)):
        if ("SpinPolarized" in contents[i]):
            if ("true") in contents[i]:
               spins=2
            else:
               spins=1
            break
    print("Spin components: "+str(spins))
    
    for i in range(len(contents)):
        if ("AtomicSpecies" in contents[i]):
            #print "Found AtomicCoord..."
            for j in range(i+1,len(contents)):
                print(str(j-i)+"  "+str(contents[j].split()))
                numberlegend[str(j-i)]=str(contents[j].split()[3])
                if ("AtomicCoordinatesAndAtomicSpecies" in contents[j+1]):
                #    print "Found!"
                    break
            else:
                continue
            break

    for i in range(len(contents)):
        if ("ChemicalSpeciesLabel" in contents[i]):
            for j in range(i+1,len(contents)):
                specieslegend[str(contents[j].split()[0])]=str(contents[j].split()[2])
                if ("ChemicalSpeciesLabel" in contents[j+1]):
                    break
            else:
                 continue
            break
    # print numberlegend
    # print specieslegend

    
    for i in range(len(contents)):
        if ("(Ang)" in contents[i])  and ("outcoor" in contents[i]):
                atoms=[]
                for j in range(i+1,i+numatm+1):
                    atoms.append(atom())
                    atoms[-1].rvec=[float(contents[j].split()[k]) for k in range(0,3)]
                    atoms[-1].name=specieslegend[numberlegend[str(contents[j].split()[4])]]
                time_pos.append(atoms)
        elif ("(Bohr)" in contents[i]) and ("outcoor" in contents[i]):
                atoms=[]
                for j in range(i+1,i+numatm+1):
                    atoms.append(atom())
                    #print([contents[j].split()[k] for k in range(0,3)])
                    atoms[-1].rvec=np.multiply([float(contents[j].split()[k]) for k in range(0,3)], 0.529177249)  # Convert Bohr to Angstrom
                    atoms[-1].name=specieslegend[numberlegend[str(contents[j].split()[4])]]
                time_pos.append(atoms)

# Approximately two lines per atom times number of spins + overhead of a few lines
    approx_mulliken_size=numatm*(3+spins*2)

    for i in range(len(contents)):
        if ("mulliken: Atomic and Orbital Populations:" in contents[i]):
            mulls, orblegend= parse_mulliken(contents[i:i+approx_mulliken_size],numatm,basissize,spins,outfile)
            time_mulliken.append(mulls)
    #print([a.rvec for a in time_pos[0]])#time_pos is list of lists with coordinates of atom at each timestep
    return time_pos, time_mulliken, orblegend, specieslegend, numberlegend
        

