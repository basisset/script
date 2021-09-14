#!/usr/bin/env python3   
import os, sys
import numpy as np
import shutil
# import matplotlib.pyplot as plt
from numpy import linalg as LA
from itertools import combinations

assert sys.version_info >= (3,0) # checks to be atleast Python3
dummy = 0
class atom:
    def __init__(self):
        self.rvec=[]
        self.name=""
        self.specnum=0
        self.neigbors=[]

    def distance(self,center=np.asarray([0.0,0.0,0.0])):#distance from point in space
        return np.linalg.norm(np.subtract(self.rvec,center))
        # return float(np.sqrt(self.rvec[0]**2+self.rvec[1]**2+self.rvec[2]**2))

    def in_cluster(self,maxrad,center=np.asarray([0.0,0.0,0.0]),minrad=0.0):
        return (self.distance(center) <= float(maxrad) and self.distance(center) >= float(minrad))

    def __repr__(self): # prints the attributes of atom in a nicer way
        return f"Atom: {self.name}, rvec: {self.rvec}, specnum: {self.specnum}, neighbor: {self.neigbors}"

def deleteContent(fName): #Clearing the previous input file
            with open(fName, "w"):
                pass
def Date_and_Time():
    from time import gmtime, strftime #Current date and time
    t = strftime("%Y-%m-%d %H:%M:%S", gmtime())
    return t


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

def Is_Float(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

#Parsing .ANI file
def parse_ANI(filename):
    global dummy
    with open(filename, 'r') as f:
        contents = f.readlines()
    atoms=[]
    time_serie=[]
    for i in range(len(contents)):
        if (Is_Int(contents[i])):
            atoms_in_timestep=int(contents[i].split()[0])	   
            for j in range(i+2,i+2+int(atoms_in_timestep)):
                atoms.append(atom())
                atoms[-1].rvec=[float(contents[j].split()[k]) for k in range(1,4)] # coords of atom in ANI file
                atoms[-1].name=contents[j].split()[0] # name of atom in ANI file
            time_serie.append(atoms)
            atoms=[]
    return atoms_in_timestep, time_serie


def distR(D):
    N = np.loadtxt(D, dtype=np.float, delimiter=',')     		
    Q = [np.linalg.norm(a-b) for a, b in combinations(N, 2)]
    return Q	 


def dist_timestep(timestep,atom1,atom2):
    return np.linalg.norm(np.subtract(timestep[atom2].rvec,timestep[atom1].rvec))


def bond_broken(dist,T=150):
    B=[]	
    for num in range(0,T+1): # + 1 so that we get the bond-integrity at 75 fs also
		
        try:
            a = np.sqrt((np.sum(dist[num]-dist[0]))**2)-0.5
        except:
            a = a # will this work for simulations that broke before T=150?
        b = 10*a
        c = np.exp(b)
        d = 1+c
        e = 1/d	
        B.append(e)
    
    return np.asarray(B)
    
def write_xyz_anim(filename,timesteps,skipstep=1):
    f = open(filename,'w')
    for i, step in enumerate(timesteps):
        if (np.mod(i,skipstep)<0.5):
            f.write(str(len(step))+"\n")
            f.write('Timestep: '+str(i*skipstep)+"\n")
            for atm in step:
                f.write(str(atm.name)+" "+str(atm.rvec[0])+" "+str(atm.rvec[1])+" "+str(atm.rvec[2])+"\n")

#if len(sys.argv)==1:
#    print(str(sys.stderr)+ "Arg 1: amino acid, Arg 2: start-index of geometry file, Arg 3: end-index of geometry file, Arg 4: ionization stage "+str(sys.argv[0]))
#    exit(1);


#Script START
def write_new_ANIfile(stdout):
    os.chdir(stdout)
    f = open('output', 'r')
    lines = f.readlines()
    t = False
    printing = False ###What is first_timestep.ANI and why do you need it?
    fle = open('first_timestep.ANI', 'w')

    fle.write('   ' + '33\n\n' )

    for line in lines:

        if '%block AtomicCoordinatesAndAtomicSpecies' in line and 'endblock' not in line:  
            if printing == True:
                t = True
        if t == True:
            if 'block' not in line:
               
                if Is_Float(line.split()[0]):
                    if '1' == line.split()[3]:
                        atomicID = 'O'
                    if '2' == line.split()[3]:
                        atomicID = 'H'
                    if '3' == line.split()[3]:
                        atomicID = 'N'
                    if '4' == line.split()[3]:
                        atomicID = 'C'
                    if '5' == line.split()[3]:
                        atomicID = 'S'
                    coordinates = line.split()
                    fle.write(atomicID + '       ' + coordinates[0] + '  ' +  coordinates[1] + '  ' + coordinates[2] + '\n')

        if '%endblock ChemicalSpeciesLabel' in line:

            printing = True
        
        if '%endblock AtomicCoordinatesAndAtomicSpecies' in line:
            break
        
    fle = open('first_timestep.ANI','r')
    f2 = open('list2.ANI','w')
    for x in fle.readlines():
        f2.write(x)
    fle.close()
    f2.close()

    for element in os.listdir(stdout):
        if element.endswith('.ANI') and 'alanine' in element:
            print(element)
            f=open(element)
            break
    #f=open('*.ANI')  
    f1=open('list2.ANI','a')
    for x in f.readlines():
        f1.write(x)
    f.close()
    f1.close()

#natoms, md_verlet = parse_ANI('/home/ibrahim/phd/PEPTIDE-PEK/peptides/peptides/trialanine/alpha_helix_ionization/1/DFT_ionization_simulation/startgeometry1_ionization0/list2.ANI') # read any .ANI file you have to get the structure


natoms, md_verlet = parse_ANI('/home/pamela/projects/methylcysteine/SIESTA/startstruct/metcyst.ANI') # read any .ANI file you have to get the structure

neiglist =  get_neighborlist(md_verlet[0],2.0)  # Checks which atoms that are in range of 2 angstrom from a specific atom for time zero to make a list. Change the value if you want a larger range. 


geometry = 1 
lastgeometry = 3 # here we are analyzing 10 different trajectories (11 here for lastgeometry since the "range" function in python does not take the last value). 
 
#os.chdir('/home/ibrahim/phd/PEPTIDE-PEK/peptides/peptides/trialanine/alpha_helix_ionization/') # os.chdir changes where we "are" to where your siesta data is
os.chdir('/home/pamela/projects/methylcysteine/SIESTA/charge/') # os.chdir changes where we "are" to where your siesta data is

output_path = os.getcwd() + '/output' #Get current working path + new output folder

try:
    # create an output folder where the data will be stored
    os.mkdir(output_path) 
except:
    # remove if there is already a folder, and create a new one.
    shutil.rmtree(output_path)
    os.mkdir(output_path)

acid = 'metcyst' # the name of the system you are studying

for main, atmlst in enumerate(neiglist): 
    for neig in atmlst:
        print("Analyzing "+ str(acid) + " atoms: "+ str(main)+ " to "+ str(neig))
        mean_integrity=[]
        std_integrity=[]
        counter = 0
        bond_integrity=[]
        for geostage in range(geometry,lastgeometry): 
            
                #Sim = '/home/ibrahim/phd/PEPTIDE-PEK/peptides/peptides/trialanine/alpha_helix_ionization/{0}/DFT_ionization_simulation/startgeometry{0}_ionization{1}'.format(geostage,ionstage)
            Sim = '/home/pamela/projects/methylcysteine/SIESTA/charge/SIM-{0}'.format(geostage)
            os.chdir(Sim) # Sim is there directory of the specific run of the geometry "geostage" and ionization "ionstage". os.chdir takes us there
            natoms, md_verlet = parse_ANI('./metcyst.ANI') # read the .ANI file in this trajectory. Might be called something else for you..
                #write_new_ANIfile(os.getcwd())
                #natoms, md_verlet = parse_ANI('./list2.ANI')
              
                # Trajectory parsed, needs to be done:
                # 1)   Calculate bond integrity between all adjacent atoms
                # 2)   Save to some data structure
                # 3)   Average between trajectories (after ionstage loop)
                # 4)   Make plots, check standard deviation

            dists=[dist_timestep(timestep,main,neig) for timestep in md_verlet] # check distances between the atoms for each timestep
            bond_integrity.append(np.asarray(bond_broken(dists))) # calculate bond integrity, 
            os.chdir("..") 
                # now repeat for the next geometry
            

        mean_integrity.append(np.transpose(np.asarray(bond_integrity)).mean(1)) # calculate the mean over the 10 (or more/less depending on what you have) trajectories
        std_integrity.append(np.transpose(np.asarray(bond_integrity)).std(1)) # calculate the standard deviation over the 10 trajecories
        
        # SAVE DATA

        mean_data_name='{0}_mean_bond_integrity_{1}{2}-{3}{4}_num_{5}_hirshrun.txt'.format(acid,md_verlet[0][main].name,main,md_verlet[0][neig].name,neig, counter)
        np.savetxt(output_path + '/' + mean_data_name,np.transpose(mean_integrity))        
        std_integrity=np.asarray(std_integrity)
        std_data_name='{0}_stdev_bond_integrity_{1}{2}-{3}{4}_num_{5}_hirshrun.txt'.format(acid,md_verlet[0][main].name,main,md_verlet[0][neig].name,neig, counter)

        np.savetxt(output_path + '/' + std_data_name,np.transpose(std_integrity))
        mean_integrity=[]
        std_integrity=[]

