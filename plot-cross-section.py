import numpy as np
from scipy.interpolate import interp1d
#import sys
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

#argv = sys.argv[1:]

def main(argv):
	

    flux_data = np.loadtxt('B1_Pt_400.txt')
    
    # PFOA + BEMP
    element = np.array(['P','F','C','O','N','Total']) #Label for the elements
    nr_atoms = [1,15,21,2,4]

    # Cytochrome C
    #element = np.array(['Fe','S','O','N','C','Total'])
    #nr_atoms = [1,2,6,8,42]
    
    colorkey = np.array(['yellow','gray','red','blue','black']) # Only important for cross-sections per element
	#tot = np.zeros(len(np.loadtxt(argv[0]))) 
    
    n = len(argv)
    tot_intensity = np.zeros(1248)
    plt.figure(figsize=[8.0,4.8])
	#plt.xlim([-10,10])
    
    interp_func = interp1d(flux_data[:,0],flux_data[:,1])
    interp_x = np.arange(250,1498,1)
    newarr_flux = interp_func(interp_x)
    norm_flux = newarr_flux/np.max(newarr_flux)
    #norm_flux = newarr_flux.copy()
    for i in range(n):
        data = np.loadtxt(argv[i])
        x = data[:,0]
        y = data[:,1]*nr_atoms[i]
        
        interp_func = interp1d(x,y)
        newarr = interp_func(interp_x)
        tot_intensity += newarr
        

        #plt.plot(interp_x,newarr,color=colorkey[i],label=element[i]) #to plot cross-section of each atom	
    #tot_intensity = tot_intensity/np.max(tot_intensity)
    
    plt.plot(interp_x,tot_intensity,color=colorkey[-1],label=element[-1])
    plt.plot(interp_x,norm_flux,label='Flux')
    #plt.plot(interp_x,norm_flux*tot_intensity,label='Max cross-section') # to get max cross section with respect to flux
    plt.grid()

    plt.xlabel('Photon Energy [eV]', fontsize='x-large')
    plt.ylabel('Intensity [Arb.u]', fontsize='x-large')
    plt.legend(loc=2,edgecolor='none',fontsize='large')
    plt.tick_params(
    axis='y',          # changes apply to the y-axis
    which='both',      # both major and minor ticks are affected
    left=False,      # ticks along the left edge are off
    labelleft=False) # labels along the bottom edge are off
    #plt.ylim([0,1.75])
    print('Done')
    plt.savefig('tune-mix-cross-section.png',format='png')
    plt.show()


if __name__ == "__main__":
    import sys
    main(sys.argv[1:])
