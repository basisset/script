mass={}
mass['Br']=80.0
mass['I']=127.0
mass['S']=32.0 #065
mass['O']=16.0 #994
mass['N']=14.0 #007
mass['C']=12.0 #011
mass['H']=1.00

def get_mass(frag):
    total_mass=0.0
    for i, c in enumerate(frag):#takes only one letter at a time, fails for Br
        print(f'mass_spec i,c: {i},{c}')
        if ('^' not in c):
            if c.isalnum():
                if c.isalpha():
 #            print(i,c,mass[c])
                    mass_of_current_ion=mass[c]
                    total_mass = total_mass + mass_of_current_ion
#            print("Ion: ", c)
#            print(mass_of_current_ion)
                else:
#            print(i,'multiplicity: ',c)
                    total_mass = total_mass + mass_of_current_ion*(float(c)-1.0) #
#            print("Multiplicity: ", c)
        else:
            break    
    return total_mass

def get_charge(frag):
    charge=1.0
    for i, c in enumerate(frag):
        if ('^' in c):
            for j in range(i,len(frag)):
                if frag[j].isalnum():
                    charge=float(frag[j])
                    break
        else:
            continue
        break
    return charge
    
def m_over_z(frag):
    m=get_mass(frag) 
    z=get_charge(frag) 
    return m/z

#F1='COOH'
#F2='S2COHC3H2'
#full='COOHCHNH3CH2S2C2CHNH2COOH'
#F1='NH2CCH2SSCH2CHNH2'
#F1='NH2CCH2SSCH2CHNH2'
#F1='NH2CHCH2SSCH2CHNH2'
#F2='NH2CH2COOH'
#F3='SS'
#F4='CH2CHCO'
#F5='CH2SSCH2'
#F6='C4S2N2H8'
#F7='CH2CHNH2COOH'
#F8='CH2CHCOOH'
#F9='NH$_2$CH$_2$COOH$^{+}$'
#F10="NH$_2$CHCH$_2$SSCH$_2$CHNH$_2^{2+}"



#a = get_mass(full)
#b = get_mass(F2)
#print('full: ', full, get_mass(full))
#print('F1: ', F1, get_mass(F1))
#print('F2: ', F2, get_mass(F2))
#print('F3: ', F3, get_mass(F3))
##print('F4: ', F4, get_mass(F4))
##print('F5: ', F5, get_mass(F5))
#print('F6: ', F6, get_mass(F6))
#print('F7: ', F7, get_mass(F7))
#print('F8: ', F8, get_mass(F8))
#print('F9: ', F9, m_over_z(F9))
#print('F10: ', F10, m_over_z(F10))
#print('F11: ', F11, get_mass(F11))
#print('F12: ', F12, get_mass(F12))
#print('F13: ', F13, get_mass(F13))
#print('F14: ', F14, get_mass(F14))
#print('F15: ', F15, get_mass(F15))
#print('F16: ', F16, get_mass(F16))
