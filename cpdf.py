# -*- coding: utf-8 -*-
"""
Created on Thu Feb 24 13:30:58 2022

@author: kayanatm
"""

import numpy as np
import math
import os
import pandas as pd
import matplotlib.pyplot as plt

#------------------------------------------------------------------------------
#Geometry file handling
xyz_file_path = os.getcwd()
xyz_file_name ='Glycerol.xyz'
xyz = pd.read_csv(os.path.join(xyz_file_path,xyz_file_name), sep = '\s+', dtype='str',names = ['Atoms', 'x', 'y', 'z'])
xyz.iloc[0,0]
xyz[["x", "y", "z"]] = xyz[["x", "y", "z"]].apply(pd.to_numeric)

#------------------------------------------------------------------------------

#Constants and other important parameters
KE = 50 #kinetic energy of electrons kev
planck_c = 4.135667516e-15 #planck constant
m0 = 9.10938215e-31 #rest mass of electron
speed_light = 299792458 #Speed of light
lambda_e = (planck_c*speed_light/(math.sqrt(KE*1000*(KE*1000+2*m0*speed_light**2*6.241506363e+18)))); 
#Defining the detector parameters this can be adapted from the experiment
xcen = 391.55 #center of the diffraction
ycen = 393.03 #center of the diffraction
pixel_size_det = 4.8e-6#size of the pixel
shape_det = [900,900]
distanceDet = 0.5 # in m
wavel = lambda_e
print('Wave length of electrons = {} m'.format(wavel)) #for 50kv electrons
N_atoms = len(xyz); print('Number of atoms in the molecule = {}'.format(N_atoms))
pairs = int(N_atoms*(N_atoms-1)/2); print('Number of pairs of atomic connections in the molecule = {}'.format(pairs)); 
r = np.zeros(pairs)
#------------------------------------------------------------------------------
#s-matrix creation
sM = np.zeros(shape_det)
x = np.arange(1,shape_det[0])
y = np.arange(1,shape_det[1])
[xx,yy] = np.meshgrid(x,y)
#Calculates the distance of each px to the center of the image in meters and in s units.
d = []

for x,y in zip(xx,yy):
    temp_holder = ((x-xcen)**2)+((y-ycen)**2)

    d.append(np.sqrt(temp_holder)*pixel_size_det)

theta = np.arctan(np.divide(d,distanceDet))

four_pi_lam  = 4*math.pi/wavel             

s = np.multiply(np.sin(np.divide(theta,2)),four_pi_lam)*1e-10

print("Minimum s value = {}".format(s.min()))
print("Maximum s value = {}".format(s.max()))
#--------------------------------------------------------------------------------

#Dependant Functions
#------------------------------------------------------------------------------
def lobato_scat(zi):
    """The function takes in the atomic number and returns the scattering (Lobato) function constant values for 
        that specific atom
        Input = Z (Atomic number)
        Output = Scattering vector dictionary with constants (form factors"""
    lobato = {'8':['O', [2.994740452423624e+01, -7.761012662552783e+01, 9.988177646231442e+01, -5.121270055056731e+01, 8.196189544460320e-03],[1.302839878800107e+00, 1.157941052583095e+00, 1.009885493380251e+00, 9.433279714332660e-01, 4.331976113218256e-02]],
              '6':['C', [1.244660886213433e+02, -2.203528570789638e+02,1.952353522804791e+02,-9.810793612697997e+01,1.420230412136232e-02],[2.421208492560056e+00, 2.305379437524258e+00, 2.048519321065642e+00, 1.933525529175474e+00, 7.689768184783397e-02]],
              '1':['H',[6.473848488352918e-03, -4.901925767802290e-01, 5.732841603908765e-01, -3.794033014839905e-01, 5.544264747740791e-01], [2.785198853791489e+00, 2.776204283306448e+00, 2.775385910506251e+00, 2.767593028672588e+00, 2.765118976429275e+00]]
        
    }
    
    
    a_vals = lobato[str(zi)][1]
    b_vals = lobato[str(zi)][2]
    
    return a_vals, b_vals
    
def f_x_lobato(s,zi):
    "Scattering factor calculation Lobato method"
    
    a,b = lobato_scat(zi)
    
    q = np.divide(s,math.pi*2)
    
    q_sq = np.multiply(q,q)
    
    f=0
    for i in range(5):
        f = f + (a[i]*(2+b[i]*q_sq))/(1+b[i]*q_sq)**2
    
    return f

def scattering_int_lobato(xyz,dim,N_atoms,s):
    """
    Calculates the scattering intensity
    Input parameters: xyz = geometry data
                      dim = detector size in pixels 
                      N_atoms = Number of atoms in the molecule
                      s = scattering vector calculated in the main program body
    output: Iat = Atomic scattering intensity
            Imol = Molecular scattering intensity
    
    """
    IAt = np.zeros(dim)
    Imol = np.zeros(dim)

    pi = 0
    Z = find_z(xyz)
    print('Scattering simulation running....')
 
    for i in range(N_atoms):
        for j in range(i+1):
            zi = Z[i]
            zj = Z[j]
            #print('Zi = {}, Zj = {}'.format(zi,zj))
            
            #fi,fj are the elastic scattering amplitude of atoms
            fi = f_x_lobato(s,zi)
            fj = f_x_lobato(s,zj)
            
            if (i == j):
                IAt = IAt + np.multiply(fi,fj)
                #print('rij = 0')
                
            else:
                #We can divide the Imol(s) into 3 parts, first the scatering factor
                #second the cos (phase) factor and 3rd the sin factor
                
                #scat factor calculation
                scat_fact = np.multiply(fi,fj)
                
                #sin factor calculation 
                r1 = xyz.iloc[i,1]-xyz.iloc[j,1]
                r2 = xyz.iloc[i,2]-xyz.iloc[j,2]
                r3 = xyz.iloc[i,3]-xyz.iloc[j,3]
                r_squared = r1**2+r2**2+r3**2
                rij = math.sqrt(r_squared);#print('rij = {}'.format(rij))
                sin_fact = np.divide(np.sin(np.multiply(rij,s)),np.multiply(rij,s))            
                
                #cos_factor calculation 
                
                #total
                Imol = Imol + np.multiply(scat_fact,sin_fact)
            pi=pi+1        
            #print('Pairs = {}'.format(pi))
    return IAt, Imol               

def find_z(xyz):
    """ This functon retuns the atomic numbers of the elements in the molecule into a list
        (Currently only functional for H,C and O)
    """
    z = []
    for i in range(len(xyz)):
        if xyz.iloc[i,0] == 'H':
            z.append(1)
        elif xyz.iloc[i,0] == 'C':
            z.append(6)
        elif xyz.iloc[i,0] == 'O':
            z.append(8)
    return z


#---------------------------------------------------------------------------
#Execution and plotting
Iatom_loba, Imol_loba = scattering_int_lobato(xyz=xyz, dim = [899,899], s=s,N_atoms=14);
Itot_loba = Imol_loba + Iatom_loba
sMs = np.multiply(np.divide(Imol_loba,Iatom_loba),s)

plt.imshow(Itot_loba)
plt.colorbar()
plt.imsave(fname='Itot_glycerol_isolated.png', arr=Itot_loba, cmap='gray', format='png')
#-------------------------------------------------------------------------------