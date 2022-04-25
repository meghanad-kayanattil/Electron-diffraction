# -*- coding: utf-8 -*-
"""
Created on Wed Apr 13 13:15:47 2022

@author: kayanatm
"""

import math
import numpy as np
import pyFAI.detectors
from pyFAI.azimuthalIntegrator import AzimuthalIntegrator
from sklearn import metrics


def create_s_matrix(xcen = 391.55, ycen = 393.03):
    pixel_size_det = 4.8e-6#size of the pixel
    shape_det = [900,900]
    distanceDet = 0.5 # in m
    
    KE = 50 #kinetic energy of electrons kev
    planck_c = 4.135667516e-15 #planck constant
    m0 = 9.10938215e-31 #rest mass of electron
    speed_light = 299792458 #Speed of light
    wavel = (planck_c*speed_light/(math.sqrt(KE*1000*(KE*1000+2*m0*speed_light**2*6.241506363e+18))));
    

    x = np.arange(1,shape_det[0]+1)
    y = np.arange(1,shape_det[1]+1)
    [xx,yy] = np.meshgrid(x,y)

    #Calculates the distance of each px to the center of the image in meters and in s units.
    d = []

    for x,y in zip(xx,yy):
        temp_holder = ((x-xcen)**2)+((y-ycen)**2)

        d.append(np.sqrt(temp_holder)*pixel_size_det)

    theta = np.arctan(np.divide(d,distanceDet))

    four_pi_lam  = 4*math.pi/wavel             

    s = np.multiply(np.sin(np.divide(theta,2)),four_pi_lam)*1e-9
    
    return s

def rmatrix(xcen = 391.55, ycen = 393.03, size = 900):
    
    x = np.arange(1,size+1)
    y = np.arange(1,size+1)
    [xx,yy] = np.meshgrid(x,y)

    #Calculates the distance of each px to the center of the image in meters and in s units.
    d = []

    for x,y in zip(xx,yy):
        temp_holder = ((x-xcen)**2)+((y-ycen)**2)

        d.append(np.sqrt(temp_holder))
    
    return d

def diffraction_to_epdf(img, KE=50):
    #Defining experiment parameters
    #KE = 50 #kinetic energy of electrons kev
    planck_c = 4.135667516e-15 #planck constant
    m0 = 9.10938215e-31 #rest mass of electron
    speed_light = 299792458 #Speed of light
    lambda_e = (planck_c*speed_light/(math.sqrt(KE*1000*(KE*1000+2*m0*speed_light**2*6.241506363e+18)))); 
    print("Electron beam wavelength = {} m".format(lambda_e))

    #Defining the detector parameters this can be adapted from the experiment
    xcen = 391.55 #center of the diffraction
    ycen = 393.03 #center of the diffraction
    pixel_size = 4.8e-5 #size of the pixel
    
    distanceDet = 0.5 # in m 
    wavel = lambda_e
    
    #Defining the detector in pyFAI
    

    detector = pyFAI.detectors.Detector(pixel1=pixel_size, pixel2=pixel_size,max_shape=(img.shape[0],img.shape[1]))
    ai = AzimuthalIntegrator(dist=distanceDet, detector=detector, wavelength=wavel)
    ai.setFit2D(directDist=500,centerX=xcen,centerY=ycen)

    print(ai)
    #azimuthal integration using the intergate1d method
    #the default radial unit is “q_nm^1”, so the scattering vector length expressed in inverse nanometers. 
    #To be able to calculate q, one needs to specify the wavelength used
    res = ai.integrate1d(img, 100)

    #The ds in stuart's code is calculated as ds = (2*np.pi*pixelSize*radialDist[0])/(wavel*distanceDet)
    #here I believe that is done by 2*pi*pixel_size/(wavel/distanceDet)
    s = res[0]*0.1 
    sMs = res[1]
    
    return s, sMs

def function_inside(s,sMs,r):
    k = 0.02
    s2 = np.multiply(s,s)
    sinfun = np.sin(np.multiply(s,r))
    expfun = np.exp(np.multiply(-k,s2))
    out = sMs*np.multiply(sinfun,expfun)
    return out

def sms_to_mrdf(s,sMs):
    r = np.linspace(0,8, 100)
    mrdf = []


    for ri in r:
        y = function_inside(s, sMs, ri)
        area = metrics.auc(s,y)
        mrdf.append(area)
    return mrdf

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
        elif xyz.iloc[i,0] == 'Cl':
            z.append(17)
    return z

def f_x_kirk(s,zi):
    """Scattering factor calculation kirkland method
    
    
    """
    p = kirk_scat(zi)
    q = np.divide(s,math.pi*2)
    
    q_sq = np.multiply(q,q)
    
    f = np.divide(p['a1'],q_sq+p['b1'])+np.divide(p['a2'],q_sq+p['b2'])+np.divide(p['a3'],q_sq+p['b3'])
    
    f = f + p['c1']*np.exp(-p['d1']*q_sq) + p['c2']*np.exp(-p['d2']*q_sq) + p['c3']*np.exp(-p['d3']*q_sq)
    
    return f



def scattering_int(xyz,dim,N_atoms,s):
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
 
    for i in range(N_atoms):
        for j in range(i+1):
            zi = Z[i]
            zj = Z[j]
            #print('Zi = {}, Zj = {}'.format(zi,zj))
            
            #fi,fj are the elastic scattering amplitude of atoms
            fi = f_x_kirk(s,zi)
            fj = f_x_kirk(s,zj)
            
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
            print('Pairs = {}'.format(pi))
    return IAt, Imol