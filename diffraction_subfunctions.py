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