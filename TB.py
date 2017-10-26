#!/usr/bin/env python
"""\

A simple semiconductor tight-binding program. SPS* basis.

Parameters taken from Vogl, Hjalmarson and Dow, A Semiempirical Tight-Binding
Theory of the Electronic Structure of Semiconductors, J. Phys. Chem. Sol.
44 (5), pp 365-378 (1983).
"""

# The Materials Database from Harry's paper.

import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import eigvalsh
from collections import namedtuple
import itertools

Material = namedtuple('Material',['name','Esa','Epa','Esc','Epc','Essa','Essc',
                                  'Vss','Vxx','Vxy','Vsapc','Vscpa','Vssapc',
                                  'Vsscpa','bf_per_atom'])

C    = Material('C'    , -4.5450,  3.8400, -4.5450, 3.8400, 11.3700, 11.3700, 
                -22.7250, 3.8400, 11.6700, 15.2206, 15.2206, 8.2109, 8.2109,
                5)
Si   = Material('Si'   , -4.2000,  1.7150, -4.2000, 1.7150,  6.6850,  6.6850,  
                -8.3000, 1.7150,  4.5750,  5.7292,  5.7292, 5.3749, 5.3749,
                5)
Ge   = Material('Ge'   , -5.8800,  1.6100, -5.8800, 1.6100,  6.3900,  6.3900,  
                -6.7800, 1.6100,  4.9000,  5.4649,  5.4649, 5.2191, 5.2191,
                5)
Sn   = Material('Sn'   , -5.6700,  1.3300, -5.6700, 1.3300,  5.9000,  5.9000,  
                -5.6700, 1.3300,  4.0800,  4.5116,  4.5116, 5.8939, 5.8939,
                5)
SiC  = Material('SiC'  , -8.4537,  2.1234, -4.8463, 4.3466,  9.6534,  9.3166, 
                -12.4197, 3.0380,  5.9216,  9.4900,  9.2007, 8.7138, 4.4051,
                5)
AlP  = Material('AlP'  , -7.8466,  1.3169, -1.2534, 4.2831,  8.7069,  7.4231,  
                -7.4535, 2.3749,  4.8378,  5.2451,  5.2775, 5.2508, 4.6388,
                5)
AlAs = Material('AlAs' , -7.5273,  0.9833, -1.1627, 3.5867,  7.4833,  6.7267,  
                -6.6642, 1.8780,  4.2919,  5.1106,  5.4965, 4.5216, 4.9950,
                5)
AlSb = Material('AlSb' , -6.1714,  0.9807, -2.0716, 3.0163,  6.7607,  6.1543,  
                -5.6448, 1.7199,  3.6648,  4.9121,  4.2137, 4.3662, 3.0739,
                5)
GaP  = Material('GaP'  , -8.1124,  1.1250, -2.1976, 4.1150,  8.5150,  7.1850,  
                -7.4709, 2.1516,  5.1369,  4.2771,  6.3190, 4.6541, 5.0950,
                5)
GaAs = Material('GaAs' , -8.3431,  1.0414, -2.6569, 3.6686,  8.5914,  6.7386,  
                -6.4513, 1.9546,  5.0779,  4.4800,  5.7839, 4.8422, 4.8077,
                5)
GaSb = Material('GaSb' , -7.3207,  0.8554, -3.8993, 2.9146,  6.6354,  5.9846,  
                -6.1567, 1.5789,  4.1285,  4.9601,  4.6675, 4.9895, 4.2180,
                5)
InP  = Material('InP'  , -8.5274,  0.8735, -1.4826, 4.0465,  8.2635,  7.0665,  
                -5.3614, 1.8801,  4.2324,  2.2265,  5.5825, 3.4623, 4.4814,
                5)
InAs = Material('InAs' , -9.5381,  0.9099, -2.7219, 3.7201,  7.4099,  6.7401,  
                -5.6052, 1.8398,  4.4693,  3.0354,  5.4389, 3.3744, 3.9097,
                5)
InSb = Material('InSb' , -8.0157,  0.6738, -3.4643, 2.9162,  6.4530,  5.9362,  
                -5.5193, 1.4018,  3.8761,  3.7880,  4.5900, 3.5666, 3.4048,
                5)
ZnSe = Material('ZnSe' ,-11.8383,  1.5072,  0.0183, 5.9928,  7.5872,  8.9928,  
                -6.2163, 3.0054,  5.9942,  3.4980,  6.3191, 2.5891, 3.9533,
                5)
ZnTe = Material('ZnTe' , -9.8150,  1.4834,  0.9350, 5.2666,  7.0834,  8.2666,  
                -6.5765, 2.7951,  5.4670,  5.9827,  5.8199, 1.3196, 0.0000,
                5)

def phase_factors(kxyz):
    from math import sin,cos,pi
    kxp = kxyz[0]*pi/2
    kyp = kxyz[1]*pi/2
    kzp = kxyz[2]*pi/2
    g0 =  cos(kxp)*cos(kyp)*cos(kzp) - sin(kxp)*sin(kyp)*sin(kzp)*1j
    g1 = -cos(kxp)*sin(kyp)*sin(kzp) + sin(kxp)*cos(kyp)*cos(kzp)*1j
    g2 = -sin(kxp)*cos(kyp)*sin(kzp) + cos(kxp)*sin(kyp)*cos(kzp)*1j
    g3 = -sin(kxp)*sin(kyp)*cos(kzp) + cos(kxp)*cos(kyp)*sin(kzp)*1j
    return g0,g1,g2,g3

def get_Ha(id): return np.diag([id.Esa, id.Epa, id.Epa, id.Epa, id.Essa])
def get_Hc(id): return np.diag([id.Esc, id.Epc, id.Epc, id.Epc, id.Essc])

def get_Hac(id,kxyz):
    g0,g1,g2,g3 = phase_factors(kxyz)
    return np.array([
        [id.Vss*g0,    id.Vscpa*g1, id.Vscpa*g2, id.Vscpa*g3,  0],
        [-id.Vsapc*g1, id.Vxx*g0,   id.Vxy*g3,   id.Vxy*g2,   -id.Vssapc*g1],
        [-id.Vsapc*g2, id.Vxy*g3,   id.Vxx*g0,   id.Vxy*g1,   -id.Vssapc*g2],
        [-id.Vsapc*g3, id.Vxy*g2,   id.Vxy*g1,   id.Vxx*g0,   -id.Vssapc*g3],
        [0,            id.Vsscpa*g1,id.Vsscpa*g2,id.Vsscpa*g3,0]])


def get_H(id,kxyz):
    Hac = get_Hac(id,kxyz)
    # The Julia/Matlab [Ha Hac;Hac Hc] is better still    
    return np.bmat([[get_Ha(id),Hac],[Hac.conjugate(),get_Hc(id)]])

def get_kpoints(n):
    L = (0.5,0.5,0.5)
    G = (0,0,0)
    X = (1,0,0)
    K = (1,1,0)
    LG = kinterpolate(L,G,n+1) 
    GX = kinterpolate(G,X,n+1) 
    KG = kinterpolate(K,G,n+1)
    return LG+GX[1:]+KG

def kinterpolate(k1,k2,n):
    return [(k2[0]*i+k1[0]*(1-i),k2[1]*i+k1[1]*(1-i),k2[2]*i+k1[2]*(1-i))
            for i in np.linspace(0,1,n)]

def band_labels(n):
    plt.axvline(x=n,color='k')
    plt.axvline(x=2*n,color='k')
    plt.axvline(x=2*n+1,color='k')
    plt.xticks((0,n,2*n,2*n+1,3*n+1),('L',r"$\Gamma$",'X','U',r"$\Gamma$"))
    return

def bandpts(id,n=25):
    data = []
    kpts = get_kpoints(n)
    for kxyz in kpts:
        H = get_H(id,kxyz)
        E = eigvalsh(H).real
        data.append(E)
    return np.array(data)
    
def bandplt(name,data,figsize=(8,6)):
    nk,nplot = data.shape
    n = (nk-2)//3
    plt.figure(figsize=figsize)
    for i in range(nplot):
        plt.plot(data[:,i])
    band_labels(n)
    plt.axis(xmax=3*n+1)
    plt.title("Band structure for %s" % name)
    plt.ylabel("E (eV)")
    plt.savefig("%s-band-harry.png" % name)
    return
    
def band(id,n=25): bandplt(id.name,bandpts(id,n))

    
if __name__ == '__main__':
    print phase_factors((0,0,0))
