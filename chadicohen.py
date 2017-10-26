#!/usr/bin/env python

# chadicohen.py Richard P. Muller, 12/99

# This program is the tight-binding program for Diamond/Zincblende
# structures that is presented in Chadi and Cohen's paper
# "Tight-Binding Calculations of the Valence Bands of Diamond and
# Zincblende Crystals", Phys. Stat. Soli. (b) 68, 405 (1975).  This
# program is written for diamond and zincblende structures only.

# Copyright 1999, Richard P. Muller and William A. Goddard, III
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# Here are some sample band gaps (from Kittel) that may aid in fitting:
# C    i 5.4       GaAs d 1.52 
# Si   i 1.17      GaP  i 2.32 
# Ge   i 0.744     GaN    3.5  (not from Kittel) 
# Sn   d 0.00      InP  d 1.42 
#                  InAs d 0.43

import sys,getopt,os 
from Numeric import zeros,Complex        # Make sure that NumPy is in your
from LinearAlgebra import eigenvalues    #  PYTHONPATH variable
from math import sqrt,pi,cos,sin

# This should be the only variable you have to change:
#path_to_gnuplot = 'c:\Gnuplot3.7\wgnuplot.exe' # On my windows98 box
#path_to_gnuplot = '/usr/local/bin/gnuplot'
path_to_gnuplot = 'gnuplot'

def error_and_exit(line):
    print line
    sys.exit()
    return

def help_and_exit():
    help()
    sys.exit()
    return

def help():
    print "chadicohen.py: Tight-binding band structure of II-VI, III-V,"
    print "and IV semiconductors. Based on Chadi/Cohen's approach."
    print ""
    print "usage: chadicohen.py [options]"
    print ""
    print "Options:"
    print "-n #  The number of points in each Brillouin zone region "
    print "      (default=10)"
    print "-h    Print this help screen and exit"
    print "-s    Structure to compute; currently supported are:"
    print "      C    Diamond"
    print "      Si   Silicon"
    print "      Ge   Germanium"
    print "      GaAs Gallium Arsenide"
    print "      ZnSe Zinc Selenide"
    print "-P    output a postscript version of the plot"
    print "-G    output a GIF of the plot"
    print ""
    print "Caveats:"
    print "(1) The parameters in the code are simply taken from Chadi/Cohen."
    print "    No checking is performed to make sure that they work for "
    print "    the case of interest"
    print "(2) This program assumes that Gnuplot is installed, and is "
    print "    started by the command \"gnuplot\". If this isn't the "
    print "    case on your system edit path_to_gnuplot accordingly"
    print "(3) This program assumes that /usr/bin/env python can find"
    print "    python on your system. If not, edit the first line of this"
    print "    file accordingly."
    print "(4) This program assumes that the Numeric Extensions to Python"
    print "    (see ftp://ftp-icf.llnl.gov/pub/python) are installed,"
    print "    and are in your $PYTHONPATH."
    print ""
    print "References:"
    print "D.J. Chadi and M.L. Cohen, \"Tight Binding Calculations"
    print "of the Valence Bands of Diamond and Zincblende Crystals.\""
    print "Phys. Stat. Sol. (b) 68, 405 (1975)" 
    return
def get_k_points(n):
    # Define a set of k points along the Brillouin zone boundary from
    # L (0.5,0.5,0.5) to Gamma (0,0,0) to X (1.,0.,0.). Space the points
    # evenly, based on the scaling parameter n.
    kpoints = []
    step = 0.5/float(n)
    kx,ky,kz = 0.5,0.5,0.5    # Start at the L point (1/2,1/2,1/2)
    kpoints.append((kx,ky,kz))
    for i in range(n): # Move to the Gamma point (0,0,0)
        kx,ky,kz = kx-step,ky-step,kz-step
        kpoints.append((kx,ky,kz))
    for i in range(n): # Now go to the X point (1,0,0)
        kx = kx+2.*step
        #ky = ky+2.*step        
        kpoints.append((kx,ky,kz))
    kx = ky = 1.0 # Jump to the U,K point
    kz = 0.0
    kpoints.append((kx,ky,kz))
    for i in range(n): # Now go back to Gamma
        kx = kx - 2.*step
        ky = ky - 2.*step
        kpoints.append((kx,ky,kz))
    return kpoints        

def sort_eigenvalues(E):
    # This is trickier than it sounds, since NumPy doesn't define
    # sort on Complex numbers. Convert to a normal python array of
    # reals, and sort.
    enarray = []
    for en in E: enarray.append(en.real)
    enarray.sort()
    return enarray

# ---------------Top of main program------------------
# program defaults:

n = 10
structure = 'Si'

# Get command line options:
opts, args = getopt.getopt(sys.argv[1:],'nhs:PG')
postscript = 0
gif = 0
for (key,value) in opts:
    if key == '-n':   n = eval(value)
    if key == '-h':   help_and_exit()
    if key == '-s':   structure = value
    if key == '-P':   postscript = 1
    if key == '-G':   gif = 1

# K points (these must be multiplied by 2*pi/a)
kpoints = get_k_points(n)

# Tight binding parameters; these are in eV:
if structure == 'C':
    e_s_c = e_s_a = 0.0     # Arbitrary;
    e_p_c = e_p_a = 7.40 - e_s_c
    v_ss = -15.2
    v_sc_p = v_sa_p = 10.25
    v_xx = 3.0
    v_xy = 8.30
elif structure == 'Si':
    e_s_c = e_s_a = 0.0     # Arbitrary
    e_p_c = e_p_a = 7.20 - e_s_c
    v_ss = -8.13
    v_sc_p = v_sa_p = 5.88
    v_xx = 3.17
    v_xy = 7.51
elif structure == 'Ge':
    e_s_c = e_s_a = 0.0     # Arbitrary
    e_p_c = e_p_a = 8.41 - e_s_c
    v_ss = -6.78
    v_sc_p = v_sa_p = 5.31
    v_xx = 2.62
    v_xy = 6.82
elif structure == 'GaAs':
    e_s_c = -6.01
    e_s_a = -4.79
    e_p_c = 0.19
    e_p_a = 4.59
    v_ss = -7.00
    v_sc_p = 7.28
    v_sa_p = 3.70
    v_xx = 0.93
    v_xy = 4.72
elif structure == 'ZnSe':
    e_s_c = -8.92
    e_s_a = -0.28
    e_p_c = 0.12
    e_p_a = 7.42
    v_ss = -6.14
    v_sc_p = 5.47
    v_sa_p = 4.73
    v_xx = 0.96
    v_xy = 4.38
else:
    error_and_exit('Program can\'t cope with structure %s' % structure)

gfile = open('chadicohen.tmp','w')

for (kx,ky,kz) in kpoints:
    kxp,kyp,kzp = kx*pi/2.,ky*pi/2.,kz*pi/2.# The a's cancel here

    g0_real = cos(kxp)*cos(kyp)*cos(kzp)
    g0_imag = -sin(kxp)*sin(kyp)*sin(kzp)
    g1_real = -cos(kxp)*sin(kyp)*sin(kzp)
    g1_imag = sin(kxp)*cos(kyp)*cos(kzp)
    g2_real = -sin(kxp)*cos(kyp)*sin(kzp)
    g2_imag = cos(kxp)*sin(kyp)*cos(kzp)
    g3_real = -sin(kxp)*sin(kyp)*cos(kzp)
    g3_imag = cos(kxp)*cos(kyp)*sin(kzp)
    
    # "s" stands for "star": the complex conjugate
    g0,g0s = g0_real+g0_imag*1j,g0_real-g0_imag*1j
    g1,g1s = g1_real+g1_imag*1j,g1_real-g1_imag*1j
    g2,g2s = g2_real+g2_imag*1j,g2_real-g2_imag*1j
    g3,g3s = g3_real+g3_imag*1j,g3_real-g3_imag*1j
    
    H = zeros((8,8),Complex)

    # Make the diagonal elements
    H[0,0] = e_s_c
    H[1,1] = e_s_a
    H[2,2] = H[3,3] = H[4,4] = e_p_c
    H[5,5] = H[6,6] = H[7,7] = e_p_a

    # Make the off-diagonal parts
    H[1,0] = v_ss*g0s
    H[0,1] = v_ss*g0

    H[2,1] = -v_sa_p*g1
    H[1,2] = -v_sa_p*g1s
    H[3,1] = -v_sa_p*g2
    H[1,3] = -v_sa_p*g2s
    H[4,1] = -v_sa_p*g3
    H[1,4] = -v_sa_p*g3s

    H[5,0] = v_sc_p*g1s
    H[0,5] = v_sc_p*g1
    H[6,0] = v_sc_p*g2s
    H[0,6] = v_sc_p*g2
    H[7,0] = v_sc_p*g3s
    H[0,7] = v_sc_p*g3

    H[5,2] = v_xx*g0s
    H[2,5] = v_xx*g0
    H[6,2] = v_xy*g3s
    H[2,6] = v_xy*g3
    H[7,2] = v_xy*g2s
    H[2,7] = v_xy*g2

    H[5,3] = v_xy*g3s
    H[3,5] = v_xy*g3
    H[6,3] = v_xx*g0s
    H[3,6] = v_xx*g0
    H[7,3] = v_xy*g1s
    H[3,7] = v_xy*g1

    H[5,4] = v_xy*g2s
    H[4,5] = v_xy*g2
    H[6,4] = v_xy*g1s
    H[4,6] = v_xy*g1
    H[7,4] = v_xx*g0s
    H[4,7] = v_xx*g0

    enarray = eigenvalues(H)
    enarray = sort_eigenvalues(enarray)
    for en in enarray:
        gfile.write("%13.8f" % en)
    gfile.write("\n")

gfile.close()

label_y = 0.05
L_x = 0
gamma_x = n
X_x = 2*n
K_x = 2*n+1
gamma2_x = 3*n+1

pfile = open('chadicohen.gnu','w')
pfile.write('set data style linespoints\n')
pfile.write('set noxtics\n')
pfile.write('set ylabel "E(eV)"\n')
pfile.write('set xlabel "k points"\n')
pfile.write('set label "L" at %d,graph %f\n' % (L_x,label_y))
pfile.write('set label "G" at %d,graph %f\n' % (gamma_x,label_y))
pfile.write('set label "X" at %d,graph %f\n' % (X_x,label_y))
pfile.write('set label "K" at %d,graph %f\n' % (K_x,label_y))
pfile.write('set label "G" at %d,graph %f\n' % (gamma2_x,label_y))
pfile.write('set arrow from %d,graph 0 to %d, graph 1 nohead\n' %\
            (gamma_x,gamma_x))
pfile.write('set arrow from %d,graph 0 to %d, graph 1 nohead\n' %\
            (X_x,X_x))
pfile.write('set arrow from %d,graph 0 to %d, graph 1 nohead\n' %\
            (K_x,K_x))
pfile.write('set nokey\n') 
pfile.write('set title "Band Structure for cubic %s"\n' % structure)
pfile.write('plot ')
for i in range(8):
    pfile.write('"chadicohen.tmp" using %d' % (i+1))
    if i != 7:
        pfile.write(',')
    else:
        pfile.write('\n')
pfile.write('pause -1\n')
if postscript:
    pfile.write('set term post\n')
    pfile.write('set output "chadicohen.ps"\n')
    pfile.write('replot\n')
if gif:
    pfile.write('set term gif\n')
    pfile.write('set output "chadicohen.gif"\n')
    pfile.write('replot\n')
pfile.close()
os.system('%s chadicohen.gnu' % path_to_gnuplot) 
