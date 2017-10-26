#!/usr/bin/env python

# harrison.py Richard P. Muller, 12/99

# This program is the tight-binding program for Diamond/Zincblende
# structures that is presented in Harrison's "Electronic Structure and
# the Properties of Solids," Dover, 1989, p. 77.  This program is
# written for diamond and zincblende structures only.

# Copyright 1999, Richard P. Muller and William A. Goddard, III
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# General program notes:
# The band gap here is greatly over-estimated (by almost exactly a
# factor of two, which suggests an obvious workaround).
# No one seems to publish tight-binding band gaps, so maybe they all
# do a bad job. The width of the
# valence band compares well with Chadi and Cohen. They, however, see
# much more dispersion. I suppose this is an exercise in fitting.

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

def help_and_exit():
    help()
    sys.exit()
    return

def help():
    print "harrison.py: Tight-binding band structure of II-VI, III-V,"
    print "and IV semiconductors. Based on Harrison's version of "
    print "Chadi/Cohen's approach."
    print ""
    print "usage: harrison.py [options]"
    print ""
    print "Options:"
    print "-c #  The atomic number for the semiconductor cation (default=14)"
    print "-a #  The atomic number for the semiconductor anion (default=14)"
    print "-n #  The number of points in each Brillouin zone region "
    print "      (default=10)"
    print "-h    Print this help screen and exit"
    print "-P    output a postscript version of the plot"
    print "-G    output a GIF of the plot"
    print ""
    print "Caveats:"
    print "(1) The parameters in the code are simply taken from Harrison."
    print "    No checking is performed to make sure that they work for "
    print "    the case of interest"
    print "(2) Similarly, no checking is done to insure that the species"
    print "    you input make sense in a zincblende structure. I.e., you could"
    print "    input GaAl and the program would give you a (meaningless)"
    print "    answer"
    print "(3) This program assumes that Gnuplot is installed, and is "
    print "    started by the command \"gnuplot\". If this isn't the "
    print "    case on your system edit path_to_gnuplot accordingly"
    print "(4) This program assumes that /usr/bin/env python can find"
    print "    python on your system. If not, edit the first line of this"
    print "    file accordingly."
    print "(5) This program assumes that the Numeric Extensions to Python"
    print "    (see ftp://ftp-icf.llnl.gov/pub/python) are installed,"
    print "    and are in your $PYTHONPATH."
    print ""
    print "References:"
    print "D.J. Chadi and M.L. Cohen, \"Tight Binding Calculations"
    print "of the Valence Bands of Diamond and Zincblende Crystals.\""
    print "Phys. Stat. Sol. (b) 68, 405 (1975)" 
    print ""
    print "W.A. Harrison, \"Electronic Structure and the Properties"
    print "of Solids: The Physics of the Chemical Bond.\""
    print "Dover Publications, Inc., NY, 1989"
    return

# Parameters from Harrison pp. 552-553. Zeroes indicate that the program
# shouldn't work with these values and that it should stop
e_s = [
    0.0, 0.0, 0.0, 0.0, 8.17,          # 4
    12.54, 17.52, 23.04, 29.14, 0.0,   # 9
    0.0, 0.0, 6.86, 10.11, 13.55,      # 14
    17.10, 20.80, 0.0, 0.0, 0.0,       # 19
    0.0, 0.0, 0.0, 0.0, 0.0,           # 24
    0.0, 0.0, 0.0, 0.0, 0.0,           # 29
    8.40, 11.37, 14.38, 17.33, 20.32,  # 34
    0.0, 0.0, 0.0, 0.0, 0.0,           # 39
    0.0, 0.0, 0.0, 0.0, 0.0,           # 44
    0.0, 0.0, 0.0, 7.70, 10.12,        # 49
    12.50, 14.80, 17.11, 0.0, 0.0,     # 54
    0.0, 0.0, 0.0, 0.0, 0.0,           # 59
    0.0, 0.0, 0.0, 0.0, 0.0,           # 64
    0.0, 0.0, 0.0, 0.0, 0.0,           # 69
    0.0, 0.0, 0.0, 0.0, 0.0,           # 74
    0.0, 0.0, 0.0, 0.0, 0.0,           # 79
    7.68]                              # 80 (Hg)

e_p = [
    0.0, 0.0, 0.0, 0.0, 4.14,          # 4
    6.64, 8.97, 11.47, 14.13, 0.0,     # 9
    0.0, 0.0, 2.99, 4.86, 6.52,        # 14
    8.33, 10.27, 0.0, 0.0, 0.0,        # 19
    0.0, 0.0, 0.0, 0.0, 0.0,           # 24
    0.0, 0.0, 0.0, 0.0, 0.0,           # 29
    3.38, 4.90, 6.36, 7.91, 9.53,      # 34
    0.0, 0.0, 0.0, 0.0, 0.0,           # 39
    0.0, 0.0, 0.0, 0.0, 0.0,           # 44
    0.0, 0.0, 0.0, 3.38, 4.69,         # 49
    5.94, 7.24, 8.59, 0.0, 0.0,        # 54
    0.0, 0.0, 0.0, 0.0, 0.0,           # 59
    0.0, 0.0, 0.0, 0.0, 0.0,           # 64
    0.0, 0.0, 0.0, 0.0, 0.0,           # 69
    0.0, 0.0, 0.0, 0.0, 0.0,           # 74
    0.0, 0.0, 0.0, 0.0, 0.0,           # 79
    3.48]                              # 80

r0 = [
    0.0, 0.0, 0.0, 0.0, 1.54,          # 4
    1.54, 1.54, 1.54, 1.54, 0.0,       # 9
    0.0, 0.0, 2.35, 2.35, 2.35,        # 14
    2.35, 2.35, 0.0, 0.0, 0.0,         # 19
    0.0, 0.0, 0.0, 0.0, 0.0,           # 24
    0.0, 0.0, 0.0, 0.0, 0.0,           # 29
    2.44, 2.44, 2.44, 2.44, 2.44,      # 34
    0.0, 0.0, 0.0, 0.0, 0.0,           # 39
    0.0, 0.0, 0.0, 0.0, 0.0,           # 44
    0.0, 0.0, 0.0, 2.80, 2.80,         # 49
    2.80, 2.80, 2.80, 0.0, 0.0,        # 54
    0.0, 0.0, 0.0, 0.0, 0.0,           # 59
    0.0, 0.0, 0.0, 0.0, 0.0,           # 64
    0.0, 0.0, 0.0, 0.0, 0.0,           # 69
    0.0, 0.0, 0.0, 0.0, 0.0,           # 74
    0.0, 0.0, 0.0, 0.0, 0.0,           # 79
    3.00]                              # 80 (this parameter was fudged by RPM)

label = [
    '', '', '', '', 'Be',         # 4
    'B', 'C', 'N', 'O', '',       # 9
    '', '', 'Mg', 'Al', 'Si',     # 14
    'P', 'S', '', '', '',         # 19
    '', '', '', '', '',           # 24
    '', '', '', '', '',           # 29
    'Zn', 'Ga', 'Ge', 'As', 'Se', # 34
    '', '', '', '', '',           # 39
    '', '', '', '', '',           # 44
    '', '', '', 'Cd', 'In',       # 49
    'Sn', 'Sb', 'Te', '', '',     # 54
    '', '', '', '', '',           # 59
    '', '', '', '', '',           # 64
    '', '', '', '', '',           # 69
    '', '', '', '', '',           # 74
    '', '', '', '', '',           # 79
    'Hg']                         # 80
    

# It would be faster to put these in the code; I'm defining them as
#  functions to make it easier to understand what I'm doing.
def convert_angstroms_to_bohr(d): return d/0.52918
def convert_eV_to_hartrees(en): return en*0.03675
def convert_hartrees_to_eV(en): return en*27.212
def dot(t1,t2): #dot product of two 3-tuples
    return t1[0]*t2[0]+t1[1]*t2[1]+t1[2]*t2[2]

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
    for i in range(len(enarray)):
        enarray[i] = convert_hartrees_to_eV(enarray[i])
    enarray.sort()
    return enarray

# ---------------Top of main program------------------
# program defaults:
atno_c = 14 # Note on convention: "_c" stands for "cation", and "_a"
atno_a = 14 #  stands for "anion"
n = 10

postscript = 0
gif = 0
# Get command line options:
opts, args = getopt.getopt(sys.argv[1:],'n:hc:a:PG')
for (key,value) in opts:
    if key == '-n':   n = eval(value)
    if key == '-h':   help_and_exit()
    if key == '-c':   atno_c = eval(value)
    if key == '-a':   atno_a = eval(value)
    if key == '-P':   postscript = 1
    if key == '-G':   gif = 1

# K points (these must be multiplied by 2*pi/a)
kpoints = get_k_points(n)

# Parameters (from Harrison, pp 551-552).

#Note: The nxx_* parameters were given in Harrison without units. I'm
# assuming that they are dimensionless, since when I try and convert
# them into hartrees I get garbage.
nss_sigma = -1.40
nsp_sigma = 1.84
npp_sigma = 3.24
npp_pi = -0.81

#Get element-based variables
e_s_c = -e_s[atno_c]
e_s_a = -e_s[atno_a]
e_p_c = -e_p[atno_c]
e_p_a = -e_p[atno_a]
d = 0.5*(r0[atno_c] + r0[atno_a])
label0 = label[atno_c]
label1 = label[atno_a]

# Check to make sure that all parameters are defined; zero indicates that
# this is an element that the program shouldn't work for (i.e. one that isn't
# in column II, III, IV, V, or VI).
if e_s_c == 0. or e_s_a == 0. or e_p_c == 0. or e_p_a == 0. \
   or r0[atno_c] == 0. or r0[atno_a] == 0.:
    print "Bad parameter; values are:"
    print e_s_c,e_s_a,e_p_c,e_p_a,r0[atno_c],r0[atno_a]
    sys.exit()

#convert all of these to the proper units
e_s_c = convert_eV_to_hartrees(e_s_c)
e_s_a = convert_eV_to_hartrees(e_s_a)
e_p_c = convert_eV_to_hartrees(e_p_c)
e_p_a = convert_eV_to_hartrees(e_p_a)

d = convert_angstroms_to_bohr(d)

Vss_sigma = nss_sigma/(d*d) 
Vsp_sigma = nsp_sigma/(d*d) 
Vpp_sigma = npp_sigma/(d*d) 
Vpp_pi = npp_pi/(d*d)

Ess = Vss_sigma
Esp = -Vsp_sigma/sqrt(3.)
Exx = Vpp_sigma/3. + 2.*Vpp_pi/3.
Exy = Vpp_sigma/3. - Vpp_pi/3.

d1 = (0.25,0.25,0.25)     # The d* variables are in units of a
d2 = (0.25,-0.25,-0.25)   
d3 = (-0.25,0.25,-0.25)   
d4 = (-0.25,-0.25,0.25)  

gfile = open('harrison.tmp','w')

for kpoint in kpoints:
    kx,ky,kz = kpoint
    kd1 = 2.*pi*dot(kpoint,d1) # The a's cancel here
    kd2 = 2.*pi*dot(kpoint,d2) 
    kd3 = 2.*pi*dot(kpoint,d3) 
    kd4 = 2.*pi*dot(kpoint,d4) 

    g0_real = cos(kd1) + cos(kd2) + cos(kd3) + cos(kd4)
    g0_imag = sin(kd1) + sin(kd2) + sin(kd3) + sin(kd4)
    g1_real = cos(kd1) + cos(kd2) - cos(kd3) - cos(kd4)
    g1_imag = sin(kd1) + sin(kd2) - sin(kd3) - sin(kd4)
    g2_real = cos(kd1) - cos(kd2) + cos(kd3) - cos(kd4)
    g2_imag = sin(kd1) - sin(kd2) + sin(kd3) - sin(kd4)
    g3_real = cos(kd1) - cos(kd2) - cos(kd3) + cos(kd4)
    g3_imag = sin(kd1) - sin(kd2) - sin(kd3) + sin(kd4)
    
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
    H[1,0] = Ess*g0s
    H[0,1] = Ess*g0

    H[2,1] = -Esp*g1
    H[1,2] = -Esp*g1s
    H[3,1] = -Esp*g2
    H[1,3] = -Esp*g2s
    H[4,1] = -Esp*g3
    H[1,4] = -Esp*g3s

    H[5,0] = Esp*g1s
    H[0,5] = Esp*g1
    H[6,0] = Esp*g2s
    H[0,6] = Esp*g2
    H[7,0] = Esp*g3s
    H[0,7] = Esp*g3

    H[5,2] = Exx*g0s
    H[2,5] = Exx*g0
    H[6,2] = Exy*g3s
    H[2,6] = Exy*g3
    H[7,2] = Exy*g2s
    H[2,7] = Exy*g2

    H[5,3] = Exy*g3s
    H[3,5] = Exy*g3
    H[6,3] = Exx*g0s
    H[3,6] = Exx*g0
    H[7,3] = Exy*g1s
    H[3,7] = Exy*g1

    H[5,4] = Exy*g2s
    H[4,5] = Exy*g2
    H[6,4] = Exy*g1s
    H[4,6] = Exy*g1
    H[7,4] = Exx*g0s
    H[4,7] = Exx*g0

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

pfile = open('harrison.gnu','w')
pfile.write('set data style linespoints\n')
pfile.write('set noxtics\n')
pfile.write('set ylabel "E(eV)"\n')
pfile.write('set xlabel "kpoints"\n')
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
pfile.write('set title "Band Structure for cubic %s:%s"\n' % (label0,label1))
pfile.write('plot ')
for i in range(8):
    pfile.write('"harrison.tmp" using %d' % (i+1))
    if i != 7:
        pfile.write(',')
    else:
        pfile.write('\n')
pfile.write('pause -1\n')
if postscript:
    pfile.write('set term post\n')
    pfile.write('set output "harrison.ps"\n')
    pfile.write('replot\n')
if gif:
    pfile.write('set term gif\n')
    pfile.write('set output "harrison.gif"\n')
    pfile.write('replot\n')
pfile.close()
os.system('%s harrison.gnu' % path_to_gnuplot) 
