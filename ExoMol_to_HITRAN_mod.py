# -*- coding: utf-8 -*-
from __future__ import print_function
# ExoMol_to_HITRAN.py
# Christian Hill (christian.hill@ucl.ac.uk)
# v2.0, 2/12/15
# v1.0, 4/7/14

# ExoMol_to_HITRAN is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ExoMol_to_HITRAN is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with ExoMol_to_HITRAN.  If not, see <http://www.gnu.org/licenses/>

"""
A script to convert an ExoMol dataset into (160-character)
HITRAN format.

Run this program as:
python ExoMol_to_HITRAN.py <isotopologue> <dataset name> <numin> <numax> <T>
for example,
python ExoMol_to_HITRAN.py 14N-1H3 BYTe 1000. 3000. 296

Pass -e or --enforce-ordering to require input transitions to be ordered by
wavenumber.

"""
import pdb as pdb
import sys, os
import bz2
import numpy as np
from glob import glob
import argparse
try:
    from exomol_def import ExoMolDef
except ImportError:
    pass

parser = argparse.ArgumentParser(description='Convert from the'
            ' ExoMol data format to a format compatible with the'
            ' HITRAN2004+ format.')
parser.add_argument('--def', metavar='<def_file>', default=None,
        help='the isotopologue ExoMol def file', action='store',
        dest='def_file')
parser.add_argument('iso', metavar='<isotopologue>', nargs='?',
        help='the isotopologue name as e.g. 1H2-16O', default=None)
parser.add_argument('dataset', metavar='<data-set>', default=None, nargs='?',
       help='the ExoMol data set to take the data from, e.g. BT2')
parser.add_argument('numin', metavar='<numin>', type=float,
        help='the minimum wavenumber to output (in cm-1)')
parser.add_argument('numax', metavar='<numax>', type=float,
        help='the maximum wavenumber to output (in cm-1)')
parser.add_argument('T', metavar='<T>', type=float,
        help='the temperature (in K) at which to calculate the'
             ' line intensity, S')
parser.add_argument('-d', '--data_dir', dest='data_dir',
        action='store', default='.',
        help='path to the directory holding the ExoMol data'
             ' files for this isotopologue')
parser.add_argument('-s', '--smin', dest='Smin', action='store',
                    type=float, default=0., help='Minimum line'
             ' intensity, S, to output')
parser.add_argument('-a', '--abundance', dest='abundance', action='store',
                    type=float, default=1., help='Isotopologue abundance'
                    'scaling factor')
parser.add_argument('-Q', '--pf', dest='Q', action='store',
                    type=float, default=None, help='Partition sum')
parser.add_argument('-e', '--enforce-ordering', dest='enforce_ordering',
                    action='store_true', default=False,
                    help='Require input transitions to be ordered by'
                         ' increasing wavenumber')

args = parser.parse_args()
iso, dataset, numin, numax, T, data_dir, Smin, abundance = (args.iso,
    args.dataset, args.numin, args.numax, args.T, args.data_dir,
    args.Smin, args.abundance)    #Set variables equal to those given as arguments
def_file = args.def_file          #Set def_file

def read_def_file(def_filename):
    if def_filename is None:
        def_filename = '{0}-{1}.def'.format(iso, dataset)
    def_filepath = os.path.join(data_dir, def_filename)
    edef = ExoMolDef(def_filepath)
    edef.read_def()
    return edef

if def_file is None:
    if iso is None or dataset is None:
        print('No def file specified, but no iso or dataset provided either.')
        sys.exit(1)
else:
    if iso is not None or dataset is not None:
        print('A def file was specified, but iso or dataset provided too.')
        sys.exit(1)
    edef = read_def_file(def_file)
    iso = edef.iso_slug
    dataset = edef.dataset

def read_states(states_name):
    if not os.path.exists(states_name):
        states_name = states_name + '.bz2'
    print('reading states from', states_name, '... ', end='')
    states = np.loadtxt(states_name, usecols=(1,2,10),    #E_state, gtot, v
                        dtype=[('E', 'f8'), ('g', 'u4'), ('v','f2')])
    print('Done')
    return states

def get_trans_names(numin, numax):
    if len(all_trans_names) == 1:
        # all the transitions are in one file
        return all_trans_names
    trans_names = []
    for trans_name in all_trans_names:
        # strip all paths and extensions
        ext = None
        filestem = os.path.basename(trans_name)
        filestem, ext = os.path.splitext(filestem)
        if ext == '.bz2':
            filestem, ext = os.path.splitext(filestem)

        fields = filestem.split('__')
        s_nurange = fields[2].split('-')
        fnumin, fnumax = float(s_nurange[0]), float(s_nurange[1])
        if fnumax < numin:
            continue
        if fnumin > numax:
            break
        trans_names.append(trans_name)
    print('The matching transitions files found:', trans_names)
    return trans_names

def read_trans(trans_name):
    if trans_name.endswith('bz2'):
        fi = bz2.BZ2File(trans_name)
    else:
        fi = open(trans_name)
    for line in fi:
        fields = line.split()
        # get the state energies and degeneracies from their stateIDs
        stateIDp = int(fields[0])
        stateIDpp = int(fields[1])
        # NB ExoMol counts from 1, but Python counts from 0
        Ep, gp, vp = states[stateIDp-1]     #Access upper level energy, degeneracy
        Epp, gpp, vpp = states[stateIDpp-1]  #Access lower level energy, degeneracy
        # Einstein A coefficient
        A = float(fields[2])
        # calculate transition wavenumber
        nu = Ep - Epp
        yield (nu, Epp, gp, gpp, vp, vpp, A)

def read_Q(q_name):
    _T = []; _Q = []
    with open(q_name, 'r') as f:
        for line in f.readlines():
            fields = line.split()
            _T.append(float(fields[0]))
            _Q.append(float(fields[1]))
    return _T, _Q
            
def use_provided_Q():
    if args.Q is None:
        print('No .pf file found and no partition sum provided.')
        sys.exit(1)
    return args.Q

def get_Q(q_name, T):
    try:
        _T, _Q = read_Q(q_name)
    except FileNotFoundError:
        return use_provided_Q()

    if T < _T[0] or T > _T[-1]:
            # requested temperature out of range
            return None

    # horribly inefficient way to bracket the temperature requested
    # but the partition function list is only small, so who cares?
    for i in range(len(_T)):
        if _T[i] >= T:
            Thi = _T[i]
            Tlo = _T[i-1]
            break
    # linear interpolation
    return _Q[i-1] + (_Q[i] - _Q[i-1])/(Thi - Tlo) * (T - Tlo)

def calc_S(nu0, Epp, gp, A, Q, T):
    c2oT = c2 / T
    S = A * gp * np.exp(-c2oT * Epp) * (1 - np.exp(-c2oT * nu0))\
             / pic8 / nu0**2 / Q
    return S

#http://cedadocs.ceda.ac.uk/998/1/hitran_format.html#:~:text=The%20HITRAN%20format%20is%20the,in%20formatted%20ASCII%20text%20files.
def get_parline(nu, S, A, Epp, gp, gpp, Vp, Vpp):
    """
    Return a 160-character, HITRAN2004+ formatted line containing the
    wavenumber, line intensity, Einstein A coefficient, lower-state energy,
    and upper and lower state degeneracies, but blank space elsewhere.

    """
    Qp = Qpp = ' '*15
    return '601%12.6f%10.3e%10.3e.00000.000%10.4f'\
           '0.00-.000000%15s%15s%15s%15s%19s%7.1f%7.1f'\
         % (nu, S, A, Epp, Vp, Vpp, Qp, Qpp, ' '*19,
            float(gp), float(gpp))
    #mol number of 60 is just a placeholder for SiO for now, since HITRAN goes up to 55

#Define file name(s) of trans files
exomol_base = '%s__%s' % (iso, dataset)
print('exomol_base: ',exomol_base)   #base name for files, e.g. 28Si-16O__SiOUVenIR
states_name = os.path.join(data_dir, '%s.states' % exomol_base)
print('states_name: ',states_name)   #Name of staes file, e.g. 28Si-16O__SiOUVenIR.states
trans_names_glob = os.path.join(data_dir, '%s*.trans' % exomol_base)
all_trans_names = glob(trans_names_glob)
bz2_trans_names_glob = os.path.join(data_dir, '%s*.trans.bz2' % exomol_base)
bz2_trans_names = glob(bz2_trans_names_glob)
for bz2_trans_name in bz2_trans_names:
    trans_name = os.path.splitext(bz2_trans_name)[0]
    if trans_name not in all_trans_names:
        all_trans_names.append(bz2_trans_name)
all_trans_names.sort()
trans_names = get_trans_names(numin, numax)

#Define partition function for given temperature
q_name = os.path.join(data_dir, '%s.pf' % exomol_base) #Define file name of partition function file
Q = get_Q(q_name, T)  #Determine Q at given T 
print('Q(%d K) = %f' % (T, Q))  #Print to command line

#Define some constants
# Planck's const (J s); speed of light (m s-1); Boltzmann's const (J K-1)
h = 6.62607015e-34; c = 299792458; kB = 1.380649e-23
c2 = h * c * 100 / kB     # second radiation constant, cm.K
pic8 = 8 * np.pi * c * 100    # 8.pi.c in cm-1 s

#Read in states file - only reads in E_state, gtot and v, as three columns of "states" variable
states = read_states(states_name)

#Define and open output file
out_name = '%s__%d-%d__%dK.par' % (exomol_base, int(numin), int(numax), int(T))
fo = open(out_name, 'w')
# number of lines calculated:
ncalc = 0
# number of lines written:
nlines = 0
print('writing HITRAN .par file for %d - %d cm-1 at T = %d K...'
            % (numin, numax, T))
if Smin:
    print('Intensity threshold, Smin =', Smin, 'cm-1/(molec.cm-2)')
last_nu = -1000.

for trans_name in trans_names:  #loop through trans files
    print('reading transitions from', trans_name)
#read_trans reads in transition, also identifies states, and returns nu, Elow, gup, glow, A
    for (nu, Epp, gp, gpp, vp, vpp, A) in read_trans(trans_name):  
        if nu < numin:
            continue

        # If we insist that the incoming transitions are ordered by wavenumber
        # we can break as soon as nu > numax.
        if args.enforce_ordering and nu > numax:
            break
        # If the incoming transitions may not be strictly ordered, wait until
        # we're some way above numax before breaking off.
        if nu > numax + 10:
            break
        if nu > numax:
            continue

        if args.enforce_ordering and nu - last_nu < -1.e-2:
            print('transitions not ordered by wavenumber.')
            print(nu, Epp, gp, gpp, A)
            print(last_nu)
            sys.exit(1)
        last_nu = nu

        S = calc_S(nu, Epp, gp, A, Q, T)  #Calculate line strength
        S *= abundance
        ncalc += 1
        if not ncalc % 100000:  #modulo (remainder) function - whenever ncalc is 100000, print current nu
            print('At %.6f cm-1' % nu)
        if S < Smin:
            continue
#        if((vp<5) & (vpp<5) & (vp>vpp)):   #Hardcoded V limit - should probably consider changing it to an arg
        if(Epp<10000):   #Hardcoded E limit - should probably consider changing it to an arg
            print(get_parline(nu, S, A, Epp, gp, gpp, vp, vpp), file=fo)   #print to par file
        nlines += 1
fo.close()

print(nlines, 'lines written to', out_name)
