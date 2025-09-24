#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
Script        : bd_tool.py
Auteur        : Roxanne BERTHIN
Créé le       : 30/07/25
Description   : Ce code permet le calcul de plusieurs fonctions utiles dans le 
                cas de l'analyse de trajectoire de dynamique brownienne ou de 
                dynamique moleculaire

Utilisation   : 
    python nom_du_script.py [options]

Arguments     :
    -i, --input     : Chemin du fichier d’input. Si aucun fichier d'input alors 
                      le script posera des questions
    -p, --plot      : permet de tracer le résultat 
    -h, --help      : Affiche ce message d’aide

Dépendances   :
    - numpy
    - matplotlib (facultatif)

Notes         : RAS

Licence       : RAS 

===============================================================================
"""

# librairie à importer
import argparse
from rdf import rdf_computation, write_rdf, plot_rdf
from force_autocorrelation import autoforce_computation, write_autoforce, plot_autoforce
from msd import msd_computation, write_msd, plot_msd
from bd_stuff import unwrap_xyz_traj, write_unwrap
import error_msg
import os

# Welcome menu
print('-------------------------------------------------------')
print('*******************************************************')
print('------------------- !! BD  TOOLS !! -------------------')
print('*******************************************************')
print('-------------------------------------------------------')

# Arguments 
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', type=str, 
                    help = 'input file containing information about the system')

parser.add_argument('-p', '--plot', action = 'store_true', 
                                help = 'save a plot in png format')
args = parser.parse_args()

bool_traj = True
bool_func = True
func_list = ['rdf', 'autoforce', 'unwrap', 'msd']

if args.input:
    inputfile = open(args.input, 'r')
    compute = inputfile.readline().rstrip()
    print('Reading file {}\n'.format(args.input))
else:
    print('No input file specified ...\n')

    print('-------------------------------------------------------')
    print('--------------------- COMPUTATION ---------------------')
    print('-------------------------------------------------------')
    print('List of functions:')
    print('------------------')
    print('-- Radial distribution function (rdf)')
    print('-- Force autocorrelation (autoforce)')
    print('-- Unwrap trajectory (unwrap)')
    print('-- Mean square displacement (msd) a venir')
    print('-- Droplet analysis (droplet) a venir')
    print('-------------------------------------------------------')
    while bool_func:
        compute = input('Which function do you want to compute?\n')
        if compute not in func_list:
            error_msg.error_func(compute)
        else:
            bool_func = False

match compute:
    case "rdf":
        if args.input:
            trajfile = inputfile.readline().rstrip()
            if not os.path.exists(trajfile):
                error_msg.error_traj('else', args.input)
            else:
                if not trajfile.endswith('.xyz'):
                    error_msg.error_traj('xyz', args.input)

            box = float(inputfile.readline())

            nbins = int(inputfile.readline())
            start_step = int(inputfile.readline())
            end_step = int(inputfile.readline())
            species_1 = inputfile.readline().rstrip()
            species_2 = inputfile.readline().rstrip()
            inputfile.close()
        else:
            while bool_traj:
                trajfile = input('Name of the trajectory file?\n')
                if not os.path.exists(trajfile):
                    error_msg.error_traj('else', args.input)

                else:
                    if not trajfile.endswith('.xyz'):
                        error_msg.error_traj('xyz', args.input)
                    else:
                        bool_traj = False

            box = float(input('Length of the box?\n'))
            nbins = int(input('Number of bins?\n'))
            start_step = int(input('Starting step for the computation of rdf?\n'))
            end_step = int(input('Ending step for the computation of rdf?\n'))
            species_1 = input('Reference specie?\n')
            species_2 = input('Observed specie?\n')
            print('-------------------------------------------------------')

        # call rdf function
        r, g_r = rdf_computation(trajfile, box, nbins, start_step,
                                 end_step, species_1, species_2)
        write_rdf(r, g_r)

        if args.plot:
            plot_rdf(r, g_r)

    case "autoforce":
        bool_forces = True
        if args.input:
            forcefile = inputfile.readline().rstrip() 
            if not os.path.exists(forcefile):
                error_msg.error_file('forces', args.input)
            freq = int(inputfile.readline())
            delta_t = float(inputfile.readline())
            n_atoms = int(inputfile.readline())

        else:
            while bool_forces:
                forcefile = input("Name of the file containing fx, fy, fz?\n")
                if not os.path.exists(forcefile):
                    error_msg.error_file('forces', args.input)
                else:
                    bool_forces = False

            freq = int(input("What is the recording frequency (in Dt unit)?\n"))
            delta_t = float(input("What is the timestep?\n"))
            n_atoms = int(input("Number of atoms?\n"))
            print('-------------------------------------------------------')

        X, Z = autoforce_computation(forcefile, freq, delta_t, n_atoms)
        
        write_autoforce(X, Z)

        if args.plot:
            plot_autoforce(X, Z)

    case "unwrap":
        if args.input:
            trajfile = inputfile.readline().rstrip()
            if not os.path.exists(trajfile):
                error_msg.error_traj('else', args.input)
            else:
                if not trajfile.endswith('.xyz'):
                    error_msg.error_traj('xyz', args.input)

            box = float(inputfile.readline())
            n_atoms = int(inputfile.readline())
            start_step = int(inputfile.readline())
            end_step = int(inputfile.readline())

        else:
            while bool_traj:
                trajfile = input('Name of the trajectory file?\n')
                if not os.path.exists(trajfile):
                    error_msg.error_traj('else', args.input)

                else:
                    if not trajfile.endswith('.xyz'):
                        error_msg.error_traj('xyz', args.input)
                    else:
                        bool_traj = False

            box = float(input('Length of the box?\n'))
            n_atoms = int(input('Number of atoms?\n'))
            start_step = int(input('Starting step for the unwrap of the trajectory?\n'))
            end_step = int(input('Ending step for the unwrap of the trajectory?\n'))
            print('-------------------------------------------------------')

        label, xx, yy, zz = unwrap_xyz_traj(trajfile, box, n_atoms, start_step, end_step)
        write_unwrap(label, xx, yy, zz)

    case "msd":
        if args.input:
            trajfile = inputfile.readline().rstrip()
            if not os.path.exists(trajfile):
                error_msg.error_traj('else', args.input)
            else:
                if not trajfile.endswith('.xyz'):
                    error_msg.error_traj('xyz', args.input)

            box = float(inputfile.readline())        
            freq = int(inputfile.readline())
            delta_t = float(inputfile.readline())
            n_atoms = int(inputfile.readline())
            start_step = int(inputfile.readline())
            end_step = int(inputfile.readline())
            selected_species = inputfile.readline()     
        else:
            while bool_traj:
                trajfile = input('Name of the trajectory file?\n')
                if not os.path.exists(trajfile):
                    error_msg.error_traj('else', args.input)

                else:
                    if not trajfile.endswith('.xyz'):
                        error_msg.error_traj('xyz', args.input)
                    else:
                        bool_traj = False

            box = float(input('Length of the box?\n'))
            freq = int(input("What is the recording frequency (in Dt unit)?\n"))
            delta_t = float(input("What is the timestep?\n"))
            n_atoms = int(input('Number of atoms?\n'))
            start_step = int(input('Starting step for the unwrap of the trajectory?\n'))
            end_step = int(input('Ending step for the unwrap of the trajectory?\n'))
            selected_species = input('On which species do you want to compute msd?\n')
            print('-------------------------------------------------------')

        print('-- Unwraping trajectory ... ---------------------------') 
        label, xx, yy, zz = unwrap_xyz_traj(trajfile, box, n_atoms, start_step, end_step) 
        print('-- Trajectory unwraped ! ------------------------------\n')
        print('-- Starting the msd computation ... -------------------')
        time, msd = msd_computation(label, xx, yy, zz, start_step, end_step, selected_species)
        write_msd(time, msd)
