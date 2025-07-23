# librairie à importer
import argparse
from rdf import rdf_computation, write_rdf, plot_rdf
import error_msg

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
func_list = ['rdf']

if args.input:
    inputfile = open(args.input, 'r')
    trajfile = inputfile.readline().rstrip()
    box = float(inputfile.readline())
    compute = inputfile.readline().rstrip()
else:
    print('No input file specified ...\n')
    print('-------------------------------------------------------')
    print('----------------- GENERAL INFORMATION -----------------')
    print('-------------------------------------------------------')

    while bool_traj:
        trajfile = input('Name of the trajectory file?\n')
        if not trajfile.endswith('.xyz'):
            error_msg.error_traj()
        else:
            bool_traj = False
    
    box = float(input('Length of the box?\n'))
    print('-------------------------------------------------------')
    print('--------------------- COMPUTATION ---------------------')
    print('-------------------------------------------------------')
    print('List of functions:')
    print('------------------')
    print('-- Radial distribution function (rdf)')
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
            nbins = int(inputfile.readline())
            start_step = int(inputfile.readline())
            end_step = int(inputfile.readline())
            species_1 = inputfile.readline().rstrip()
            species_2 = inputfile.readline().rstrip()
            inputfile.close()
        else:
            nbins = int(input('Number of bins?\n'))
            start_step = int(input('Starting step for the computation of rdf?\n'))
            end_step = int(input('Ending step for the computation of rdf?\n'))
            species_1 = input('Reference specie?\n')
            species_2 = input('Observed specie?\n')
            print('-------------------------------------------')

        # call rdf function
        r, g_r = rdf_computation(trajfile, box, nbins, start_step,
                                 end_step, species_1, species_2)
        write_rdf(r, g_r)

        if args.plot:
            plot_rdf(r, g_r)
