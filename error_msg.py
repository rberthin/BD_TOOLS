import sys

def error_traj(txt, inp):
    if txt == 'else':
        print("The trajfile specified does not exist.")
        if inp:
            sys.exit("Change the input and come back, bye!")

    elif txt == 'xyz':
        print('The trajectory name should end with xyz..\n')

#******************************************************************************************

def error_file(txt, inp):
    if txt == 'forces':
        print("This forces file does not exist in the current directory.")
        if inp:
            sys.exit("Change the input and come back, bye!")

#******************************************************************************************

def error_func(compute):
    print('The function {} do not exist. Please choose another one\n'.format(compute))
