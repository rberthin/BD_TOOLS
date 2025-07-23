import numpy as np

def read_one_xyz_frame(f, n_atoms):
    """Lit une frame dans le fichier .xyz et renvoie les atomes"""
    f.readline()  # Skip comment line
    label = []
    coords = []
    for _ in range(n_atoms):
        parts = f.readline().split()
        label.append(parts[0])
        coords.append([float(x) for x in parts[1:4]])
    return np.array(label), np.array(coords, dtype=float)