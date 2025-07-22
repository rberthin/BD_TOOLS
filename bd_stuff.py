import numpy as np

def read_one_xyz_frame(f, n_atoms):
    """Lit une frame dans le fichier .xyz et renvoie les atomes"""
    f.readline()  # Skip comment line
    data = []
    for _ in range(n_atoms):
        parts = f.readline().split()
        data.append((int(parts[0]), float(parts[1]), float(parts[2]), float(parts[3])))
    return np.array(data)

