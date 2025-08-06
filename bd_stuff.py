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

#*************************************************************
def read_all_xyz_traj(f, n_atoms, end_step):
    for _ in range(end_step):
        f.readline()
        f.readline()
        for _ in range(n_atoms):
            parts = f.readline().split()
            label.append(parts[0])
            coords.append([float(x) for x in parts[1:4]])
    # il faut reshape !! 
    return np.array(label), np.array(coords, dtype=float)

#*************************************************************
def minimum_image_convention(pos, box):
    if pos > 0.5*box:
        pos -= box
    elif pos < 0.5*box:
        pos += box
    else:
        pos = pos
    return pos

#*************************************************************
def unwrap_xyz_traj(f, n_atoms, start_step, end_step):
    label, coord_x, coord_y, coord_z = read_all_xyz_traj(f, n_atoms, end_step)
    shape = end_step - start_step
    
    xx = np.array( (n_atoms, shape) )
    yy = np.array( (n_atoms, shape) )
    zz = np.array( (n_atoms, shape) )

    for i in range(n_atoms):
        xx = coord_x[i, start_step]
        yy = coord_y[i, start_step]
        zz = coord_z[i, start_step]
        for j in range(start_step+1, end_step):
            dx = minimum_image_convention(coord_x[i, j] - coord_x[i, j-1], box)
            dy = minimum_image_convention(coord_y[i, j] - coord_y[i, j-1], box)
            dz = minimum_image_convention(coord_z[i, j] - coord_z[i, j-1], box)
            
            xx[i,j] = xx[i, j-1] + dx
            yy[i,j] = yy[i, j-1] + dy
            zz[i,j] = zz[i, j-1] + dz
    return label, xx, yy, zz

#*************************************************************

