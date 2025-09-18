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
def read_all_xyz_traj(trajfile, n_atoms, end_step):
    label = np.zeros( (n_atoms, end_step), dtype=str)
    coord_x = np.zeros( (n_atoms, end_step) )
    coord_y = np.zeros( (n_atoms, end_step) )
    coord_z = np.zeros( (n_atoms, end_step) )
    
    with open(trajfile, "r") as f:
        for i in range(end_step):
            f.readline()
            f.readline()
            for j in range(n_atoms):

                parts = f.readline().split()
                label[j, i] = str(parts[0])
                coord_x[j, i] = parts[1]
                coord_y[j, i] = parts[2]
                coord_z[j, i] = parts[3]
    
    return np.array(label), np.array(coord_x, dtype=float), np.array(coord_y, dtype=float), np.array(coord_z, dtype=float)

#*************************************************************
def minimum_image_convention(pos, box):
    if pos > 0.5*box:
        pos -= box
    elif pos < -0.5*box:
        pos += box
    else:
        pos = pos
    return pos

#*************************************************************
def unwrap_xyz_traj(trajfile, box, n_atoms, start_step, end_step):
    label, coord_x, coord_y, coord_z = read_all_xyz_traj(trajfile, n_atoms, end_step)
    shape = end_step - start_step
    
    typ = label[:, start_step:end_step]
    pos_x = coord_x[:, start_step:end_step]
    pos_y = coord_y[:, start_step:end_step]
    pos_z = coord_z[:, start_step:end_step]
    
    xx = np.zeros( (n_atoms, shape) )
    yy = np.zeros( (n_atoms, shape) )
    zz = np.zeros( (n_atoms, shape) )
    
    for i in range(n_atoms):
        xx[i, 0] = pos_x[i, 0]
        yy[i, 0] = pos_y[i, 0]
        zz[i, 0] = pos_z[i, 0]
        for j in range(1, end_step):
            dx = minimum_image_convention(pos_x[i, j] - pos_x[i, j-1], box)
            dy = minimum_image_convention(pos_y[i, j] - pos_y[i, j-1], box)
            dz = minimum_image_convention(pos_z[i, j] - pos_z[i, j-1], box)
            
            xx[i,j] = xx[i, j-1] + dx
            yy[i,j] = yy[i, j-1] + dy
            zz[i,j] = zz[i, j-1] + dz
    return typ, xx, yy, zz

#*************************************************************
def write_unwrap(typ, xx, yy, zz):
    f_out = open('unwrap_traj.xyz', 'w')
    for j in range(np.shape(xx)[1]):
        f_out.write('{0}\n'.format(np.shape(xx)[0]))
        f_out.write('Step {0}\n'.format(j))
        for i in range(np.shape(xx)[0]):
            f_out.write('{}   {:.5f}   {:.5f}   {:.5f}\n'.format(
                typ[i,j], xx[i,j], yy[i,j], zz[i, j]))

#*************************************************************

