#** Libraries
#-------------------------------------
import numpy as np
import matplotlib.pyplot as plt
from bd_stuff import read_one_xyz_frame


def rdf_computation(trajfile, box, nbins, start_step,
                    end_step, species_1, species_2):

    r_max = box / 2.0
    dr = r_max / nbins
    
    total_histogram = np.zeros(nbins)
    total_steps = 0
    
    with open(trajfile, "r") as f:
        current_step = 0
    
        while current_step < end_step:
            print('Step', current_step)
            line = f.readline()
            if not line:
                break
            n_atoms = int(line.strip())
            current_step += 1
    
            if current_step < start_step:
                for _ in range(n_atoms + 1):
                    f.readline()
                continue

            label, coords = read_one_xyz_frame(f, n_atoms)
            label = np.char.strip(label)

            A = coords[label[:] == species_1]#[:, 1:4]
            B = coords[label[:] == species_2]#[:, 1:4]

            nA = len(A)
            nB = len(B)

            if nA == 0 or nB == 0:
                continue

            histogram = np.zeros(nbins)

            density_B = nB / box**3 if species_1 != species_2 else nA / box**3

            if species_1 == species_2:
                # Calcul des paires uniques (i < j) pour éviter les doublons
                indices = np.triu_indices(nA, k=1)
                diffs = A[indices[0]] - A[indices[1]]
            else:
                diffs = A[:, None, :] - B[None, :, :]
                diffs = diffs.reshape(-1, 3)

            # Conditions périodiques (minimum image convention)
            diffs -= box * np.round(diffs / box)
            dists = np.linalg.norm(diffs, axis=1)

            # Filtrage des distances < rmax
            valid = dists < r_max
            dists = dists[valid]
            bin_indices = (dists / dr).astype(int)

            counts = np.bincount(bin_indices, minlength=nbins)
            if species_1 == species_2:
                counts *= 2  # chaque paire est unique

            # Normalisation par le volume de la coquille sphérique
            r_outer = (np.arange(1, nbins + 1)) * dr
            r_inner = (np.arange(0, nbins)) * dr
            shell_volumes = (4.0 / 3.0) * np.pi * (r_outer**3 - r_inner**3)
            with np.errstate(divide='ignore', invalid='ignore'):
                histogram = counts / (nA * density_B * shell_volumes)

            total_histogram += histogram
            total_steps += 1

    r = (np.arange(nbins) + 0.5) * dr
    g_r = total_histogram / total_steps
    return r, g_r

def write_rdf(r, g_r):
    # Écriture des résultats
    rdf_filename = "rdf.dat"
    with open(rdf_filename, "w") as out:
        for radius, value in zip(r, g_r):
            out.write(f"{radius:.5f} {value:.5f}\n")

def plot_rdf(r, g_r):
    fig = plt.figure(figsize=(6,6), tight_layout = True)
    name_file_save=r'rdf'
    name_file_save+='.png'

    plt.plot(r, g_r)
 
    plt.xlabel(r'r', fontsize = 14)
    plt.ylabel(r'g(r)', fontsize = 14)

    fig.savefig(name_file_save, bbox_inches='tight')

