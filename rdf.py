import numpy as np

def read_one_frame(f, n_atoms):
    """Lit une frame dans le fichier .xyz et renvoie les atomes"""
    f.readline()  # Skip comment line
    data = []
    for _ in range(n_atoms):
        parts = f.readline().split()
        data.append((int(parts[0]), float(parts[1]), float(parts[2]), float(parts[3])))
    return np.array(data)



trajfile = "traj_rho_0.2_eps_2.5.xyz"
box = 20.32
nbins = 500
r_max = box / 2.0
start_step = 5000
end_step = 40000
species_A = 1
species_B = 1

rdf_filename = "rdf_rho_0.2_eps_2.5.dat"


# Initialisation
total_histogram = np.zeros(nbins)
total_steps = 0

dr = r_max / nbins
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

        atoms = read_one_frame(f, n_atoms)

        A = atoms[atoms[:, 0] == species_A][:, 1:4]
        B = atoms[atoms[:, 0] == species_B][:, 1:4]

        nA = len(A)
        nB = len(B)
        histogram = np.zeros(nbins)

        if nA == 0 or nB == 0:
            continue

        density_B = nB / box**3 if species_A != species_B else nA / box**3

        if species_A == species_B:
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
        if species_A == species_B:
            counts *= 2  # chaque paire est unique

        # Normalisation par le volume de la coquille sphérique
        r_outer = (np.arange(1, nbins + 1)) * dr
        r_inner = (np.arange(0, nbins)) * dr
        shell_volumes = (4.0 / 3.0) * np.pi * (r_outer**3 - r_inner**3)
        with np.errstate(divide='ignore', invalid='ignore'):
            histogram = counts / (nA * density_B * shell_volumes)

        total_histogram += histogram
        total_steps += 1

# Écriture des résultats
with open(rdf_filename, "w") as out:
    r = (np.arange(nbins) + 0.5) * dr
    g_r = total_histogram / total_steps
    for radius, value in zip(r, g_r):
        out.write(f"{radius:.5f} {value:.5f}\n")
