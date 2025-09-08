import numpy as np
import matplotlib.pyplot as plt

def autoforce_computation(forcefile, freq, delta_t, n_atoms):
    fx, fy, fz = np.loadtxt(forcefile, unpack=True)

    nsteps = fx.shape[0] // n_atoms

    # Reshape des forces : (nsteps, natoms)
    fxx = fx.reshape(nsteps, n_atoms)
    fyy = fy.reshape(nsteps, n_atoms)
    fzz = fz.reshape(nsteps, n_atoms)

    # Stack des composantes en un seul tableau : (nsteps, 100, 3)
    forces = np.stack((fxx, fyy, fzz), axis=2)

    # Initialisation de la matrice de corrélation
    Z = np.zeros(nsteps)

    # Calcul de l'autocorrélation pour chaque particule (vectorisé sur k)
    for t in range(nsteps):
        print('Step', t)
        # vecteurs force à t=0 et t décalé
        a1 = forces[:nsteps - t]
        a2 = forces[t:]

        # Produit scalaire entre les vecteurs (axis=2)
        dots = np.sum(a1 * a2, axis=2)  # (nsteps - t, 100)

        # Moyenne sur les temps et les particules
        Z[t] = np.mean(dots)

    X = np.arange(nsteps) * freq * delta_t
    return X, Z

def write_autoforce(X, Z):
    fileout = 'autocorr_forces.out'
    np.savetxt(fileout, np.c_[X, Z])

def plot_autoforce(X, Z):
    fig = plt.figure(figsize=(6,6), tight_layout = True)
    name_file_save=r'autocorr_forces'
    name_file_save+='.png'

    plt.plot(X, Z)
    plt.xlabel(r't', fontsize = 14)
    plt.ylabel(r'< F(0) . F(t) >', fontsize = 14)

    fig.savefig(name_file_save, bbox_inches='tight')
    print('Plot saved in autocorr_forces.png ! ')

