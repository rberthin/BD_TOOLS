import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

def msd_computation(label, xx, yy, zz, delta_t, freq, start_step, end_step, selected_species):
    name_type, count_type = np.unique(label[:,0], return_counts=True)
    index = np.where(name_type == selected_species.strip())
    nselected = count_type[index][0]
    
    size_sample = end_step - start_step
    x_type = np.zeros((nselected, size_sample))
    y_type = np.zeros((nselected, size_sample))
    z_type = np.zeros((nselected, size_sample))
    
    for i in range(np.shape(label)[1]): # loop on step
        count = 0
        for j in range(np.shape(label)[0]): # loop on atoms
            if str(label[j,i]) == selected_species.strip():
                x_type[count, i] = xx[j,i]
                y_type[count, i] = yy[j,i]
                z_type[count, i] = zz[j,i]
                count += 1
    
    coords = np.stack([x_type, y_type, z_type], axis=-1)
    
    msd = np.zeros( (size_sample) )

    for k in range(1, size_sample):
        print('Step', k, '/', size_sample)
        disp = coords[:, k:] - coords[:, :-k]   # shape (nselected, size_sample-k, 3)
        sqdist = np.sum(disp**2, axis=-1)       # shape (nselected, size_sample-k)
        msd[k] = np.mean(sqdist)                # moyenne directe

    time = np.arange(size_sample)
    diff, intercept, r, p, se = stats.linregress(time[1:int(size_sample/2)], msd[1:int(size_sample/2)])
    
    print('Diffusion coefficient estimated : {}\n'.format(round(diff/(delta_t*freq), 3)))
    
    return time, msd

def write_msd(time, msd):
    fileout = 'msd.out'
    np.savetxt(fileout, np.c_[time[1:], msd[1:]])

def plot_msd(time, msd):
    fig = plt.figure(figsize=(6,6), tight_layout = True)
    name_file_save=r'msd'
    name_file_save+='.png'

    plt.plot(time, msd)
    plt.xlabel(r't', fontsize = 14)
    plt.ylabel(r'MSD', fontsize = 14)

    fig.savefig(name_file_save, bbox_inches='tight')
    print('Plot saved in msd.png ! ')


