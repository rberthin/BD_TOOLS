import numpy as np
import matplotlib.pyplot as plt
from collections import Counter

def msd_computation(label, xx, yy, zz, start_step, end_step, selected_species):
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
    
    msd = np.zeros( (size_sample) )
    sum_r2 = np.zeros( (size_sample) )

    for k in range(1, size_sample):
        count = 0
        for i in range(size_sample - k):
            sum_r2[:] = 0.0 
            for j in range(nselected):
                dx = x_type[j,i+k] - x_type[j,i]
                dy = y_type[j,i+k] - y_type[j,i]
                dz = z_type[j,i+k] - z_type[j,i]
                sum_r2[i] = sum_r2[i] + (dx**2 + dy**2 + dz**2)
            msd[k] = msd[k] + sum_r2[i]
            count = count + nselected
        if count > 0:
            msd[k] = msd[k] / count
        else:
            msd[k] = 0
    time = np.arange(size_sample)
    return time, msd

def write_msd(time, msd):
    fileout = 'msd.out'
    np.savetxt(fileout, np.c_[time, msd])

def plot_msd(X, Z):
    fig = plt.figure(figsize=(6,6), tight_layout = True)
    name_file_save=r'autocorr_forces'
    name_file_save+='.png'

    plt.plot(X, Z)
    plt.xlabel(r't', fontsize = 14)
    plt.ylabel(r'< F(0) . F(t) >', fontsize = 14)

    fig.savefig(name_file_save, bbox_inches='tight')
    print('Plot saved in autocorr_forces.png ! ')


