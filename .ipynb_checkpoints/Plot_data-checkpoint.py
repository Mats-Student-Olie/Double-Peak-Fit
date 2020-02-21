import math
import matplotlib.patches as mpatches
import numpy as np

pi = math.pi
sin = math.sin
cos = math.cos


def Contrast_Factor(C_h00, q, H_sqd):
    C = C_h00*(1-q*H_sqd)
    return C

def H_sqd(h,k,l):
    H_squared = (h**2*k**2+h**2*l**2+k**2*l**2)/(h**2+k**2+l**2)**2
    return H_squared

def K(pk_cntr, wavelength):
    Ka = (2*sin(0.5*pk_cntr*pi/180))/wavelength
    return Ka

def delta_K(pk_cntr, FWHM, wavelength):
    delta_Ka = (FWHM*cos(0.5*pk_cntr*pi/180))/wavelength
    return delta_Ka

def choose_cakes(cakes, z_dK, mx_dK):
    
    merge_cakes = []
    for cake in cakes:
        merge_cakes.append(cake)

    merge_dK_zrh = []
    merge_dK_mx = []

    #only use 220 & first matrix peak (for now)

    for i in merge_cakes:
        dK_zrh_cakes = z_dK[1][i]
        merge_dK_zrh.append(dK_zrh_cakes)
        #merge_dK_zrh[cake][step]

        dK_mx_cakes = mx_dK[0][i]
        merge_dK_mx.append(dK_mx_cakes)

    norm_dK_zrh1 = []
    norm_dK_mx1 = []

    for h, (i, k) in enumerate(zip(merge_dK_zrh, merge_cakes)):
        tmp_zrh = []
        norm_dK_zrh1.append(tmp_zrh)

        for j in i:
            normalised_dK_zrh = j - merge_dK_zrh[h][0]
            norm_dK_zrh1[h].append(normalised_dK_zrh)

    for h, (i, k) in enumerate(zip(merge_dK_mx, merge_cakes)):
        tmp_mx = []
        norm_dK_mx1.append(tmp_mx)

        for j in i:
            normalised_dK_mx = j - merge_dK_mx[h][0]
            norm_dK_mx1[h].append(normalised_dK_mx)

    ## calc average values

    tx_n_dK_zrh = np.transpose(norm_dK_zrh1)
    tx_n_dK_mx = np.transpose(norm_dK_mx1)
    
    av_norm_dK_zrh = np.average(tx_n_dK_zrh, axis = 1)
    av_norm_dK_mx = np.average(tx_n_dK_mx, axis = 1)
    
    return merge_cakes, norm_dK_zrh1, norm_dK_mx1, tx_n_dK_zrh, tx_n_dK_mx, av_norm_dK_zrh, av_norm_dK_mx
    
    
def plt_fwhm_vs_strain(cakes, z, m):

    merge_cakes, norm_dK_zrh1, norm_dK_mx1, tx_n_dK_zrh, tx_n_dK_mx, av_norm_dK_zrh, av_norm_dK_mx = choose_cakes(cakes, z, m)
        
    
    fig, ax0 = plt.subplots(figsize=(18, 10))

    step_range = list(range(58))    
    matrix_colours = ['lightcyan', 'powderblue', 'lightblue','skyblue', 'lightskyblue', 'steelblue']
    hydride_colours = ['papayawhip', 'blanchedalmond', 'moccasin', 'navajowhite', 'wheat', 'burlywood']

    fig.suptitle('Peak Broardening in Hydrided Zr4', fontsize=18, fontweight='bold');

    #matrix

    for i, (data, cake) in enumerate(zip(norm_dK_mx1, merge_cakes)):
        l_lab = 'FWHM (10-10) cake = ' + str(cake)
        ax0.plot(Eng_strain, data, 'o', color = matrix_colours[i], label=l_lab)
    ax0.plot(Eng_strain, medfilt(av_norm_dK_mx, 3), '--', color = 'deepskyblue', label='Av. $\Delta$FWHM Matrix (10-10)')
    ax0_y2 = ax0.twinx()
    for i, cake_val in enumerate(merge_cakes):
        m_lab = 'L_Strain Zr4 (10-10) Cake ' + str(cake_val)
        ax0_y2.plot(Eng_strain, m_strains_by_peak[0][cake_val], color = matrix_colours[i], label = m_lab)

    #hydride
    for i, (data, cake) in enumerate(zip(norm_dK_zrh1, merge_cakes)):
        l_lab = 'FWHM (220) cake = ' + str(cake)
        ax0.plot(Eng_strain, data, 'o', color = hydride_colours[i], label = l_lab)
    ax0.plot(Eng_strain, medfilt(av_norm_dK_zrh, 3), '--', color = 'orange', label = 'Av. $\Delta$FWHM Zrh (220)')
    for i, cake_val in enumerate(merge_cakes):
        s_lab = 'L_Strain ZrH (220) Cake ' + str(cake_val)
        ax0_y2.plot(Eng_strain, z_strains_by_peak[1][cake_val], color = hydride_colours[i], label = s_lab)

    ax0_y2.set_ylabel('microstrain (10$^{-6}$)', fontsize=15)
    ax0.legend(loc = 4,  bbox_to_anchor=(1.3, 0), fontsize=10)
    ax0_y2.legend(loc = 1, bbox_to_anchor=(1.3, 1),   fontsize=10)
    #title0 = ('(100)')
    #ax0.set_title(title0, fontsize=18)
    ax0.set_xlabel('Engineering Strain', fontsize=15)
    ax0.set_ylabel(r'$\Delta$ FWHM (1/nm)', fontsize=15)
    ax0.set_xlim(0.00)
    ax0_y2.set_ylim(-5000, 25000)
    ax0.set_ylim(-0.5, 1.75)
    #ax0.axvline(0.0273, ymin = 0, ymax = 0.5, ls = '--', color = 'dimgrey', lw = 1)
    #ax0.text(0.0273, 0.75, 'Onset of peak\nbroardening\nin Matrix', ha = 'center', wrap = True)
    #ax0.axvline(0.046, ymin = 0, ymax = 0.75, ls = '--', color = 'dimgrey', lw = 1)
    #ax0.text(0.046, 1.25, 'Onset of peak\nbroardening\nin Hydrides', ha = 'center', wrap = True)

    fig.tight_layout(rect=[0, 0, 1, 0.95])

    
def plt_fwhm_vs_step(cakes):
    
    merge_cakes, norm_dK_zrh1, norm_dK_mx1, tx_n_dK_zrh, tx_n_dK_mx, av_norm_dK_zrh, av_norm_dK_mx = choose_cakes(cakes)

    fig, ax0 = plt.subplots(figsize=(18, 10))

    step_range = list(range(58))    
    matrix_colours = ['lightcyan', 'powderblue', 'lightblue','skyblue', 'lightskyblue', 'steelblue']
    hydride_colours = ['papayawhip', 'blanchedalmond', 'moccasin', 'navajowhite', 'wheat', 'burlywood']

    fig.suptitle('Peak Broardening in Hydrided Zr4', fontsize=18, fontweight='bold');

    #matrix

    for i, (data, cake) in enumerate(zip(norm_dK_mx1, merge_cakes)):
        l_lab = 'FWHM (10-10) cake = ' + str(cake)
        ax0.plot(step_range, data, 'o', color = matrix_colours[i], label=l_lab)
    ax0.plot(step_range, medfilt(av_norm_dK_mx, 3), '--', color = 'deepskyblue', label='Av. $\Delta$FWHM Matrix (10-10)')
    ax0_y2 = ax0.twinx()
    for i, cake_val in enumerate(merge_cakes):
        m_lab = 'L_Strain Zr4 (10-10) Cake ' + str(cake_val)
        ax0_y2.plot(step_range, m_strains_by_peak[0][cake_val], color = matrix_colours[i], label = m_lab)

    #hydride
    for i, (data, cake) in enumerate(zip(norm_dK_zrh1, merge_cakes)):
        l_lab = 'FWHM (220) cake = ' + str(cake)
        ax0.plot(step_range, data, 'o', color = hydride_colours[i], label = l_lab)
    ax0.plot(step_range, medfilt(av_norm_dK_zrh, 3), '--', color = 'orange', label = 'Av. $\Delta$FWHM Zrh (220)')
    for i, cake_val in enumerate(merge_cakes):
        s_lab = 'L_Strain ZrH (220) Cake ' + str(cake_val)
        ax0_y2.plot(step_range, z_strains_by_peak[1][cake_val], color = hydride_colours[i], label = s_lab)

    ax0_y2.set_ylabel('microstrain (10$^{-6}$)', fontsize=15)
    ax0.legend(loc = 4,  bbox_to_anchor=(1.3, 0), fontsize=10)
    ax0_y2.legend(loc = 1, bbox_to_anchor=(1.3, 1),   fontsize=10)
    #title0 = ('(100)')
    #ax0.set_title(title0, fontsize=18)
    ax0.set_xlabel('Step', fontsize=15)
    ax0.set_ylabel(r'$\Delta$ FWHM (1/nm)', fontsize=15)
    ax0.set_xlim(0.00)
    ax0.set_ylim(-0.5, 1.75)
    ax0_y2.set_ylim(-5000, 25000)
    #ax0.axvline(23, ymin = 0, ymax = 0.5, ls = '--', color = 'dimgrey', lw = 1)
    #ax0.text(23, 0.8, 'Onset of peak\nbroardening\nin Matrix', ha = 'center', wrap = True)
    #ax0.axvline(30, ymin = 0, ymax = 0.75, ls = '--', color = 'dimgrey', lw = 1)
    #ax0.text(30, 1.26, 'Onset of peak\nbroardening\nin Hydrides', ha = 'center', wrap = True)

    fig.tight_layout(rect=[0, 0, 1, 0.95])

