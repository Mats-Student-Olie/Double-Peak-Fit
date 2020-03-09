import math
import matplotlib.patches as mpatches
import numpy as np
import matplotlib.pyplot as plt
from ipywidgets import FloatProgress
from IPython.display import display
from scipy.signal import medfilt
from ipywidgets import HBox, FloatSlider
import glob, os
import ipywidgets as widgets
import time
import os
import pickle

pi = math.pi
sin = math.sin
cos = math.cos


def unpack_pickles(dic_names, pickles_dir):    
    
    list_matrix_peaks = []
    list_matrix_peak_centers = []
    list_matrix_peak_fwhm = []
    list_matrix_peak_amplitude = []

    list_ZrH_peaks = []
    list_ZrH_peak_centres = []
    list_ZrH_peak_fwhm = []
    list_ZrH_peak_amplitude = []

    for pickles in dic_names:
        with open(pickles_dir + '\\' + pickles, 'rb') as x:
            data = pickle.load(x)
        #print(len(data))
        #enter names of reflections in dictionary you intend to analyse
        reflections=pickles[pickles.find('('):-7]
        cake = pickles.split('_')[1]

        # Seperate reflections from a filename
        def find(ch,string1):
            pos = []
            for i in range(len(string1)):
                if ch == string1[i]:
                    pos.append(i)
            return(pos)

        allpeaks = []
        for i,j in zip(find('(', pickles), find(')', pickles)):
            allpeaks.append(pickles[i:j+1])

        numpeaks=len(allpeaks)
        #print("Reflections: {: >25}      Cake = {: >3}       Peaks in Pickle = {: >3}      Peaks List: {}".format(reflections, cake, numpeaks, allpeaks))

    ### Extracting peak centres

        firstFile = 0
        lastFile = len(data)
        #print('Number of Strain Steps: {0}'.format(lastFile))
        step = 1

    ### isolate hydride reflections

        #n = number of reflections in peak_list, peak will be the reflections in the peak_list
        for n, peak in enumerate(allpeaks):

            #set up dictionaries to store the peak information
            m_peak_centres={}
            m_peak_fwhm={}
            m_peak_amplitude={}    
            ZrH_peak_centres={}
            ZrH_peak_fwhm={}
            ZrH_peak_amplitude={}

            if len(peak) > 5:
                m_peak_centres[peak]=[]
                m_peak_fwhm[peak]=[]
                m_peak_amplitude[peak]=[]

            if len(peak) <= 5:
                ZrH_peak_centres[peak]=[]
                ZrH_peak_fwhm[peak]=[]
                ZrH_peak_amplitude[peak]=[]


            image_number=[]

    #### Modify list length so data has required shape ###

            for i in range(len(data)-1):

                if numpeaks > 1:
                    prefix = 'pv_' + str(n+1)
                if numpeaks == 1:
                    prefix = ''

                if len(peak) <= 5:

                    ttheta = data[i][reflections]['params_values'][prefix + 'center']
                    ZrH_peak_centres[peak].append(ttheta)

                    ttheta = data[i][reflections]['params_values'][prefix + 'amplitude']
                    ZrH_peak_amplitude[peak].append(ttheta)

                    ttheta = data[i][reflections]['params_values'][prefix + 'fwhm']
                    ZrH_peak_fwhm[peak].append(ttheta)



                if len(peak) > 5:
                    ttheta = data[i][reflections]['params_values'][prefix + 'center']
                    m_peak_centres[peak].append(ttheta)

                    ttheta = data[i][reflections]['params_values'][prefix + 'amplitude']
                    m_peak_amplitude[peak].append(ttheta)

                    ttheta = data[i][reflections]['params_values'][prefix + 'fwhm']
                    m_peak_fwhm[peak].append(ttheta)

            #print(peak_centres[reflection])

            if len(peak) <= 5:

                list_ZrH_peaks.append(cake + '_' + peak)

                x = list(ZrH_peak_centres.values())[0]
                fwhm = list(ZrH_peak_fwhm.values())[0]
                amp = list(ZrH_peak_amplitude.values())[0]
                if len(list(ZrH_peak_centres.values())) > 1:
                    raise Exception
                list_ZrH_peak_centres.append(x)
                list_ZrH_peak_fwhm.append(fwhm)
                list_ZrH_peak_amplitude.append(amp)


            else:

                list_matrix_peaks.append(cake + '_' + peak)

                y = list(m_peak_centres.values())[0]
                m_fwhm = list(m_peak_fwhm.values())[0]
                m_amp = list(m_peak_amplitude.values())[0]
                list_matrix_peak_centers.append(y)           
                list_matrix_peak_fwhm.append(m_fwhm)
                list_matrix_peak_amplitude.append(m_amp)
    
    return list_ZrH_peaks, list_ZrH_peak_centres, list_ZrH_peak_fwhm, list_ZrH_peak_amplitude, list_matrix_peaks, list_matrix_peak_centers, list_matrix_peak_fwhm, list_matrix_peak_amplitude


##data is now stored in two lists
##Peak name and cake are in 'list_blah_peaks'
##the location (in two theta) are in 'list_blah_peak_centres' which is a nested list containing all strain levels within top level list

##re-arrange data into format: [peak][step][cake-strain]

#Create list of cake numbers
#Create list of peak names

def order_peaks_n_cakes(list_matrix_peaks, list_ZrH_peaks):
    m_cakes_list = []
    m_peaks_list = []
    z_cakes_list = []
    z_peaks_list = []

    for i in list_matrix_peaks:
        peak = i.split('_')[1]
        m_peaks_list.append(peak)

    for i in list_matrix_peaks:
        cake = i.split('_')[0]
        m_cakes_list.append(cake)

    for i in list_ZrH_peaks:
        peak = i.split('_')[1]
        z_peaks_list.append(peak)

    for i in list_ZrH_peaks:
        peak = i.split('_')[0]
        z_cakes_list.append(peak)
        
    return m_cakes_list, m_peaks_list, z_cakes_list, z_peaks_list




#Seperate into the peaks

#list[first peak, second peak, third... ]

def all_mx_by_peak(m_peaks_list, m_cakes_list, m_strains_list, list_matrix_peak_fwhm, list_matrix_peak_amplitude, list_matrix_peak_centers):

    m_strains_by_peak = [[],[],[],[]]
    m_cakes_by_peak = [[],[],[],[]]
    m_peaks_by_peak = [[],[],[],[]]
    m_fwhm_by_peak = [[],[],[],[]]
    m_amp_by_peak = [[],[],[],[]]
    m_cntr_by_peak = [[],[],[],[]]

    for i in range(len(m_strains_by_peak)):
        for j in range(0,36):
            j = []
            m_strains_by_peak[i].append(j)
            m_cakes_by_peak[i].append(j)
            m_peaks_by_peak[i].append(j)
            m_fwhm_by_peak[i].append(j)
            m_amp_by_peak[i].append(j)
            m_cntr_by_peak[i].append(j)

    for i,j,k,l,m,n in zip(m_peaks_list, m_cakes_list, m_strains_list, list_matrix_peak_fwhm, list_matrix_peak_amplitude, list_matrix_peak_centers):
        index = int(j)    
        if i == '(01-10)':
            m_strains_by_peak[0][index] = k
            m_cakes_by_peak[0][index] = j
            m_peaks_by_peak[0][index] = i
            m_fwhm_by_peak[0][index] = l
            m_amp_by_peak[0][index] = m
            m_cntr_by_peak[0][index] = n

        if i == '(11-20)':
            m_strains_by_peak[1][index] = k
            m_cakes_by_peak[1][index] = j
            m_peaks_by_peak[1][index] = i
            m_fwhm_by_peak[1][index] = l
            m_amp_by_peak[1][index] = m
            m_cntr_by_peak[1][index] = n

        if i == '(02-20)':
            m_strains_by_peak[2][index] = k
            m_cakes_by_peak[2][index] = j
            m_peaks_by_peak[2][index] = i
            m_fwhm_by_peak[2][index] = l
            m_amp_by_peak[2][index] = m
            m_cntr_by_peak[2][index] = n

        if i == '(01-13)':
            m_strains_by_peak[3][index] = k
            m_cakes_by_peak[3][index] = j
            m_peaks_by_peak[3][index] = i
            m_fwhm_by_peak[3][index] = l
            m_amp_by_peak[3][index] = m
            m_cntr_by_peak[3][index] = n

    #print('cake order', len(m_cakes_by_peak), len(m_cakes_by_peak[0]), '\n', m_cakes_by_peak[0], '\n')

    ##Pad empty lists

    for i, j in enumerate(m_strains_by_peak):
        for k, l in enumerate(j):
            if l == []:
                m_strains_by_peak[i][k] = np.zeros(233)

    for i, j in enumerate(m_cakes_by_peak):
        for k, l in enumerate(j):
            if l == []:
                m_cakes_by_peak[i][k] = k

    for i, j in enumerate(m_fwhm_by_peak):
        for k, l in enumerate(j):
            if l == []:
                m_fwhm_by_peak[i][k] = np.zeros(233)

    for i, j in enumerate(m_amp_by_peak):
        for k, l in enumerate(j):
            if l == []:
                m_amp_by_peak[i][k] = np.zeros(233)

    for i, j in enumerate(m_cntr_by_peak):
        for k, l in enumerate(j):
            if l == []:
                m_cntr_by_peak[i][k] = np.zeros(233)

    ## Data now stored in 4 peak lists in order

    #print(len(m_strains_by_peak), len(m_strains_by_peak[0]), len(m_strains_by_peak[0][14]), m_strains_by_peak[0][14][56], '\n')
    #print(len(m_strains_by_peak), len(m_strains_by_peak[1]), len(m_strains_by_peak[1][0]), m_strains_by_peak[1][15][:], '\n')
    #print(len(m_strains_by_peak), len(m_strains_by_peak[2]), len(m_strains_by_peak[2][0]), m_strains_by_peak[2][15][:], '\n')
    #print(len(m_strains_by_peak), len(m_strains_by_peak[3]), len(m_strains_by_peak[3][0]), m_strains_by_peak[3][15][:], '\n')
    
    return m_strains_by_peak, m_cakes_by_peak, m_peaks_by_peak, m_fwhm_by_peak, m_amp_by_peak, m_cntr_by_peak




##re-arrange data into format: [peak][step][cake-strain]

def arrange_mx_pk_step(m_strains_by_peak, m_cakes_by_peak, m_peaks_by_peak, m_fwhm_by_peak, m_amp_by_peak, m_cntr_by_peak):
    
    m_strains2 = []
    m_fwhm2 = []
    m_amp2 = []
    m_cntr2 = []

    for peak in range(4): #generate 4 lists
        tmp1_strain = []
        tmp1_fwhm = []
        tmp1_amp = []
        tmp1_cntr = []

        for i in range(233): # generate 58 lists
            tmp2_strain = []
            tmp2_fwhm = []
            tmp2_amp = []
            tmp2_cntr = []
            tmp1_strain.append(tmp2_strain) #append 58 lists to tmp1
            tmp1_fwhm.append(tmp2_fwhm)
            tmp1_amp.append(tmp2_amp)
            tmp1_cntr.append(tmp2_cntr)

            for cake in range(36): # generate 36 lists
                #print(peak, cake, i, m_strains_by_peak[peak][cake][i])
                strain = m_strains_by_peak[peak][cake][i]
                tmp1_strain[i].append(strain) #append 36 lists to tmp1[i]

                fwhm = m_fwhm_by_peak[peak][cake][i]
                tmp1_fwhm[i].append(fwhm)

                amp = m_amp_by_peak[peak][cake][i]
                tmp1_amp[i].append(amp)

                cntr = m_cntr_by_peak[peak][cake][i]
                tmp1_cntr[i].append(cntr)

        m_strains2.append(tmp1_strain)
        m_fwhm2.append(tmp1_fwhm)
        m_amp2.append(tmp1_amp)
        m_cntr2.append(tmp1_cntr)
    return m_strains2, m_fwhm2, m_amp2, m_cntr2


##objective: re-arrange zrh data into format: [peak][cake][step]

def all_zrh_by_peak(zrh_pk_names, list_ZrH_peaks, list_ZrH_peak_centres, list_ZrH_peak_fwhm, list_ZrH_peak_amplitude):

    peak_cake_pos = []
    peak_centres = []
    z_fwhm = []
    z_amp = []

    for i in zrh_pk_names:
        new_list = []
        new_list2 = []
        tmp_fwhm = []
        tmp_amp = []
        for j,k,l,m in zip(list_ZrH_peaks, list_ZrH_peak_centres, list_ZrH_peak_fwhm, list_ZrH_peak_amplitude):
            if i in j:
                new_list.append(j.split('_')[0]) # append name 'cake_peak'
                new_list2.append(k) #append list of centres
                tmp_fwhm.append(l)
                tmp_amp.append(m)
        peak_cake_pos.append(new_list) #list of 'cake value' for peaks
        peak_centres.append(new_list2) #list of data matching shape of list: peak_cake_pos
        z_fwhm.append(tmp_fwhm)
        z_amp.append(z_amp)

    return peak_cake_pos, peak_centres, z_fwhm, z_amp




#Seperate zrh into the peaks

def arrange_zrh_by_peak(z_peaks_list, z_cakes_list, z_strains_list, list_ZrH_peak_fwhm, list_ZrH_peak_amplitude, list_ZrH_peak_centres):

    #Seperate zrh into the peaks

    #list[first peak, second peak, third... ]

    z_strains_by_peak = [[],[],[],[]]
    z_cakes_by_peak = [[],[],[],[]]
    z_peaks_by_peak = [[],[],[],[]]
    z_fwhm_by_peak = [[],[],[],[]]
    z_amp_by_peak = [[],[],[],[]]
    z_cntr_by_peak = [[],[],[],[]]

    for i in range(len(z_strains_by_peak)):
        for j in range(0,36):
            j = []
            z_strains_by_peak[i].append(j)
            z_cakes_by_peak[i].append(j)
            z_peaks_by_peak[i].append(j)
            z_fwhm_by_peak[i].append(j)
            z_amp_by_peak[i].append(j)
            z_cntr_by_peak[i].append(j)

    for i,j,k,l,m,n in zip(z_peaks_list, z_cakes_list, z_strains_list, list_ZrH_peak_fwhm, list_ZrH_peak_amplitude, list_ZrH_peak_centres):
        index = int(j)    
        if i == '(111)':
            z_strains_by_peak[0][index] = k
            z_cakes_by_peak[0][index] = j
            z_peaks_by_peak[0][index] = i
            z_fwhm_by_peak[0][index] = l
            z_amp_by_peak[0][index] = m
            z_cntr_by_peak[0][index] = n

        if i == '(220)':
            z_strains_by_peak[1][index] = k
            z_cakes_by_peak[1][index] = j
            z_peaks_by_peak[1][index] = i
            z_fwhm_by_peak[1][index] = l
            z_amp_by_peak[1][index] = m
            z_cntr_by_peak[1][index] = n

        if i == '(311)':
            z_strains_by_peak[2][index] = k
            z_cakes_by_peak[2][index] = j
            z_peaks_by_peak[2][index] = i
            z_fwhm_by_peak[2][index] = l
            z_amp_by_peak[2][index] = m
            z_cntr_by_peak[2][index] = n

    #print('cake order', len(m_cakes_by_peak), len(m_cakes_by_peak[0]), '\n', m_cakes_by_peak[0], '\n')

    ##Pad empty lists

    for i, j in enumerate(z_strains_by_peak):
        for k, l in enumerate(j):
            if l == []:
                z_strains_by_peak[i][k] = np.zeros(233)

    for i, j in enumerate(z_cakes_by_peak):
        for k, l in enumerate(j):
            if l == []:
                z_cakes_by_peak[i][k] = k

    for i, j in enumerate(z_fwhm_by_peak):
        for k, l in enumerate(j):
            if l == []:
                z_fwhm_by_peak[i][k] = np.zeros(233)

    for i, j in enumerate(z_amp_by_peak):
        for k, l in enumerate(j):
            if l == []:
                z_amp_by_peak[i][k] = np.zeros(233)

    for i, j in enumerate(z_cntr_by_peak):
        for k, l in enumerate(j):
            if l == []:
                z_cntr_by_peak[i][k] = np.zeros(233)

    ## Data now stored in 4 peak lists in order
    
    return z_strains_by_peak, z_cakes_by_peak, z_fwhm_by_peak, z_amp_by_peak, z_cntr_by_peak



##re-arrange zrh data into format: [peak][step][cake-strain]
def arrange_zrh_pk_step(z_strains_by_peak, z_fwhm_by_peak, z_amp_by_peak, z_cntr_by_peak):

    z_strains2 = []
    z_fwhm2 = []
    z_amp2 = []
    z_cntr2 = []

    for peak in range(len(z_strains_by_peak)): #generate 4 lists
        tmp1_strain = []
        tmp1_fwhm = []
        tmp1_amp = []
        tmp1_cntr = []

        for i in range(len(z_strains_by_peak[peak][0])): # generate 58 lists
            tmp2_strain = []
            tmp2_fwhm = []
            tmp2_amp = []
            tmp2_cntr = []
            tmp1_strain.append(tmp2_strain) #append 58 lists to tmp1
            tmp1_fwhm.append(tmp2_fwhm)
            tmp1_amp.append(tmp2_amp)
            tmp1_cntr.append(tmp2_cntr)

            for cake in range(len(z_strains_by_peak[peak])): # generate 36 lists
                strain = z_strains_by_peak[peak][cake][i]
                tmp1_strain[i].append(strain) #append 36 lists to tmp1[i]

                fwhm = z_fwhm_by_peak[peak][cake][i]
                tmp1_fwhm[i].append(fwhm)

                amp = z_amp_by_peak[peak][cake][i]
                tmp1_amp[i].append(amp)

                cntr = z_cntr_by_peak[peak][cake][i]
                tmp1_cntr[i].append(cntr)

        z_strains2.append(tmp1_strain)
        z_fwhm2.append(tmp1_fwhm)
        z_amp2.append(tmp1_amp)
        z_cntr2.append(tmp1_cntr)

    return z_strains2, z_fwhm2, z_amp2, z_cntr2




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

def steps_slide_bar():
    step_slider = widgets.IntSlider(value=0, min = 0, max = 234, step = 5)
    return step_slider

def cake_menu():
    
    cake_range = np.linspace(0,35, 36, dtype = int)
    cake_range_str = []

    for i in cake_range: 
        j = str(i)
        cake_range_str.append(j)

    Choose_cakes_menu = widgets.SelectMultiple(options=cake_range_str, value=['0'], description='Choose Cakes', disabled=False)
    
    return Choose_cakes_menu

def choose_cakes(cakes, z_dK, mx_dK):
    
    merge_cakes = []
    for cake in cakes:
        merge_cakes.append(int(cake))

    merge_dK_zrh = []
    merge_dK_mx = []

    #only use 220 & first matrix peak (for now)

    for i in merge_cakes:
        dK_zrh_cakes = z_dK[1][i]
        merge_dK_zrh.append(dK_zrh_cakes)
        #merge_dK_zrh[cake][step]

        dK_mx_cakes = mx_dK[1][i]
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
   
    
def plt_fwhm_vs_strain(cakes, z_dK, mx_dK, Eng_strain, z_strains_by_peak, m_strains_by_peak):

    merge_cakes, norm_dK_zrh1, norm_dK_mx1, tx_n_dK_zrh, tx_n_dK_mx, av_norm_dK_zrh, av_norm_dK_mx = choose_cakes(cakes, z_dK, mx_dK)
        
    
    fig, ax0 = plt.subplots(figsize=(18, 10))

    step_range = list(range(233))    
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
        ax0_y2.plot(Eng_strain, m_strains_by_peak[1][cake_val], color = matrix_colours[i], label = m_lab)

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
    ax0.set_ylim(-0.5, 2.5)
    #ax0.axvline(0.0273, ymin = 0, ymax = 0.5, ls = '--', color = 'dimgrey', lw = 1)
    #ax0.text(0.0273, 0.75, 'Onset of peak\nbroardening\nin Matrix', ha = 'center', wrap = True)
    #ax0.axvline(0.046, ymin = 0, ymax = 0.75, ls = '--', color = 'dimgrey', lw = 1)
    #ax0.text(0.046, 1.25, 'Onset of peak\nbroardening\nin Hydrides', ha = 'center', wrap = True)

    fig.tight_layout(rect=[0, 0, 1, 0.95])

    #fig.savefig(r'C:\Users\mbgnwob2\Dropbox (The University of Manchester)\2. Project\Python Script\Single Peak Fitting Script ORIGINAL - Copy\S2_Out_plots\W-H Plots\FWHM_&_L-Strain_vs_E-Strain_TDm3')
    
def plt_fwhm_vs_step(cakes, z_dK, mx_dK, Eng_strain, z_strains_by_peak, m_strains_by_peak):
    
    merge_cakes, norm_dK_zrh1, norm_dK_mx1, tx_n_dK_zrh, tx_n_dK_mx, av_norm_dK_zrh, av_norm_dK_mx = choose_cakes(cakes, z_dK, mx_dK)

    fig, ax0 = plt.subplots(figsize=(18, 10))

    step_range = list(range(len(Eng_strain)))
    matrix_colours = ['lightcyan', 'powderblue', 'lightblue','skyblue', 'lightskyblue', 'deepskyblue']
    matrix_dark = ['cornflowerblue','royalblue','dodgerblue','steelblue','lightslategrey','slateblue']
    hydride_colours = ['papayawhip', 'blanchedalmond', 'moccasin', 'navajowhite', 'wheat', 'burlywood']
    hydride_dark = ['peachpuff', 'coral', 'sandybrown', 'peru', 'chocolate', 'burlywood']

    fig.suptitle('Peak Broardening in Hydrided Zr4', fontsize=18, fontweight='bold');

    #matrix
    for i, (data, cake) in enumerate(zip(norm_dK_mx1, merge_cakes)):
        l_lab = 'FWHM (10-10) cake = ' + str(cake)
        ax0.plot(step_range, data, 'o', color = matrix_colours[i], label=l_lab)
    ax0.plot(step_range, medfilt(av_norm_dK_mx, 3), '--', color = 'midnightblue', linewidth=4, label='Av. $\Delta$FWHM Matrix (10-10)')
    ax0_y2 = ax0.twinx()
    for i, cake_val in enumerate(merge_cakes):
        m_lab = 'L_Strain Zr4 (10-10) Cake ' + str(cake_val)
        ax0_y2.plot(step_range, m_strains_by_peak[1][cake_val], color = matrix_dark[i], label = m_lab)

    #hydride
    for i, (data, cake) in enumerate(zip(norm_dK_zrh1, merge_cakes)):
        l_lab = 'FWHM (220) cake = ' + str(cake)
        ax0.plot(step_range, data, 'o', color = hydride_colours[i], label = l_lab)
    ax0.plot(step_range, medfilt(av_norm_dK_zrh, 3), '--', color = 'sienna', linewidth=4, label = 'Av. $\Delta$FWHM Zrh (220)')
    for i, cake_val in enumerate(merge_cakes):
        s_lab = 'L_Strain ZrH (220) Cake ' + str(cake_val)
        ax0_y2.plot(step_range, z_strains_by_peak[1][cake_val], color = hydride_dark[i], label = s_lab)

    ax0_y2.set_ylabel('microstrain (10$^{-6}$)', fontsize=15)
    ax0.legend(loc = 4,  bbox_to_anchor=(1.3, 0), fontsize=10)
    ax0_y2.legend(loc = 1, bbox_to_anchor=(1.3, 1),   fontsize=10)
    #title0 = ('(100)')
    #ax0.set_title(title0, fontsize=18)
    ax0.set_xlabel('Step', fontsize=15)
    ax0.set_ylabel(r'$\Delta$ FWHM (1/nm)', fontsize=15)
    ax0.set_xlim(0.00)
    ax0.set_ylim(-0.5, 2.5)
    ax0_y2.set_ylim(-5000, 25000)
    #ax0.axvline(23, ymin = 0, ymax = 0.5, ls = '--', color = 'dimgrey', lw = 1)
    #ax0.text(23, 0.8, 'Onset of peak\nbroardening\nin Matrix', ha = 'center', wrap = True)
    #ax0.axvline(30, ymin = 0, ymax = 0.75, ls = '--', color = 'dimgrey', lw = 1)
    #ax0.text(30, 1.26, 'Onset of peak\nbroardening\nin Hydrides', ha = 'center', wrap = True)

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    
    fig.savefig(r'C:\Users\mbgnwob2\Dropbox (The University of Manchester)\2. Project\Python Script\Single Peak Fitting Script ORIGINAL - Copy\S3_Out_plots\FWHM_&_L-Strain_vs_Step_LDm3')

def zrh_fwhm_plt_by_steps(step, pp_range_rad, z_fwhm2, zrh_pk_names, cake_nums):
    fig, ax = plt.subplots(figsize=(10, 10))
    ax = plt.subplot(111, projection='polar')
    ax.plot(pp_range_rad, z_fwhm2[0][step], label = zrh_pk_names[0])
    ax.plot(pp_range_rad, z_fwhm2[1][step], label = zrh_pk_names[1])
    ax.plot(pp_range_rad, z_fwhm2[2][step], label = zrh_pk_names[2])
    #ax.plot(pp_range_rad, z_fwhm2[3][0], label = zrh_pk_names[3])
    ax.legend(loc = 0, bbox_to_anchor = (1.2, 1.2))
    ax.set_rlabel_position(35)  # Move radial labels away from plotted line
    ax.set_theta_direction(-1)
    ax.set_theta_zero_location('N', offset=-0)
    ax.set_xticks(np.pi/180. * np.linspace(0, 360, 36, endpoint=False))
    ax.set_rmax(0.1)
    ax.set_rmin(0)
    ax.set_xticklabels(cake_nums)
    ax.set_rorigin(-0.02)
    plt.show()
    
def zrh_amp_plt_by_steps(step, pp_range_rad, z_amp2, zrh_pk_names, cake_nums):
    fig, ax = plt.subplots(figsize=(10, 10))
    ax = plt.subplot(111, projection='polar')
    ax.plot(pp_range_rad, z_amp2[0][step], label = zrh_pk_names[0])
    ax.plot(pp_range_rad, z_amp2[1][step], label = zrh_pk_names[1])
    ax.plot(pp_range_rad, z_amp2[2][step], '.', label = zrh_pk_names[2])
    #ax.plot(pp_range_rad, z_fwhm2[3][0], label = zrh_pk_names[3])
    ax.legend(loc = 0, bbox_to_anchor = (1.2, 1.2))
    ax.set_rlabel_position(35)  # Move radial labels away from plotted line
    ax.set_theta_direction(-1)
    ax.set_theta_zero_location('N', offset=-0)
    ax.set_xticks(np.pi/180. * np.linspace(0, 360, 36, endpoint=False))
    ax.set_rmax(700)
    ax.set_rmin(0)
    ax.set_xticklabels(cake_nums)
    #ax.set_rorigin(-0.02)
    plt.show()
    
def mx_amp_plt_by_steps(step, pp_range_rad, m_amp2, matrix_peak_names, cake_nums):
    fig, ax = plt.subplots(figsize=(10, 10))
    ax = plt.subplot(111, projection='polar')
    ax.plot(pp_range_rad, m_amp2[0][step], label = matrix_peak_names[0])
    ax.plot(pp_range_rad, m_amp2[1][step], label = matrix_peak_names[1])
    #ax.plot(pp_range_rad, m_fwhm2[2][step], label = matrix_peak_names[2])
    #ax.plot(pp_range_rad, m_strains2[3][step], label = matrix_peak_names[3])
    ax.legend(loc = 0, bbox_to_anchor = (1.2, 1.2))
    ax.set_rlabel_position(35)  # Move radial labels away from plotted line
    ax.set_theta_direction(-1)
    ax.set_theta_zero_location('N', offset=-0)
    ax.set_xticks(np.pi/180. * np.linspace(0, 360, 36, endpoint=False))
    ax.set_rmax(12200)
    ax.set_rmin(0)
    ax.set_xticklabels(cake_nums)
    #ax.set_rorigin(-2000)
    plt.show()
    
def mx_fwhm_plt_by_steps(step, pp_range_rad, m_fwhm2, matrix_peak_names, cake_nums):
    fig, ax = plt.subplots(figsize=(10, 10))
    ax = plt.subplot(111, projection='polar')
    ax.plot(pp_range_rad, m_fwhm2[0][step], label = matrix_peak_names[0])
    ax.plot(pp_range_rad, m_fwhm2[1][step], label = matrix_peak_names[1])
    #ax.plot(pp_range_rad, m_fwhm2[2][step], label = matrix_peak_names[2])
    #ax.plot(pp_range_rad, m_strains2[3][step], label = matrix_peak_names[3])
    ax.legend(loc = 0, bbox_to_anchor = (1.2, 1.2))
    ax.set_rlabel_position(35)  # Move radial labels away from plotted line
    ax.set_theta_direction(-1)
    ax.set_theta_zero_location('N', offset=-0)
    ax.set_xticks(np.pi/180. * np.linspace(0, 360, 36, endpoint=False))
    ax.set_rmax(0.055)
    ax.set_rmin(0)
    ax.set_xticklabels(cake_nums)
    #ax.set_rorigin(-2000)
    plt.show()

    
    
    
## Calc K and dK over the full strain range
# H_sqrd & C are constant BUT K and dK will vary

def WH_parameters_K_dK(z_cntr2, z_fwhm2, m_cntr2, m_fwhm2, wavelength):

    K_zrh = []
    dK_zrh = []

    K_mx = []
    dK_mx = []


    # generate 3 empty lists for "peaks"
    for i in range(3):
        peak = int(i)

        K_zrh_tmp = []
        dK_zrh_tmp = []

        K_zrh.append(K_zrh_tmp)
        dK_zrh.append(dK_zrh_tmp)

        # generate and append 36 cake level lists into each peak level list
        for k in range(36):
            cake = int(k)

            K_zrh_tmp2 = []
            dK_zrh_tmp2 = []
            K_zrh[peak].append(K_zrh_tmp2)
            dK_zrh[peak].append(dK_zrh_tmp2)    

            for j in range(len(z_cntr2[0])):
                step = int(j)
                pk_cntr_zrh = z_cntr2[peak][step][cake]
                pk_fwhm_zrh = z_fwhm2[peak][step][cake]

                k_zirk_h = K(pk_cntr_zrh, wavelength)
                dk_zirk_h = delta_K(pk_cntr_zrh, pk_fwhm_zrh, wavelength)

                K_zrh[peak][cake].append(k_zirk_h)       
                dK_zrh[peak][cake].append(dk_zirk_h)


    # generate 4 empty lists for "peaks"
    for i in range(4):
        peak = int(i)

        K_mx_tmp = []
        dK_mx_tmp = []

        K_mx.append(K_mx_tmp)
        dK_mx.append(dK_mx_tmp)

        # generate and append 36 cake level lists into each peak level list
        for k in range(36):
            cake = int(k)

            K_mx_tmp2 = []
            dK_mx_tmp2 = []
            K_mx[peak].append(K_mx_tmp2)
            dK_mx[peak].append(dK_mx_tmp2)    

            for j in range(233):
                step = int(j)
                pk_cntr_mx = m_cntr2[peak][step][cake]
                pk_fwhm_mx = m_fwhm2[peak][step][cake]
                k_mx = K(pk_cntr_mx, wavelength)
                dk_mx = delta_K(pk_cntr_mx, pk_fwhm_mx, wavelength)

                K_mx[peak][cake].append(k_mx)
                dK_mx[peak][cake].append(dk_mx)
                
    return K_zrh, dK_zrh, K_mx, dK_mx



##calc K*sqrt(C)

def K_root_C(K_zrh, C_zrh, K_mx, C_mx):
    K_rootC_zrh = []
    K_rootC_mx = []

    for h, i, j in zip(range(3), K_zrh, C_zrh):
        peak = int(h)
        tmp = []
        K_rootC_zrh.append(tmp)

        for k, l in enumerate(i):
            cake = int(k)

            tmp2 = []
            K_rootC_zrh[peak].append(tmp2)

            for m in l:
                rt_j_zrh = math.sqrt(j)
                K_sqrtC_zrh = m*rt_j_zrh
                K_rootC_zrh[peak][cake].append(K_sqrtC_zrh)

    for h, i, j in zip(range(4), K_mx, C_mx):
        peak = int(h)
        tmp = []
        K_rootC_mx.append(tmp)

        for k, l in enumerate(i):
            cake = int(k)

            tmp2 = []
            K_rootC_mx[peak].append(tmp2)

            for m in l:
                rt_j_mx = math.sqrt(j)
                K_sqrtC_mx = m*rt_j_mx
                K_rootC_mx[peak][cake].append(K_sqrtC_mx)
                
    return K_rootC_zrh, K_rootC_mx




    
def W_H_Plots(cake, step, zrh_pk_names, K_rootC_zrh, dK_zrh, matrix_peak_names, K_rootC_mx,dK_mx):

    fig, (ax0, ax1) = plt.subplots(nrows=1, ncols=2, sharey=False, figsize=(20, 8))

    for i, j in enumerate(zrh_pk_names):
        ax0.plot(K_rootC_zrh[i][cake][step], dK_zrh[i][cake][step], 'o', label=j)
        ax0.legend(loc = 2,  fontsize=10)
        title0 = ('Line Broardening in ZrH-Delta')
        ax0.set_title(title0, fontsize=18)
        ax0.set_xlabel(r'K$\sqrt{C}$ (1/nm)', fontsize=15)
        ax0.set_ylabel(r'FWHM (1/nm)', fontsize=15)
        ax0.set_ylim(0, 15)

    for i, j in enumerate(matrix_peak_names[:-1]):
        ax1.plot(K_rootC_mx[i][cake][step], dK_mx[i][cake][step], 'o', label=j)
        ax1.legend(loc = 2, fontsize=10)
        title1 = ('Line Broardening in Zr-Alpha')
        ax1.set_title(title1, fontsize=18)
        ax1.set_xlabel(r'K$\sqrt{C}$ (1/nm)', fontsize=15)
        ax1.set_ylabel(r'FWHM (1/nm)', fontsize=15)

    #fig.savefig(r'C:\Users\mbgnwob2\Dropbox (The University of Manchester)\2. Project\Python Script\Single Peak Fitting Script ORIGINAL - Copy\S2_Out_plots\W-H Plots\W-H Plots')
