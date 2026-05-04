#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 19 20:05:13 2023

@author: mollygravett
"""
#to add: 
# - arg parser so can print off all FFTs and x-projs or not
# - layer line analysis on its own?
 

from EMAN2 import *
import pandas as pd
import glob
from scipy.signal import find_peaks, peak_prominences
import numpy as np
import matplotlib.pyplot as plt
import csv

def pad_with(vector, pad_width, iaxis, kwargs):
    pad_value = kwargs.get('padder', 0)
    vector[:pad_width[0]] = pad_value
    vector[-pad_width[1]:] = pad_value

def av_ffts(ffts):
    fft = ffts[0]
    for item in range(1,len(ffts)):
        fft_next = ffts[item]
        fft = fft+fft_next
    return(fft)#/len(ffts) #can change to average if want to average rather than sum, leanne sums ffts

def check_pf(stack, references, MT, mask):
    top_classes = []
    top_scores = []
    shape = np.shape(stack[0].numpy())
    pf_array = np.zeros((12,len(stack),shape[0], shape[1]))
    for j in range(len(stack)):
        img = stack[j]
        img.process('xform.center')
        img.process('filter.highpass.gauss', {'cutoff_abs':0.01})
        score = []
        best=(-2,-1)
        for item in range(12):
            pf = item+(11-((item+(item%2))/2))
            rotated = img.process('xform.applysym',{'sym':f'c{pf}'})
            #rotated.write_image(f'pf{pf}_class{item}_MT{MT}_ptcl{j}.mrc')
            rotated = rotated.process('normalize')
            reference = references[item].process('xform.center')
            reference = reference.process('normalize')
            #aligned1 = rotated.align("translational", reference, {'masked':True, 'useflcf':1, 'maxshift':2})#, {'useflcf':1})#.process('normalize'), {'maxshift':5}) #not sure normalisation necessary
            #check1 = aligned1
            #check1.write_image(f'translate_class{item}_MT{MT}_ptcl{j}.mrc')
            aligned = rotated.align("rotational", reference)# , {'useflcf':1} .process('normalize'), {'maxshift':5}) #not sure normalisation necessary
            pf_array[item, j] = aligned.numpy()
            #check = aligned
            #check.write_image(f'rotate_class{item}_MT{MT}_ptcl{j}.mrc')
            c = aligned.cmp("ccc",reference,{"negative":0,"mask":mask})
            best=max(best,(c,item))
            score.append(c)
        #print(score)
        top_score = best[0]
        top_score_index = best[1]+1
        top_classes.append(top_score_index)
        top_scores.append(top_score)
        to_write = [MT,j, top_score_index, top_score]
        writer_h.writerow(to_write)
        av_pf_array = np.sum(pf_array, axis = 1)
    return top_classes, top_scores, av_pf_array

def rot_av(stack, pf, mask):
    middle = int(len(stack)/2)
    ref = stack[middle].process('xform.applysym',{'sym':f'c{pf}'})
    new = ref
    for img in range(middle):
        av = stack[img].process('xform.applysym',{'sym':f'c{pf}'})
        img_nxt = av.align("rotate_translate", ref)#, {"mask":mask})
        new = new+img_nxt
    for img in range(middle, len(stack)):
        av = stack[img].process('xform.applysym',{'sym':f'c{pf}'})
        img_nxt = av.align("rotate_translate", ref)#, {"mask":mask})
        new = new+img_nxt
    return(new)

def get_mts(ptcl_list):
    mt_list = []
    for item in ptcl_list:
        item_list = item.split(' ')
        mt_number = item_list[7]
        mt_list.append(mt_number)
    return(mt_list)

def skippedimage(imagename):
    try:
        data=EMData(imagename)
        return False
    except RuntimeError:
        return True

def get_window(index, lw, hi, size):
    lwin, hwin = index - lw, index + hi
    if lwin < 0 :
        lwin = 0
    if hwin > size:
        hwin = size
    return lwin, hwin

def window_fft(ffts, win_len=[2,3]):
    mt_len = len(ffts)
    windowed_avgs = []
    for ptcl_index in range(mt_len):
        low, hi = get_window(ptcl_index, win_len[0], win_len[1], mt_len)
        avg = Averagers.get('mean')
        for i in range(low, hi):
            avg.add_image(ffts[i])
        avg_img = avg.finish()
        windowed_avgs.append(avg_img)
    return windowed_avgs
    

def peak_analysis(fft, fourier_size, pixel_size, hack=[0]):
    """Distance between layer line peaks"""
    binned = []
    for column in fft.T:
        counts, bin_edges = np.histogram(column, bins=256)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        mean = np.sum(bin_centers * counts) / np.sum(counts)
        binned.append(mean)
    reduced_data = np.array(binned)
    peaks, _ = find_peaks(reduced_data) #find peaks    
    middle = int(len(reduced_data)/2)
    start = int((fourier_size*pixel_size)/33)
    end = int((fourier_size*pixel_size)/50)
    window_1 = reduced_data[middle-start:(middle-end)+1]
    window_2 = reduced_data[middle+end:middle+start+1]
    peaks_w1, _ = find_peaks(window_1)
    peaks_w2, _ = find_peaks(window_2)
    heights_w1 = window_1[peaks_w1]
    heights_w2 = window_2[peaks_w2]
    prominences_w1 = heights_w1*peak_prominences(window_1, peaks_w1)[0]
    prominences_w2 = heights_w2*peak_prominences(window_2, peaks_w2)[0]
    top_prominence_w1 = np.max(prominences_w1)
    sorted_w1 = np.sort(prominences_w1)[:-1]
    upper_1 = np.mean(sorted_w1) + (2*np.std(sorted_w1))
    loc_1 = peaks_w1[np.where(prominences_w1==top_prominence_w1)]+(middle-start)
    top_prominence_w2 = np.max(prominences_w2)
    sorted_w2 = np.sort(prominences_w2)[:-1]
    upper_2 = np.mean(sorted_w2) + (2*np.std(sorted_w2))
    loc_2 = peaks_w2[np.where(prominences_w2==top_prominence_w2)]+(middle+end)
    list_loc2 = np.array([loc_1,loc_2])
    diff = list_loc2-middle
    prominences = peak_prominences(reduced_data, peaks)[0] 
    av_prom = np.mean(prominences)
    sd_prom = np.std(prominences)
    top_prominence = np.sort(prominences)[-1:]
    top_prominences = np.sort([top_prominence_w1,top_prominence_w2,top_prominence])
    signal_power=np.square(reduced_data).mean()
    noise_power=np.square(reduced_data).std()
    #snr = 10*np.log10(signal_power/noise_power)
    #print(signal_power, noise_power, snr)
    signal=np.mean(reduced_data)
    noise=np.std(reduced_data)
    snr_pre1=signal_power/noise_power
    snr_pre2 = np.log10(signal_power/noise_power)
    snr = 10*np.log10(signal_power/noise_power)
    snr2 = np.where(signal==0,0,signal/noise) #https://github.com/scipy/scipy/issues/9097
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 6))
    fft_em = EMNumPy.numpy2em(fft)
    fft_em.write_image(f"fft_MT{hack[0]}.mrc")
    hack[0]+= 1
    ax1.imshow(fft, cmap='gray')
    ax1.set_title('Fourier transform')
    ax2.plot(reduced_data) 
    ax2.set_title('Intensity plot')
    ax2.set_xlabel('Distance (pixels)')
    ax2.set_ylabel('Grey value')
    middle = peaks[np.where(prominences==top_prominence)]
    list_loc = [loc_1,loc_2,middle]
    for loc in list_loc:
        ax2.plot(loc, reduced_data[loc],'o')
    return(list_loc, fig, av_prom, av_prom+(2*sd_prom), top_prominences, snr_pre1, snr_pre2, snr, snr2)

def check_symmetry(list_loc):
    """Check if equal prominence"""
    if len(list_loc[0])==2:
        #print('sym', list_loc)
        return('symmetrical')
    if len(list_loc[0])==1:
        #print('asym', list_loc)
        return('asymmetrical')
    else:
        return('peak unclear')

def lattice_spacing(list_loc, symmetry, pixel_size, fourier_size):
    """Calculates lattice spacing from layer lines"""
    if symmetry == 'symmetrical':
        diff = list_loc[0]-list_loc[2]
    if symmetry == 'asymmetrical':
        diff = list_loc[:2]-list_loc[2]
    #print(diff)
    distances = (fourier_size*pixel_size)/abs(diff)
    #print(distances)
    counter = 0
    for item in distances:
        if 38 <= item <= 45:
            counter = counter+0
        else:
            counter = counter+1
    if counter == 0:
        return(distances)
    else:
        return('Peaks not in usual range')

def calculate(fft, img_name, pixel_size, save_correct_peaks = False):
    """Overall function calculating lattice spacing from FFT"""
    distance_array = []
    errors = []
    fourier_size = 2400 #edit this - padded box size can edit to size you like
    locations, fig2, av_prom, sd_prom,top_prom,SNR_pre1, SNR_pre2, SNR, SNR2 = peak_analysis(fft, fourier_size, pix_size)
    symmetry = check_symmetry(locations)
    lattice_size = lattice_spacing(locations, symmetry, pixel_size, fourier_size)
    if type(lattice_size) == str:
        error = [img_name,lattice_size]
        errors.append(error)
        fig2.suptitle(img_name)
        fig2.savefig(img_name.replace('.mrc','fail_test.png'))
        plt.cla()
        plt.close(fig2)
    else:
        distance_array.append([img_name, float(lattice_size[0]), float(lattice_size[1]), av_prom,sd_prom,top_prom, SNR_pre1, SNR_pre2, SNR, SNR2])
        lattice = float(lattice_size[0])
        if save_correct_peaks == True:
            fig2.suptitle(img_name)
            fig2.savefig(img_name.replace('.mrc','success_test.png'))
            plt.close(fig2)
            plt.cla()
    writer_g.writerows(distance_array)
    if len(errors)>0:
        writer_f.writerows(errors)
    else:
        return lattice

def FFT_per_MT(df, apix):
    #get mask
    m = EMData("/relionpath/MaskCreate/job004/mask.mrc") #edit this - path to unbinned mask
    zmask = m.process('misc.directional_sum',{'axis':'z'})
    references = []
    for idx in range(12):
        references.append(EMData("./MT_reference_stack_masked_rotational_average_center_normalized_5p56_unbinned.hdf", idx)) #edit this - path to references
    for MT, group in df.groupby('rlnHelicalTubeID'):
        ffts = []
        z_projs = []
        ptcl_list= []
        for r_idx, row in group.iterrows():
            ptcl = row['ptcl_number']
            #print(ptcl,MT)
            image = f'../WT_tomo01_subtomo.{ptcl}' #edit this  - tomoname
            if skippedimage(image) == True:
                continue
            else:
                data=EMData(image)
            try:
                rot = 0#row['rlnAngleRot']
                tilt = row['rlnAngleTilt']
                psi = row['rlnAnglePsi']
            except KeyError:
                rot = 0#row['rlnAngleRotPrior']
                tilt = row['rlnAngleTiltPrior']
                psi = row['rlnAnglePsiPrior']
            #get rotations
            t = Transform()
            t.set_rotation({"type":"spider","phi":rot,"theta":tilt,"psi":psi})
            i = t.inverse()
            data.transform(i)
            #mask subtomo
            masked = data*m     
            #project in x
            x_proj = masked.process('misc.directional_sum',{'axis':'x'})
            x_proj_arr = x_proj.numpy()
            padded_box = 2400 
            box_size = 300 #edit this - box_size
            x_proj_pad = np.pad(x_proj_arr, int((padded_box-box_size)/2), pad_with)
            #x_proj.write_image(f"x_proj_{MT}_{ptcl}.mrc")
            x_proj_pad_em = EMNumPy.numpy2em(x_proj_pad)
            x_proj_pad_em.write_image(f"x_proj_{padded_box}_{MT}_{ptcl}.mrc")
            #do fft
            fft = x_proj_pad_em.do_fft()
            #shift fft
            fft_shifted = fft.process("xform.fourierorigin.tocenter")
            #get amplitudes
            amps=fft_shifted.get_fft_amplitude()
            #get absolute values
            ab_values=amps.process("math.absvalue")
            #write fourier trasnform to file
            #ab_values.write_image(f"fft_{MT}_{ptcl}.mrc")
            ffts.append(ab_values)
            #hel_distances.append(float(row['rlnHelicalTrackLengthAngst']))
            z_proj = masked.process('misc.directional_sum',{'axis':'z'})
            z_projs.append(z_proj)
            ptcl_list.append(ptcl)
        average = av_ffts(ffts) #actually adding up
        class_nums,scores, av_pf_array = check_pf(z_projs, references, MT, zmask)
        cmap = plt.get_cmap("tab20")  # Get tab20 colormap
        unique_y_values = range(1,13)
        color_map = {val: cmap(i / 20) for i, val in enumerate(unique_y_values)}
        colors = [color_map[val] for val in class_nums]
        plt.figure(figsize=(12, 6)) 
        for i in range(len(ptcl_list)):
            plt.plot(ptcl_list[i], class_nums[i], marker="o", linestyle="None", color=colors[i])
        labels = ["",
                  "11 pf minus", "11 pf plus", 
                  "12 pf minus", "12 pf plus", 
                  "13 pf minus", "13 pf plus", 
                  "14 pf minus", "14 pf plus", 
                  "15 pf minus", "15 pf plus", 
                  "16 pf minus", "16 pf plus",
                  ""]
        plt.yticks(range(14), labels) 
        plt.ylabel('PF num and direction')
        plt.xlabel('Ptcl number')
        plt.savefig(f'pfnum_direction_plot_MT{MT}.png')
        plt.close()
        imgs = av_pf_array[::2] #[np.array(img.numpy()) for img in av_pf_array[::2]]
        titles = ['c11 average','c12 average','c13 average','c14 average','c15 average','c16 average']
        global_min = min(img.min() for img in imgs)
        global_max = max(img.max() for img in imgs)
        fig, axes = plt.subplots(2, 3, figsize=(12, 8))
        for ax, img, title in zip(axes.ravel(), imgs, titles):
            ax.imshow(img, cmap="gray", vmin=global_min, vmax=global_max)
            ax.set_title(title)
            ax.axis("off")
        fig.savefig(f'rotational_averages_MT{MT}.png')
        plt.close(fig)
        np_fft = average.numpy()
        print(np.shape(np_fft)[0])
        av_name = f"fft_overall_average_MT_{MT}.mrc"
        crop_start = int((len(np_fft[0])/2)-(box_size))
        crop_end = int((len(np_fft[0])/2)+(box_size))
        np_fft_cropped = np_fft[crop_start:crop_end,crop_start:crop_end]
        calculate(np_fft_cropped,av_name,apix, save_correct_peaks=True)


file = open('../../WT_tomo01_unbinned_resaved_PtsAdded_particles.star','r') #edit this path to particles.star
ptcls = file.readlines()

def isfloat(num):
    try:
        float(num)
        return True
    except ValueError:
        return False

data = []
header = []
counter = 0
for item in ptcls:
    if item.startswith('_'):
        header.append(item.split(' ')[0][1:])
    else:
        value = item.split('\t')[0]
        if isfloat(value) == True:
            item_list = item.split('\t')
            counter = counter+1
            for num in range(len(item_list)):
                point = item_list[num]
                if isfloat(point) == True:
                    item_list[num] = float(point)
                else:
                    pass
            item_list.append("%03d" % counter) #edit this  - change to %03d or %02d depending on digits of tomox_stack-ptcl.mrc
            data.append(item_list)
        else:
            continue
header.append('ptcl_number')
df = pd.DataFrame(data, columns=header)
pix_size = 2.21 #edit this - pixel size

f = open('Error_test.csv', 'w', encoding='UTF8', newline='')
g = open('Distances_test.csv', 'w', encoding='UTF8', newline='')
h = open('PF_direction_cross_correlation.csv', 'w', encoding='UTF8', newline='')

header_suc = ['Image', 'Distance 1', 'Distance 2', 'Av prominence', 'SD Prominences','Prominences', 'SNR_pre1', 'SNR_pre2', 'SNR', 'SNR2']
header_err = ['Image', 'Error']
header_pf = ['MT', 'Ptcl number','PF class', 'Correlation coefficient']
writer_g = csv.writer(g)
writer_g.writerow(header_suc)
writer_f = csv.writer(f)
writer_f.writerow(header_err)
writer_h = csv.writer(h)
writer_h.writerow(header_pf)
FFT_per_MT(df, pix_size)
g.close()
f.close()
h.close()