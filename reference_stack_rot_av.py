import glob
from EMAN2 import *

m = EMData("relionpath/MaskCreate/job004/mask.mrc") #edit this - path to mask
stack_filename = 'MT_reference_stack_masked_rotational_average_center_normalized_2p21_box300_unbinned.hdf' #edit this - output name

for item in range(6):
    pf = item+11
    if item == 2 or item == 3:
        stack_loc1 = (2*item)+1
        stack_loc2 = 2*item
    else:
        stack_loc1 = 2*item
        stack_loc2 = (2*item)+1
    data=EMData(f'{pf}pf_syn_ref_tubulin_only_lpf20_2p21Apix_box_300_inverted_normalize.mrc') #edit this - output of invert but fstring for pf number
    masked = data*m
    z_proj = masked.process('misc.directional_sum',{'axis':'z'})
    rot_av = z_proj.process('xform.applysym',{'sym':f'c{pf}'})
    rot_av.write_image(stack_filename, stack_loc1)
    flip = Transform()
    flip.set_rotation({"type":"XYZ","ztilt":0,"ytilt":180,"xtilt":0})#y 180
    rot_av.transform(flip)
    rot_av.write_image(stack_filename, stack_loc2)
