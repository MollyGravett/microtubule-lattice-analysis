#based on alister burt's script https://gist.github.com/alisterburt/eec5ba8c0acbaf3f12eb1d21cf2d53b5
#this if for use with imod .mod models of filaments that have already had points added with PEET
#it takes the points model and generates angular priors based on the vectors between the neighbouring points
#the rot angles are randomised per microtubule and the contour id from imod is taken as the relion helical id
#the helical track length is calculated
from typing import Tuple
import eulerangles
import imodmodel
import numpy as np
import pandas as pd
import starfile as starfile
from pydantic import BaseModel
import glob
import random


class Dipole(BaseModel):
    center: Tuple[float, float, float]
    north: Tuple[float, float, float]
    arb: Tuple[float, float, float]


def relion_euler_angles_from_dipole(dipole: Dipole) -> Tuple[float, float, float]:
    arbitrary_vec = np.array(dipole.arb) #random vector
    z = np.array(dipole.north) - np.array(dipole.center) #vector between neighbouring particles
    z /= np.linalg.norm(z) #unit vector
    y = np.cross(z, arbitrary_vec) #vector perpendicular to arb vec and z
    y /= np.linalg.norm(y) #unit vector y
    x = np.cross(y, z) #vector perpendicular to z and y
    x /= np.linalg.norm(x) #unit vector x
    rotation_matrix = np.zeros((3, 3)) #empty array
    rotation_matrix[:, 0] = x #add unit vectors to rot matrix
    rotation_matrix[:, 1] = y
    rotation_matrix[:, 2] = z
    rotation_matrix = np.linalg.pinv(rotation_matrix)  # relion expects ref2particle
    euler_angles = eulerangles.matrix2euler( #convert rotation matrix to euler
        rotation_matrices=rotation_matrix,
        axes='zyz',
        intrinsic=True,
        right_handed_rotation=True,
    )
    return np.squeeze(euler_angles)


def mod_file_to_rln3_df(model_file) -> pd.DataFrame:
    df = imodmodel.read(model_file)
    ptcl_no = len(df)
    bin_size = 1 #edit this -  bin size of picking 
    original_pix_size = 2.21 #edit this -  pixel size collected at
    binned_pix_ang = original_pix_size * bin_size #calculates unbinned vox size
    dipoles = []
    distances = []
    half = int(ptcl_no/2)
    #print(half)
    num_mts = df.contour_idx.nunique()
    tube_ids = df["contour_idx"].unique()
    rand_range = random.sample(range(num_mts), num_mts)
    counter_half_1 = 0
    for rand in rand_range:
        tube_id = tube_ids[rand]
        mt_ptcl_no = df['contour_idx'].value_counts()[tube_id]
        added = counter_half_1 + mt_ptcl_no
        if added <= half:
            counter_half_1 = added
            df.loc[df['contour_idx'] == tube_id, 'subset'] = 1
        else:
            df.loc[df['contour_idx'] == tube_id, 'subset'] = 2
    #print(df['subset'].value_counts()[1])
    #print(df['subset'].value_counts()[2])
    for contour_idx, group in df.groupby('contour_idx'):
        arbitrary_vec = [] #creates random vector
        for num in range(3):
            rand = random.uniform(-1,1)
            arbitrary_vec.append(rand)
        dist = 0
        distances.append(dist) #first distance 0
        for item in range(len(group)-1):
            c = group.iloc[item, :][['x', 'y', 'z']] #point 1
            n = group.iloc[item+1, :][['x', 'y', 'z']] #point 2
            length = np.linalg.norm(np.array(n) - np.array(c)) #distance between neighbouring points
            dist= dist+(length*binned_pix_ang) #unbinned distance from first point
            distances.append(dist)
            dipoles.append(Dipole(center=tuple(c), north=tuple(n), arb=tuple(arbitrary_vec)))
            if item == len(group)-2:
                dipoles.append(Dipole(center=tuple(c), north=tuple(n), arb=tuple(arbitrary_vec))) #match final particle angle to previous
    #xyz = np.array([(d.center[0], d.center[1], d.center[2]) for d in dipoles])
    euler_angles = np.array([relion_euler_angles_from_dipole(d) for d in dipoles])
    #print('dipoles', len(dipoles))
    #print('euler_angles', len(euler_angles))
    #print('x', len(df.loc[:,"x"]))
    tomoname = model_file.split('_unbinned')[0] #edit this -  works if your naming convention tomoname_MT*.mrc
    df_new = pd.DataFrame(
        {
            'rlnCoordinateX': df.loc[:,"x"]*bin_size,
            'rlnCoordinateY': df.loc[:,"y"]*bin_size,
            'rlnCoordinateZ': df.loc[:,"z"]*bin_size,
            'rlnHelicalTubeID': df.loc[:,"contour_idx"].astype(int)+1,
            'rlnAngleRotPrior': euler_angles[:, 0],
            'rlnAngleTiltPrior': euler_angles[:, 1],
            'rlnAnglePsiPrior': euler_angles[:, 2],
            'rlnHelicalTrackLengthAngst': distances,
            'rlnAnglePsiFlipRatio': [0.500000]*ptcl_no,
            'rlnTomoName': tomoname,
            'rlnRandomSubset': df["subset"].astype(int)
        }
    )
    return df_new


if __name__ == '__main__':
    mod_files = list(glob.glob('*PtsAdded_resaved.mod')) #edit this -  needs to be correct search in file
    dfs = []
    for mod_file in mod_files:
        df = mod_file_to_rln3_df(mod_file)
        dfs.append(df)
    df = pd.concat(dfs)
    starfile.write(df, 'WT_tomo01_unbinned_PtsAdded_resaved_particles.star') #edit this -  rename star file
