import mrcfile
import numpy as np
import glob

list_files = glob.glob('*box_300.mrc') #edit this - wildcard reference names

for mrc in list_files:
    image = mrcfile.open(mrc, mode="r")
    data = image.data.copy()  # Copy original data
    data = -data  # Multiply all values by -1
    max = np.max(data)
    data = data+max+0.2
    mrc_out = mrcfile.new(f'{mrc[:-4]}_inverted_normalize.mrc', overwrite=True)
    mrc_out.set_data(data.astype(np.float32))
    mrc_out.voxel_size = 2.17
    mrc_out.close()
