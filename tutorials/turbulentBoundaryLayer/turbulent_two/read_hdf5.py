'''
Script for reading and inspecting the files activewall.hdf5 and passivewall.hdf5.
'''

import h5py

# --- SETTINGS ---
input_file_path = "./activewall.hdf5"
# or
# input_file_path = "./passivewall.hdf5"

# --- Load HDF5 file ---
with h5py.File(input_file_path, 'r') as f:
    
    # get global attributes (roots attributes)
    for attr in f.attrs:
        print(f"{attr}: {f.attrs[attr]}")

    print("\n")
    print("Keys in the HDF5 file:")
    print(f.keys())

    for i, key in enumerate(f.keys()):
        
        print(f"---{key}---")
        # print(f"-attributes-")
        for attr in f[key].attrs:
            print(f"{attr}: {f[key].attrs[attr]}")

        print(f"-datasets-")
        for keys2 in f[key].keys():
            print(f"{keys2}")
        print("-----------------")