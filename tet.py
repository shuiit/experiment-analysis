import h5py

filename = 'H:/My Drive/dark 2022/csv_dark/2023_08_09_60ms.hdf5'
with h5py.File(filename, "r") as f:
    print("Keys: %s" % f.keys())