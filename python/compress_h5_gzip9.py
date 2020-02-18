#!/usr/bin/env python3

import sys
import h5py
import numpy as np
import os

DEST_DIR_NAME = "_dest/"

if __name__ == '__main__':
    if len(sys.argv) == 3:
        h5_file_name = sys.argv[1]
        dataset_name = sys.argv[2]
        print(h5_file_name)
    else:
        print("Usage: {} h5_file dataset".format(sys.argv[0]))
        sys.exit(1)

    if not os.path.exists(DEST_DIR_NAME):
        os.mkdir(DEST_DIR_NAME)

    # Read.
    with h5py.File(h5_file_name, "r") as f:
        dset = f[dataset_name]
        data = np.zeros(dset.shape)
        dset.read_direct(data)
        print("  chunks: {}".format(dset.chunks))
        chunks = dset.chunks

    # Save.
    with h5py.File(DEST_DIR_NAME + h5_file_name, "w") as f:
        # No compression.
        #dset = f.create_dataset(dataset_name, data=data, chunks=chunks)
        # LZF.
        #dset = f.create_dataset(dataset_name, data=data, chunks=chunks, compression='lzf')
        # GZIP 9.
        #dset = f.create_dataset(dataset_name, data=data, chunks=chunks, compression='gzip', compression_opts=9)
        # GZIP 9, auto-chunking.
        dset = f.create_dataset(dataset_name, data=data, chunks=True, compression='gzip', compression_opts=9)

    # Check.
    with h5py.File(DEST_DIR_NAME + h5_file_name, "r") as f:
        dset2 = f[dataset_name]
        data2 = np.zeros(dset2.shape)
        dset2.read_direct(data2)
        print("  chunks in new file: {}".format(dset2.chunks))
        if not np.array_equal(data, data2):
            print("ERROR: Data mismatch in file {}".format(DEST_DIR_NAME + h5_file_name))
            sys.exit(1)
