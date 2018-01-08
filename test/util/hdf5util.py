import h5py
from h5py import h5z
import numpy as np

DEFLATE_LEVEL = 3
MAX_CHUNK_SIZE_1D = 64 * 64 * 2
MAX_CHUNK_SIZE_2D = 64



def calc_shape_size(shape):
    size = 1
    for n in shape:
        size *= n
    return size

def calc_chunk_shape(data_shape):
    chunk_shape = list(data_shape)
    rank = len(data_shape)
    if rank == 1:
        chunk_shape[0] = min(data_shape[0], MAX_CHUNK_SIZE_1D)
        return tuple(chunk_shape)
    elif rank == 2:
        if data_shape[0] == 1:
            chunk_shape[1] = min(data_shape[1], MAX_CHUNK_SIZE_1D);
            return tuple(chunk_shape)
        elif data_shape[1] == 1:
            chunk_shape[0] = min(data_shape[0], MAX_CHUNK_SIZE_1D);
            return tuple(chunk_shape)
        elif data_shape[0] >= MAX_CHUNK_SIZE_2D and data_shape[1] >= MAX_CHUNK_SIZE_2D:
            chunk_shape[0] = min(data_shape[0], MAX_CHUNK_SIZE_2D);
            chunk_shape[1] = min(data_shape[1], MAX_CHUNK_SIZE_2D);
            return tuple(chunk_shape)
        else:
            data_size = calc_shape_size(chunk_shape)
            while data_size > MAX_CHUNK_SIZE_1D:
                i = np.argmax(chunk_shape)
                chunk_shape[i] /= 2
                print '2:', chunk_shape
                data_size = calc_shape_size(chunk_shape)
            return tuple(chunk_shape)
    else:
        raise ValueError('Invalid rank.', rank)

def read_to_ndarray(file_path, dataset_name, get_attributes=False):
    with h5py.File(file_path, 'r') as f:
        dataset = f[dataset_name]

        # Converts to a numpy array.
        data = dataset[...]

        if get_attributes:
            attributes = {}
            for key in dataset.attrs.iterkeys():
                attributes[key] = dataset.attrs[key]
            return data, attributes
        else:
            return data

def write_ndarray(data, file_path, dataset_name, **attributes):
    with h5py.File(file_path, 'w') as f:
        compress = False
        if h5z.filter_avail(h5z.FILTER_DEFLATE):
            filter_info = h5z.get_filter_info(h5z.FILTER_DEFLATE)

            if ((filter_info & h5z.FILTER_CONFIG_ENCODE_ENABLED) and
                (filter_info & h5z.FILTER_CONFIG_DECODE_ENABLED)):

                compress = True

        if compress:
            chunk_shape = calc_chunk_shape(data.shape)
            dataset = f.create_dataset(dataset_name,
                                       data=data,
                                       chunks=chunk_shape,
                                       compression='gzip',
                                       compression_opts=DEFLATE_LEVEL)
        else:
            dataset = f.create_dataset(dataset_name, data=data)

        if len(attributes) > 0:
            for key in attributes:
                dataset.attrs[key] = attributes[key]
