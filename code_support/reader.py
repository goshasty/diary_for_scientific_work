import numpy as np


def reader(file_name, skip_lines=0):
    '''
    read only 3 coords
    return np.array

    '''

    list_obj = list()
    with open(file_name, 'r') as fd_file:
        for line in fd_file:
            if skip_lines:
                skip_lines -= 1
                continue

            line = line.split()
            coords = np.zeros(shape=(3,))
            coords[0] = float(line[0])
            coords[1] = float(line[1])
            coords[2] = float(line[2])
            list_obj.append(coords)
    return np.array(list_obj)