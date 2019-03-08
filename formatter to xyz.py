import numpy as np

coef_scale = 300000.0 / 66.93


def convert_one_obj(data_str):
    r = float(data_str[2])
    dec = float(data_str[0]) * np.pi / 180.0
    ra = float(data_str[1]) * np.pi / 180.0
    obj_id = int(data_str[3])
    mass = float(data_str[4])

    x = r * coef_scale * np.cos(dec) * np.cos(ra)
    y = r * coef_scale * np.cos(dec) * np.sin(ra)
    z = r * coef_scale * np.sin(ra)

    return [x, y, z, obj_id, mass]


def write_info(coords, fd_output):
    new_line = str(coords[0]) + '\t' + str(coords[1]) + '\t'\
               + str(coords[2]) + '\t' + str(coords[3]) + '\t' + str(coords[4])
    fd_output.write(new_line + '\n')


def reader(file_name):
    list_objs = list()
    with open(file_name, 'r') as fd_file:
        i = 0
        j = 0
        for data_str in fd_file.readlines():
            data_str = data_str.split()
            list_objs.append(convert_one_obj(data_str))
            if list_objs[i][0] == 0.0:
                j += 1
            print(list_objs[i])
            i += 1

    file_output = 'r 0.05 0.075 xyz.txt'
    with open(file_output, 'w') as fd_file:
        for obj in list_objs:
            write_info(obj, fd_file)
    print(i, j)

reader('data/raw/0.05 0.075.txt')
