import re
import copy

import numpy as np
import matplotlib.pyplot as mpl


def solve_offsets(offsets):
    pairs = offsets.keys()
    labels = set()
    for p in pairs:
        labels.add(p[0])
        labels.add(p[1])
    labels = list(labels)
    labels.sort()

    indices = {}
    for i, lab in enumerate(labels):
        indices[lab] = i

    nrows = len(pairs) + 1
    ncols = len(labels)
    constraint_matrix = np.zeros((nrows, ncols))
    constraint_vec_east = np.zeros(nrows)
    constraint_vec_north = np.zeros(nrows)

    constraint_matrix[-1,indices['start']] = 1
    constraint_vec_east[-1] = 0  # not actually necessary
    constraint_vec_north[-1] = 0  # not actually necessary

    for i, pair in enumerate(pairs): 
        ind0 = indices[pair[0]]
        ind1 = indices[pair[1]]
        constraint_matrix[i, ind0] = -1
        constraint_matrix[i, ind1] = +1
        constraint_vec_east[i] = offsets[pair][0]
        constraint_vec_north[i] = offsets[pair][1]

    loc_east,  resid, rank, s = np.linalg.lstsq(constraint_matrix, constraint_vec_east)
    print 'errors = ', s
    loc_north, resid, rank, s = np.linalg.lstsq(constraint_matrix, constraint_vec_north)
    print 'errors = ', s

    locations = {}
    for i, label in enumerate(labels):
        locations[label] = (loc_east[i], loc_north[i])

    #print constraint_matrix
    #print constraint_vec_east
        
    return locations 


def correct_the_path(numpy_east, east1, east2):
    if np.abs(numpy_east[-1]) < 1e-6:
        numpy_east += east1
        numpy_east[-1] = east2
    else:
        numpy_east = (numpy_east / numpy_east[-1]) * (east2 - east1) + east1
    #numpy_east += east1

    return numpy_east


def rotate_declination(east, north, declination=-11):
    c = np.cos(declination*np.pi/180) 
    s = np.sin(declination*np.pi/180) 
    mat = np.array([[c,  s], 
                    [-s, c]])
    both = np.dot(mat, np.vstack((np.array(east), np.array(north))))
    return both[0,:], both[1,:]


def foo():
    with open('data.txt') as f:
        lines = f.readlines()
        f.close()

    declination = -11
    #declination = 0
    label = 'start'
    locations = [(0,0)]
    checkpoints = []
    offset_east = 0
    offset_north = 0

    offsets = {}

    for line in lines:
        if re.search('^#', line):
            continue

        if re.search('^@', line):
            old_label = copy.copy(label)
            label = line[1:].strip()
            #print label
            checkpoints.append(locations[-1])

            if offset_east != 0 or offset_north != 0:
                label_pair = (old_label, label)  # Don't sort the order
                numpy_east = np.array([p[0] for p in short_list])
                numpy_east = np.hstack((np.array([0]), np.cumsum(numpy_east)))
                numpy_north = np.array([p[1] for p in short_list])
                numpy_north = np.hstack((np.array([0]), np.cumsum(numpy_north)))
                offsets[label_pair] = (offset_east, offset_north, 
                    numpy_east, numpy_north)

                print 'pair: ', label_pair

            offset_east = 0
            offset_north = 0
            short_list = []

            continue

        try:
            pair = line.split()[0:2]
            if len(pair) == 2:
                #print pair
                heading = float(pair[0])
                heading_radians = heading * np.pi / 180
                distance = float(pair[1])
                delta_east = distance * np.sin(heading_radians)
                delta_north = distance * np.cos(heading_radians)
                #locations.append(
                #    (locations[-1][0]+deast, locations[-1][1]+dnorth))
                offset_east += delta_east
                offset_north += delta_north
                short_list.append((delta_east, delta_north))
        except:
            pass

    pairs = offsets.keys()

    #for off in offsets:
    #    print off, offsets[off]
    #    print ''

    new_locs = solve_offsets(offsets)

    mpl.figure(1, figsize=(10,8))

    #east = [x[0] for x in locations]
    #north = [x[1] for x in locations]
    #mpl.plot(east, north)

    #east = [x[0] for x in checkpoints]
    #north = [x[1] for x in checkpoints]
    #mpl.plot(east, north, 'ob')

    for pair in pairs:
        east = (new_locs[pair[0]][0], new_locs[pair[1]][0])
        north = (new_locs[pair[0]][1], new_locs[pair[1]][1])
        #mpl.plot(east, north, '-y')

        east, north = rotate_declination(east, north, declination=declination)

        off_east = offsets[pair][2]
        off_north = offsets[pair][3]
        off_east, off_north = rotate_declination(off_east, off_north, declination=declination)

        eastings = correct_the_path(off_east, east[0], east[1])
        northings = correct_the_path(off_north, north[0], north[1])
        mpl.plot(eastings, northings, '-k')

    east = [new_locs[x][0] for x in new_locs]
    north = [new_locs[x][1] for x in new_locs]
    east, north = rotate_declination(east, north, declination=declination)
    mpl.plot(east, north, 'ob')
    mpl.xlabel('paces east')
    mpl.ylabel('paces north')

    mpl.axis('equal')
    mpl.tight_layout()
    #mpl.savefig('test_no_declination.png', transparent=True)
    mpl.savefig('test_declination.png', transparent=True)
    mpl.show()


if __name__ == "__main__":
    foo()
        
