import re
import copy
import unittest

import numpy as np
import matplotlib.pyplot as mpl

# https://pypi.python.org/pypi/gpxpy
# https://github.com/tkrajina/gpxpy
# import gpxpy

def xy_to_heading_dist(easting, northing):
  """
  Convert cartesian coordinates to polar.
  heading is in degrees

       North (y)
         0
         ^
         |
  270 <--+--> 90   East (x)
         |
         V
        180
  """
  dist = np.sqrt(easting**2 + northing**2)
  heading = np.arctan2(easting, northing) * 180 / np.pi
  if heading < 0:
    heading += 360
  return (heading, dist)


class TestPolar(unittest.TestCase):
  def test1(self):
    heading, dist = xy_to_heading_dist(0, 10)
    self.assertAlmostEqual(heading, 0, places=None, delta=1e-7)
    self.assertAlmostEqual(dist, 10, places=None, delta=1e-7)

  def test2(self):
    heading, dist = xy_to_heading_dist(10, 0)
    self.assertAlmostEqual(heading, 90, places=None, delta=1e-7)
    self.assertAlmostEqual(dist, 10, places=None, delta=1e-7)

  def test3(self):
    heading, dist = xy_to_heading_dist(0, -10)
    self.assertAlmostEqual(heading, 180, places=None, delta=1e-7)
    self.assertAlmostEqual(dist, 10, places=None, delta=1e-7)

  def test4(self):
    heading, dist = xy_to_heading_dist(-10, 0)
    self.assertAlmostEqual(heading, 270, places=None, delta=1e-7)
    self.assertAlmostEqual(dist, 10, places=None, delta=1e-7)

  def test5(self):
    heading, dist = xy_to_heading_dist(10, 10)
    self.assertAlmostEqual(heading, 45, places=None, delta=1e-7)
    self.assertAlmostEqual(dist, 10*np.sqrt(2), places=None, delta=1e-7)


class ShortPath:
  def __init__(self):
    self.start_label = None
    self.finish_label = None
    self.heading_dist = []
    self.fractional_error = 0.03

  def set_start_label(self, label):
    self.start_label = label

  def set_finish_label(self, label):
    self.finish_label = label

  def set_fractional_error(self, fractional_error):
    self.fractional_error = fractional_error

  def add_step(self, degrees, distance, scale=1):
    self.heading_dist.append((degrees, distance*scale))

  def get_num_steps(self):
    return len(self.heading_dist)

  def rotate_clockwise(self, degrees_cw):
    self.heading_dist = [(deg + degrees_cw, dist) for deg, dist in self.heading_dist]

  def get_xy_offsets(self):
    delta_east = [dist * np.sin(deg * np.pi/180) for deg, dist in self.heading_dist]
    delta_north = [dist * np.cos(deg * np.pi/180) for deg, dist in self.heading_dist]
    return delta_east, delta_north

  def get_xy_path(self, with_zero=True):
    delta_east, delta_north = self.get_xy_offsets()
    if with_zero:
      delta_east = [0] + delta_east
      delta_north = [0] + delta_north
    eastings = np.cumsum(delta_east)
    northings = np.cumsum(delta_north)
    return eastings, northings

  def get_total_offset(self):
    eastings, northings = self.get_xy_offsets()
    return np.sum(eastings), np.sum(northings)

  def get_absolute_error(self):
    """
    Return path length traveled times fractional error
    """
    return np.sum([dist for heading, dist in self.heading_dist]) * self.fractional_error

  def correct_the_path(self, true_easting, true_northing):
    """
    Starting point of path is always 0, 0
    New ending point is (true_easting, true_northing)

    Error is assumed to be proportional to length of segment.
    Error is distributed proportional to length of segment,
    to minimize chi^2 for the corrected path.
    """
    delta_east, delta_north = self.get_xy_offsets()
    errors = [dist for deg, dist in self.heading_dist]
    fraction = np.array(errors) / np.sum(errors)
    raw_east = np.sum(delta_east)
    raw_north = np.sum(delta_north)
    correction_east = true_easting - raw_east
    correction_north = true_northing - raw_north
    correct_delta_east = delta_east + correction_east * fraction
    correct_delta_north = delta_north + correction_north * fraction
    self.heading_dist = [xy_to_heading_dist(e, n) for e, n in zip(correct_delta_east, correct_delta_north)]


class TestShortPath(unittest.TestCase):
  def test1(self):
    s = ShortPath()
    self.assertEqual(s.get_num_steps(), 0)
    s.add_step(0, 10)
    s.add_step(90, 10)
    s.add_step(180, 10)
    s.add_step(270, 10)
    self.assertEqual(s.get_num_steps(), 4)

    de, dn = s.get_total_offset()
    delta = 1e-7
    self.assertAlmostEqual(de, 0, places=None, delta=delta)
    self.assertAlmostEqual(dn, 0, places=None, delta=delta)

  def test2(self):
    s = ShortPath()
    s.add_step(0, 10)
    s.add_step(90, 10)
    s.add_step(180, 10)
    de, dn = s.get_total_offset()
    delta = 1e-7
    self.assertAlmostEqual(de, 10, places=None, delta=delta)
    self.assertAlmostEqual(dn, 0, places=None, delta=delta)

    s.rotate_clockwise(90)
    de, dn = s.get_total_offset()
    self.assertAlmostEqual(de, 0, places=None, delta=delta)
    self.assertAlmostEqual(dn, -10, places=None, delta=delta)

  def test3(self):
    s = ShortPath()
    s.add_step(0, 10)
    s.add_step(90, 10)
    s.add_step(0, 10)
    s.add_step(90, 10)
    delta = 1e-7
    de, dn = s.get_total_offset()
    self.assertAlmostEqual(de, 20, places=None, delta=delta)
    self.assertAlmostEqual(dn, 20, places=None, delta=delta)
    s.correct_the_path(20, 20)
    de, dn = s.get_total_offset()
    self.assertAlmostEqual(de, 20, places=None, delta=delta)
    self.assertAlmostEqual(dn, 20, places=None, delta=delta)

  def test4(self):
    s = ShortPath()
    s.add_step(0, 10)
    s.add_step(90, 10)
    s.add_step(0, 10)
    s.add_step(90, 10)
    delta = 1e-7
    s.correct_the_path(21, 22)
    de, dn = s.get_total_offset()
    self.assertAlmostEqual(de, 21, places=None, delta=delta)
    self.assertAlmostEqual(dn, 22, places=None, delta=delta)

  def test5(self):
    s = ShortPath()
    delta = 1e-7
    s.set_fractional_error(0.01)
    self.assertAlmostEqual(s.get_absolute_error(), 0, places=None, delta=delta)
    s.add_step(0, 10)
    s.add_step(90, 10)
    self.assertAlmostEqual(s.get_absolute_error(), 0.2, places=None, delta=delta)
    s.add_step(0, 10)
    s.add_step(90, 10)
    self.assertAlmostEqual(s.get_absolute_error(), 0.4, places=None, delta=delta)
    s.set_fractional_error(0.1)
    self.assertAlmostEqual(s.get_absolute_error(), 4.0, places=None, delta=delta)


def solve_paths(path_list, known_locations={'start':(0,0)}):
  """
  Take a list of ShortPaths, some with coincident endpoints, and preferably
  overcontrained, and find the least squares solution for all of the paths.
  Do not solve for any labels in known_locations; consider those fixed.
  """
  labels = [p.start_label for p in path_list] + \
    [p.end_label for p in path_list]
  labels = set(labels)
  labels.difference(known_locations.keys())  # Don't include known locations
  labels = sorted(list(labels))
  label_to_index = {label: i for i, label in enumerate(labels)}

  num_labels = len(labels)
  num_constraints = len(path_list) + 1

  if num_labels > num_constraints:
    raise Exception("Too many labels for number of constraints")

  constraint_matrix = np.zeros((num_constraints, num_labels))
  constraint_vec_east = np.zeros(num_constraints)
  constraint_vec_north = np.zeros(num_constraints)

  constraint_matrix[-1,label_to_index['start']] = 1
  constraint_vec_east[-1] = 0   # it was already zero
  constraint_vec_north[-1] = 0  # it was already zero

  for row, path in enumerate(path_list):
    abs_error = path.get_absolute_error()
    if abs_error == 0:
      raise Exception("Error: path length not allowed to be zero")
    weight = 1 / abs_error

    easting, northing = path.get_total_offset()
    constraint_vec_east[row] = easting * weight
    constraint_vec_north[row] = northing * weight

    if path.start_label in label_to_index:
      ind_start = label_to_index[path.start_label]
      constraint_matrix[row, ind_start] = -weight
    else:
      known_east, known_north = known_locations[path.start_label]
      constraint_vec_east[row] += known_east * weight
      constraint_vec_north[row] += known_north * weight

    if path.finish_label in label_to_index:
      ind_finish = label_to_index[path.finish_label]
      constraint_matrix[row, ind_finish] = weight
    else:
      known_east, known_north = known_locations[path.finish_label]
      constraint_vec_east[row] -= known_east * weight
      constraint_vec_north[row] -= known_north * weight

  loc_east,  resid, rank, s = np.linalg.lstsq(constraint_matrix, constraint_vec_east)
  print ('easting errors = ', s)
  loc_north, resid, rank, s = np.linalg.lstsq(constraint_matrix, constraint_vec_north)
  print ('northing errors = ', s)

  locations = {label: (east, north) for label, east, north in zip(labels, loc_east, loc_north)}
  locations.update(known_locations)

  corrected_paths = []
  for path in path_list:
    e1, n1 = locations[path.start_label]
    e2, n2 = locations[path.finih_label]
    path.correct_the_path(e2 - e1, n2 - n1)
    corrected_paths.append(path)

  return corrected_paths, locations


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
    print ('errors = ', s)
    loc_north, resid, rank, s = np.linalg.lstsq(constraint_matrix, constraint_vec_north)
    print ('errors = ', s)

    locations = {}
    for i, label in enumerate(labels):
        locations[label] = (loc_east[i], loc_north[i])

    #print constraint_matrix
    #print constraint_vec_east
        
    return locations 


def foo():
    with open('data.txt') as f:
        lines = f.readlines()
        f.close()

    # Default to 0 degrees declination.
    #declination = -11
    declination = 0

    # Default to angles in degrees
    use_degrees = True

    # Default to a pace count of 60 paces per 100 meters
    paces_per_100meters = 60.0

    label = 'start'
    locations = [(0,0)]
    checkpoints = []
    offset_east = 0
    offset_north = 0

    offsets = {}

    for line in lines:
        if re.search('^PACES_PER_100METERS', line):
            paces_per_100meters = float(line.split()[1])
            continue

        if re.search('^MAGNETIC_DECLINATION', line):
            declination = float(line.split()[1])
            continue

        if re.search('^MILS', line):
            use_degrees = False
            continue

        if re.search('^DEGREES', line):
            use_degrees = True

        if re.search('^#', line):
            continue

        if re.search('^@', line):
            old_label = copy.copy(label)
            # label is everything after the @ sign and before the first
            # whitespace
            label = line[1:].split()[0].strip()
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

                print ('pair: ', label_pair)

            offset_east = 0
            offset_north = 0
            short_list = []

            continue

        try:
            pair = line.split()[0:2]
            if len(pair) == 2:
                #print pair
                heading = float(pair[0])
                if use_degrees:
                  heading_radians = heading * np.pi / 180
                else:
                  # mils: 6400 mils = full circle
                  heading_radians = heading * np.pi / 3200

                # Convert paces to meters by multiplying by scale.
                scale = 100.0 / paces_per_100meters
                distance = float(pair[1]) * scale
                delta_east = distance * np.sin(heading_radians)
                delta_north = distance * np.cos(heading_radians)
                #locations.append(
                #    (locations[-1][0]+deast, locations[-1][1]+dnorth))
                offset_east += delta_east
                offset_north += delta_north
                short_list.append((delta_east, delta_north))
        except:
            print('line failed: ', line)

    pairs = offsets.keys()

    #for off in offsets:
    #    print off, offsets[off]
    #    print ''

    new_locs = solve_offsets(offsets)

    #mpl.figure(1, figsize=(10,8))
    #mpl.figure(1, figsize=(5,6))
    mpl.figure(1, figsize=(8.5,11))

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

        east, north = rotate_declination(east, north, declination=declination, scale=1)

        off_east = offsets[pair][2]
        off_north = offsets[pair][3]
        off_east, off_north = rotate_declination(off_east, off_north, declination=declination, scale=1)

        eastings = correct_the_path(off_east, east[0], east[1])
        northings = correct_the_path(off_north, north[0], north[1])
        mpl.plot(eastings, northings, '-k')

    east = [new_locs[x][0] for x in new_locs]
    north = [new_locs[x][1] for x in new_locs]
    east, north = rotate_declination(east, north, declination=declination, scale=1)
    mpl.plot(east, north, 'ob')
    mpl.xlabel('paces east')
    mpl.ylabel('paces north')

    mpl.axis('equal')

    #mpl.ylim([-800, 250])  # autoscaling wasn't working right

    mpl.tight_layout()
    #mpl.savefig('test_no_declination.png', transparent=True)
    #mpl.savefig('test_declination.png', transparent=True)
    mpl.savefig('test_declination.png', transparent=False, dpi=100)
    mpl.show()


if __name__ == "__main__":
    #foo()
    unittest.main()
        
