import argparse
import re
import copy
import unittest

import numpy as np
import matplotlib.pyplot as mpl


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
    self.origin = (0,0)
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

  def get_xy_path(self):
    delta_east, delta_north = self.get_xy_offsets()
    delta_east = [0] + delta_east
    delta_north = [0] + delta_north
    eastings = np.cumsum(delta_east) + self.origin[0]
    northings = np.cumsum(delta_north) + self.origin[1]
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
    [p.finish_label for p in path_list]
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

  loc_east,  resid, rank, s = np.linalg.lstsq(constraint_matrix, constraint_vec_east, rcond=None)
  print ('easting errors = ', s)
  loc_north, resid, rank, s = np.linalg.lstsq(constraint_matrix, constraint_vec_north, rcond=None)
  print ('northing errors = ', s)

  locations = {label: (east, north) for label, east, north in zip(labels, loc_east, loc_north)}
  locations.update(known_locations)

  corrected_paths = []
  for path in path_list:
    e1, n1 = locations[path.start_label]
    e2, n2 = locations[path.finish_label]
    path.correct_the_path(e2 - e1, n2 - n1)
    path.origin = locations[path.start_label]
    corrected_paths.append(path)

  return corrected_paths, locations


def main(args):
  with open(args.filename) as f:
    lines = f.readlines()
    f.close()

  # Default to 0 degrees declination.
  declination = 0

  # Default to angles in degrees
  use_degrees = True

  known_locations={'start':(0,0)}
  label = 'start'
  units = ''

  solve_level=0
  solved_paths = []
  path_list = []
  current_path = ShortPath()

  mpl.figure(1, figsize=args.figsize)

  for line in lines:
    if re.search('^UNITS', line):
      units = line.split()[1]
      continue

    if re.search('^PACES_PER_100METERS', line):
      scale = 100 / float(line.split()[1])
      continue

    if re.search('^SCALE', line):
      scale = float(line.split()[1])
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

    if re.search('^SOLVE', line):
      corrected_paths, new_locations = solve_paths(path_list, known_locations=known_locations)
      solved_paths += corrected_paths
      path_list = []
      known_locations = new_locations
      labels = set()
      for path in corrected_paths:
        labels.add(path.start_label)
        labels.add(path.finish_label)
        eastings, northings = path.get_xy_path()
        mpl.plot(eastings, northings, '-', color='C%d'%solve_level)
      labels = list(labels)
      x = [known_locations[label][0] for label in labels]
      y = [known_locations[label][1] for label in labels]
      mpl.plot(x, y, '+k') #, color='C%d'%solve_level)
      solve_level += 1
      for label in labels:
        easting, northing = known_locations[label]
        print('%s\t%s\t%s' % (easting, northing, label))
      continue

    if re.search('^@', line):
      # label is everything after the @ sign and before the first
      # whitespace
      label = line[1:].split()[0].strip()
      if current_path.get_num_steps() > 0:
        current_path.set_finish_label(label)
        path_list.append(current_path)
      current_path = ShortPath()
      current_path.set_start_label(label)
      continue

    try:
      pair = line.split()[0:2]
      if len(pair) == 2:
        heading = float(pair[0])
        if use_degrees:
          heading_degrees = heading + declination
        else:
          # mils: 6400 mils = full circle
          heading_degrees = heading * 360 / 6400.0 + declination

        # Convert paces to meters by multiplying by scale.
        distance = float(pair[1])
        current_path.add_step(heading_degrees, distance, scale=scale)
    except:
      print('line failed: ', line)


  mpl.xlabel(units + ' east')
  mpl.ylabel(units + ' north')
  mpl.axis('equal')
  mpl.tight_layout()
  mpl.savefig(args.output, transparent=False, dpi=args.dpi)
  mpl.show()


if __name__ == "__main__":
  #unittest.main()

  parser = argparse.ArgumentParser(description="Plot bearings and pace counts")
  parser.add_argument('filename', help='Input filename')
  parser.add_argument('-o', '--output', help='Output filaname for figure')
  parser.add_argument('-s', '--figsize', help='width,height of figure in inches', default="11,8.5")
  parser.add_argument('--dpi', type=int, help='Dots per inch', default=100)

  args = parser.parse_args()

  try:
    height, width = args.figsize.split(',')
    height = float(height)
    width = float(width)
    args.figsize = (height, width)
  except:
    print("Failed to parse figsize")
    args.figsize = (11, 8.5)

  if args.output == None:
    args.output = args.filename + '.png'

  print(args)

  main(args)



