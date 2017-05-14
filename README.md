# hiking
Create maps of hiking paths, using pace counts and compass headings

## Usage

The notes taken while hiking should be of this form:
```text
# This is a comment
@start
000 30 
@waypoint1
120 30
240 30
@start
```
Lines starting with `@` give a name to the point where you are standing.  Lines
containing two numbers give the compass bearing (with respect to magnetic north)
to the point you are about to walk to, and then the distance to that point in
paces.

To allow the code to correct for errors, hiking paths should contain loops,
where you revisit the same point at least twice.  Therefore, it is preferred to
have the `@start` waypoint show up at least twice in the list of waypoints.

The resulting map can look something like this:
![Resulting map](test_declination.png)

## To Do
- [ ] Include better error bars (fractional error on pace counts)
- [x] Let magnetic declination be read in as a keyword
- [x] Let pace count (paces per 100 meters) be read in as a keyword
- [x] Let decimal latitude and longitude to be added to a waypoint, if known
- [ ] Use a spherical Earth model to give lat/lon values to all points
- [ ] If multiple waypoints have lat/lon values, use those in the fitting process. 
- [ ] Output result in gpx format, for input into www.openorienteering.org
  software, or other mapping software.

