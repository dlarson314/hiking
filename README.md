# hiking
Create maps of hiking paths, using pace counts and compass headings

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

