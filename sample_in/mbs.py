from math import *
class Point:
	def __init__(self, x, y, z):
		self.x = x;
		self.y = y;
		self.z = z;	
	def norm(self):
		mag = sqrt(self.x*self.x + self.y*self.y + self.z*self.z)
		self.x /= mag
		self.y /= mag
		self.z /= mag
	def __add__(self, other):
		return Point(self.x+other.x, self.y+other.y, self.z+other.z)
	def mult(self,m):
		return Point(self.x*m, self.y*m, self.z*m)
	def cross(self, v):
		return Point(self.y*v.z - self.z*v.y, self.z*v.x - self.x*v.z, self.x*v.y - self.y*v.x)
	def __neg__(self):
		return Point(-self.x, -self.y, -self.z)
def sphere(p, rad, n):
		print "Sphere"
		print "Point", p.x, p.y, p.z
		print "Radius", rad
		if(n %2 == 0):
			print "Surface Spectrum 380 550 Constants .3 .2 .7 1.4 100"
		else:
			print "Surface Spectrum Magenta Constants .3 .2 .7 1.4 100"
		print

def mbs(len, n, start, vec, orth):
	if n > 6 or len < .2:
		return
	sphere(start, len, n)
	mbs(len/2, n+1, start+vec.mult(1.5*len), vec, orth)

	right = orth
	left = -orth
	forward = vec.cross(orth)
	backward = -(vec.cross(orth))

	mbs(len/4, n+1, start+right.mult(1.25*len), right, vec)
	mbs(len/4, n+1, start+left.mult(1.25*len), left, vec)
	mbs(len/4, n+1, start+forward.mult(1.25*len), forward, vec)
	mbs(len/4, n+1, start+backward.mult(1.25*len), backward, vec)

print """Plane
Point 0 0 0
Vector 0 1 0
Surface Spectrum 380 750
Constants .6 .3 0 0 100

Width 2000
Height 2500
Ambient .2

Focal -10 50 -100

Light
        Point -500      500    -1000
        Spectrum        380 750
        Intensity       1

Light
        Point 1000      1000    -1000
        Spectrum        380 750
        Intensity       .8
"""
mbs(32., 0, Point(0,32,0), Point(0,1,0), Point(1,0,0))
