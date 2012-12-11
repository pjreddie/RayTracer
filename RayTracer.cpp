#include "CImg/CImg.h"
#include <assert.h>
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
using namespace cimg_library;
using namespace std;
#define WIDTH 800
#define HEIGHT 600
#define PIXSIZE .02
#define INF 999999.9
#define FLEN 10
#define ERR 0.0001
#define AMBIENT 1
#define NUMWAVES 20
#define PI 3.14159265
#define MAXDEPTH 4
struct Point{
	double x,y,z;
	void print(){
		printf("(%f, %f, %f)\n", x,y,z);
	}
	bool operator== (const Point &p){return x == p.x && y == p.y && z == p.z;}
	Point(){
		x = 0., y = 0., z = 0.;
	}
	Point(double _x, double _y, double _z){
		x = _x;
		y = _y;
		z = _z;
	}
	double dist(Point &p){
		return sqrt(pow(p.x-x, 2) + pow(p.y-y, 2) + pow(p.z-z, 2));
	}
};
ostream &operator<<(ostream &stream, Point &p){
	stream << "(" << p.x << ", " << p.y << ", " << p.z << ")" << endl;
	return stream;
}
istream &operator>>(istream &stream, Point &p){
	stream >> p.x >> p.y >> p.z;
	return stream;
}

struct RGB{
	unsigned char r,g,b;
	RGB(unsigned char _r, unsigned char _g, unsigned char _b){
		r = _r;
		g = _g;
		b = _b;
	}
	RGB(){
		r = 0, g = 0, b = 0;
	}	
	void print(){
		cout <<(int) r << ", " <<(int) g << ", " <<(int) b << endl;
	}
};

struct Spectrum{
	vector<double> spectrum;
	Spectrum(double low, double high){
		spectrum = vector<double>(NUMWAVES, 0.0);
		double size = 370/NUMWAVES;
		for(int i = 0; i < NUMWAVES; ++i){
			if (i * size + 380 >= low && i * size + 380 <= high){
				spectrum[i] = 1.;
			}
		}
		//normalize(360);
	}
	Spectrum(vector<double> &_spectrum){
		spectrum = _spectrum;
	}
	Spectrum(){
		spectrum = vector<double>(NUMWAVES, 0.0);

	}
	Spectrum operator*(const Spectrum &s){
		vector<double> newspec(NUMWAVES, 0.0);
		for(int i = 0; i < NUMWAVES; ++i){
			newspec[i] = spectrum[i]*s.spectrum[i];
		}
		return Spectrum(newspec);
	}
	Spectrum operator*(const double d){
		vector<double> newspec(NUMWAVES, 0.0);
		for(int i = 0; i < NUMWAVES; ++i){
			newspec[i] = spectrum[i]*d;	
		}
		return Spectrum(newspec);
	}

	Spectrum operator+(const Spectrum &s){
		vector<double> newspec(NUMWAVES, 0.0);
		for(int i = 0; i < NUMWAVES; ++i){
			newspec[i] = spectrum[i]+s.spectrum[i];	
		}
		return Spectrum(newspec);
	}
	double integral(){
		double sum = 0.0;
		double size = 370/NUMWAVES;
		for(int i= 0; i < NUMWAVES; ++i){
			sum += spectrum[i]*size;
		}
		return sum;
	}
	void normalize(double i){
		double curri = integral();
		for(int j = 0; j < NUMWAVES; ++j){
			spectrum[j] = spectrum[j]*i/curri;
		}
	}
	void scaleTo(double i){
		double maxv = 0;
		for(int j = 0; j < NUMWAVES; ++j){
			maxv = max(maxv, spectrum[j]);
		}
		if(maxv > 0.0){
			for(int j = 0; j < NUMWAVES; ++j){
				spectrum[j] = spectrum[j]*i/maxv;
			}
		}
	}
	void print(){
		double maxVal = 0;
		for(int i = 0; i < NUMWAVES; ++i){
			maxVal = max(maxVal, spectrum[i]);
		}
		cout << "|------------------|\n";
		for(int i = 0; i < NUMWAVES; ++i){
			for(int j = 0; j < 20*(spectrum[i]/maxVal); ++j){
				cout << 'X';	
			}
			cout << endl;
		}
		cout << "|------------------|\n";
	}
};

Spectrum red;
Spectrum green;
Spectrum blue;

Spectrum magentaSpec(){
     Spectrum r(600,750);
     Spectrum b(380,500);
     Spectrum g(500, 600);
     Spectrum m = r+b*.5 + g*.1;
     m.scaleTo(1.);
     return m;
}   
Spectrum redSpec(){
     Spectrum r(600,750);
     Spectrum b(380,500);
     Spectrum g(500, 600);
     Spectrum m = r+b*.1 + g*.1;
     m.scaleTo(1.);
     return m;
}   
Spectrum blueSpec(){
     Spectrum r(600,750);
     Spectrum b(380,500);
     Spectrum g(500, 600);
     Spectrum m = r*.2+b+ g*.2;
     m.scaleTo(1.);
     return m;
}   
Spectrum greenSpec(){
     Spectrum r(600,750);
     Spectrum b(380,500);
     Spectrum g(500, 600);
     Spectrum m = r*.2+b*.2 + g;
     m.scaleTo(1.);
     return m;
}   


ostream &operator<<(ostream &stream, Spectrum &s){
	double maxVal = 0;
	for(int i = 0; i < NUMWAVES; ++i){
		maxVal = max(maxVal, s.spectrum[i]);
	}
	stream << "|-----Spectrum-----|\n";
	for(int i = 0; i < NUMWAVES; ++i){
		for(int j = 0; j < 20*(s.spectrum[i]/maxVal); ++j){
			stream << 'X';	
		}
		stream << endl;
	}
	stream << "|------------------|"<<endl;
	return stream;
}
istream &operator>>(istream &stream, Spectrum &s){
	string st;
	stream >> st;
	if(st.compare("Magenta")== 0){
		s = magentaSpec();
		return stream;
	}else if(st.compare("Red") == 0){
		s = redSpec();
		return stream;
	}else{
		for(int i = st.size()-1; i >= 0; --i){
			stream.putback(st[i]);
		}
	}
	double low, high;
	stream >> low >> high;
	s = Spectrum(low, high);
	return stream;
}

double gaussian(double mu, double s, double x){
	return (1/sqrt(2*PI*s*s))* exp(-(pow(x-mu, 2)/(2*s*s)));
}

Spectrum colorSpec(double mu, double s){
	double size = 370/NUMWAVES;
	vector<double> spec(NUMWAVES, 0);
	for(int i = 0; i < NUMWAVES; ++i){
		spec[i] = gaussian(mu, s, i*size + size/2 + 380);
	}
	double sum = 0;
	for(int i = 0; i < NUMWAVES; ++i) sum += size * spec[i];
	for(int i = 0; i < NUMWAVES; ++i) spec[i] = spec[i]/sum;
	return Spectrum(spec);
}


struct XYZ{
	double x;
	double y;
	double z;
	double maxv;

	XYZ(Spectrum s){
		x = (red*s).integral();
		y = (green*s).integral();
		z = (blue*s).integral();
		maxv = max(x, max(y,z));
		double intensity = s.integral();
		if (maxv != 0.0){
			x /= maxv;
			y /= maxv;
			z /= maxv;
		}
		x *= intensity;
		y *= intensity;
		z *= intensity;
	}
	XYZ(){
		x = 0;
		y = 0;
		z = 0;
		maxv = 1;
	}
	RGB rgb(){
		double r = min(x*255, 255.);
		double g = min(y*255, 255.);
		double b = min(z*255, 255.);

		return RGB((char) r,(char) g,(char) b);
	}
};

RGB operator*(const RGB &c, const double &d){
	double scale = min((255./max(c.r, max(c.g, c.b))), d);
	return RGB(c.r*scale, c.g*scale, c.b*scale);
}

struct Surface{
	Spectrum spectrum;
	double diffuse, specular, transmissive, refraction, smoothness;
	Surface(Spectrum _spectrum, double _diffuse, double _specular, double _transmissive, double _refraction, double _smoothness){
		spectrum = _spectrum;
		diffuse = _diffuse;
		specular = _specular;
		transmissive = _transmissive;
		refraction = _refraction;
		smoothness = _smoothness;
	}
	Surface(){
	}
};
ostream &operator<<(ostream &stream, Surface &s){
	stream << s.spectrum;
	stream << s.diffuse << ", " << s.specular <<", " << s.transmissive << ", " << s.refraction << ", " << s.smoothness << endl;
	return stream;
}
istream &operator>>(istream &stream, Surface &s){
	string st;
	stream >> st;
	if(st.compare("None") == 0) return stream;
	assert(st.compare("Spectrum") == 0);
	stream >> s.spectrum;
	stream >> st;
	assert (st.compare("Constants") == 0);
	stream >>s.diffuse >> s.specular >> s.transmissive >> s.refraction >> s.smoothness;
	return stream;
}

struct Vector{
	double x, y, z;
	Vector(){
		x = 0., y = 0., z = 0.;
	}
	Vector(double _x, double _y, double _z){
		x = _x;
		y = _y;
		z = _z;
	}
	Vector(Point &from, Point &to){
		x = to.x-from.x;
		y = to.y-from.y;
		z = to.z-from.z;
	}
	Vector(Point p){
		x = p.x;
		y = p.y;
		z = p.z;
	}
	double magnitude(){
		return sqrt(x*x+y*y+z*z);	
	}
	void normalize(){
		double m = magnitude();
		if (m != 0){
			x = x/m;
			y = y/m;
			z = z/m;
		}	
	}
};

ostream &operator<<(ostream &stream, Vector &p){
	stream << "(" << p.x << ", " << p.y << ", " << p.z << ")" << endl;
	return stream;
}
istream &operator>>(istream &stream, Vector &v){
	stream >> v.x >> v.y >> v.z;
	return stream;
}

Point operator+(const Point &p, const Vector &v){
	return Point(p.x + v.x, p.y + v.y, p.z + v.z);
}
Point operator-(const Point &p, const Point &q){
	return Point(p.x - q.x, p.y - q.y, p.z - q.z);
}
Point operator+(const Point &p, const Point &q){
	return Point(p.x + q.x, p.y + q.y, p.z + q.z);
}
Point operator*(const double &d, const Point &p){
	return Point(d*p.x, d*p.y, d*p.z);
}
Vector operator*(const double &d, const Vector &v){
	return Vector(d*v.x, d*v.y, d*v.z);
}
Vector operator+(const Vector &v, const Vector &t){
	return Vector(v.x+t.x, v.y+t.y, v.z+t.z);
}
Vector operator-(const Vector &v, const Vector &t){
	return Vector(v.x-t.x, v.y-t.y, v.z-t.z);
}
Vector operator-(const Vector &v){
	return Vector(-v.x, -v.y, -v.z);
}
double dot(Vector &n, Vector &m){
	return n.x*m.x + n.y*m.y + n.z*m.z;
}
double dot(Point &n, Vector &m){
	return n.x*m.x + n.y*m.y + n.z*m.z;
}
Vector cross(Vector &v, Vector &t){
	return Vector(v.y*t.z-v.z*t.y, v.z*t.x - v.x*t.z, v.x*t.y - v.y*t.x);
}
struct Line{
	Point p, q;
	Line(Point _p, Point _q){
		p = _p;
		q = _q;
	}
};


struct Ray{
	Point p;
	Vector v;
	Ray(){
		p = Point();
		v = Vector();	
	}
	Ray(Point &_p, Point &_q){
		p = _p;
		v = Vector(_q.x-p.x, _q.y-p.y, _q.z-p.z);
	}
	Ray(Point &_p, Vector &_v){
		p = _p;
		v = _v;
	}
};

struct Object{
	Surface surface;
	virtual Vector normal(Ray &r) = 0;
	virtual Point intersection(Ray &l) = 0;
};

struct Quadric : public Object{ 
	double A,B,C,D,E,F,G,H,I,J;
	Quadric(){
	A = B = C = D = E = F = G = H = I = J = 0;
	}
	Quadric(double _A, double _B, double _C, double _D, double _E, double _F, double _G, double _H, double _I, double _J){
	A = B = C = D = E = F = G = H = I = J = 0;
	}
	virtual Point intersection(Ray &r){
		r.v.normalize();
		Vector V1(r.v.x*r.v.x, r.v.y*r.v.y, r.v.z*r.v.z);
		Vector V2 = 2.*Vector(r.v.x*r.v.y, r.v.y*r.v.z, r.v.x*r.v.z);
		Vector V3 = 2.*Vector(r.v.x*r.p.x, r.v.y*r.p.y, r.v.z*r.p.z);
		Vector V4 = 2.*Vector(r.v.x*r.p.y + r.p.x*r.v.y, r.v.y*r.p.z + r.p.y*r.v.z, r.v.x*r.p.z + r.p.x*r.v.z);
		Vector V5 = 2.*Vector(r.v.x, r.v.y, r.v.z);
		Vector V6(r.p.x*r.p.x, r.p.y*r.p.y, r.p.z*r.p.z);
		Vector V7 = 2.*Vector(r.p.x*r.p.y, r.p.y*r.p.z, r.p.x*r.p.z);
		Vector V8 = 2.*Vector(r.p.x, r.p.y, r.p.z);

		Vector Q1(A,B,C);
		Vector Q2(D,E,F);
		Vector Q3(G,H,I);

		double a = dot(Q1, V1) + dot(Q2, V2);
		double b = dot(Q1, V3) + dot(Q2, V4) + dot(Q3, V5);
		double c = dot(Q1, V6) + dot(Q2, V7) + dot(Q3, V8) + J;

		double det = b*b-4*a*c;
		if(det < 0 || a == 0.) return Point(INF, INF, INF);
		double t0 = (-b + sqrt(det))/(2*a);
		double t1 = (-b -sqrt(det))/(2*a);

		if (t0 < ERR && t1 < ERR){
			return Point(INF, INF , INF);
		}
		if (t0 < ERR) t0 = INF;
		if(t1 < ERR) t1 = INF;
		t0 = min(t1, t0);
		return r.p + t0*r.v;
	}
	virtual Vector normal(Ray &r){
		double xn = 2*(A*r.p.x + D*r.p.y + F*r.p.z + G);	
		double yn = 2*(D*r.p.x + B*r.p.y + E*r.p.z + H);	
		double zn = 2*(F*r.p.x + E*r.p.y + C*r.p.z + I);	
		Vector n(xn, yn, zn);
		if(dot(n, r.v) > -ERR) return n;
		return -n;
	}
};
istream &operator>>(istream &stream, Quadric &q){
	string s;
	stream >> s;
	assert(s.compare("Constants") == 0);
	stream >> q.A >> q.B >> q.C >> q.D >> q.E >> q.F >> q.G >> q.H >> q.I >> q.J;
	stream >> s;
	assert(s.compare("Surface") == 0);
	stream >> q.surface;
	return stream;
}

struct Ellipsoid : public Quadric{
	Ellipsoid(){
	}
	Ellipsoid(Point p, double xr, double yr, double zr, Surface s){
		A = 1/(xr*xr);
		B = 1/(yr*yr);
		C = 1/(zr*zr);
		D = 0;
		E = 0;
		F = 0;
		G = -p.x/(xr*xr);
		H = -p.y/(yr*yr);
		I = -p.z/(zr*zr);
		J = p.x*p.x/(xr*xr) + p.y*p.y/(yr*yr) + p.z*p.z/(zr*zr)-1;
		surface = s;
	}
};
istream &operator>>(istream &stream, Ellipsoid &e){
	string s;
	stream >> s;
	assert(s.compare("Center") == 0);
	Point c;
	stream >> c;
	stream >> s;
	assert(s.compare("Radii") == 0);
	double x,y,z;
	stream >> x >> y >> z;
	stream >> s;
	assert(s.compare("Surface") == 0);
	Surface f;
	stream >> f;
	e = Ellipsoid(c, x,y,z, f);

	return stream;
}

struct Plane : public Object{
	Point p;
	Vector n;
	Plane(){
		p = Point();
		n = Vector();
	}
	Plane(Point _p, Vector _n, Surface _surface){
		p = _p;
		n = _n;
		surface = _surface;
	}
	virtual Vector normal(Ray &r){
		if (dot(n, r.v) > -ERR) return n;
		return -n;
	}
	virtual Point intersection(Ray &l){
		Point t = p-l.p;
		double num = dot(t, n);
		double denom = dot(l.v, n);
		if(denom == 0.0){
			if(num != 0.0){
				return Point(INF, INF, INF);
			}else{
				return l.p;
			}
		}
		double d = num/denom;
		if(d < ERR) return Point(INF, INF, INF);
		return l.p + d*l.v;	
	}
};

istream &operator>>(istream &stream, Plane &p){
	string s;
	stream >> s;
	assert(s.compare("Point") == 0);
	stream >> p.p;
	stream >> s;
	assert(s.compare("Vector") == 0);
	stream >> p.n;
	stream >> s;
	assert(s.compare("Surface") == 0);
	stream >> p.surface;
	return stream;
}

struct Point2d{
	double x, y;
	Point2d(){
		x = 0.;
		y = 0.;
	}
	Point2d(double _x, double _y){
		x = _x;
		y = _y;
	}
};
struct Line2d{
	Point2d p, q;
	Line2d(Point2d _p, Point2d _q){
		p = _p;
		q = _q;
	}
	Point2d intersection(Line2d &l){
		return Point2d(((p.x*q.y - p.y*q.x)*(l.p.x - l.q.x) - (p.x-q.x)*(l.p.x*l.q.y - l.p.y*l.q.x))/((p.x-q.x)*(l.p.y-l.q.y) - (p.y - q.y)*(l.p.x - l.q.x)),((p.x*q.y - p.y*q.x)*(l.p.y - l.q.y) - (p.y-q.y)*(l.p.x*l.q.y - l.p.y*l.q.x))/((p.x-q.x)*(l.p.y-l.q.y) - (p.y - q.y)*(l.p.x - l.q.x)));	
	}
	bool intersectSegment(Line2d &l){
		Point2d test = intersection(l);
		return (test.x <= max(p.x, q.x) && test.x >= min(p.x, q.x));
	}
};
struct Vector2d{
	double x, y;
	Vector2d(double _x, double _y){
		x = _x;
		y = _y;
	}
	Vector2d(Point2d &from, Point2d &to){
		x = to.x-from.x;
		y = to.y-from.y;
	}
	Vector2d(){
		x = 0.;
		y = 0.;
	}
};
double dot(Vector2d &v, Vector2d &t){
	return v.x*t.x + v.y*t.y;
}
struct Ray2d{
	Point2d p, q;
	Vector2d v;
	Ray2d(Point2d _p, Vector2d _v){
		p = _p;
		v = _v;
		q = Point2d(p.x + v.x, p.y + q.y);
	}
	bool intersectSegment(Line2d &l){
		if(((p.x-q.x)*(l.p.y-l.q.y) - (p.y - q.y)*(l.p.x - l.q.x)) < ERR && ((p.x-q.x)*(l.p.y-l.q.y) - (p.y - q.y)*(l.p.x - l.q.x)) > -ERR) return false;
		Point2d test(((p.x*q.y - p.y*q.x)*(l.p.x - l.q.x) - (p.x-q.x)*(l.p.x*l.q.y - l.p.y*l.q.x))/((p.x-q.x)*(l.p.y-l.q.y) - (p.y - q.y)*(l.p.x - l.q.x)),((p.x*q.y - p.y*q.x)*(l.p.y - l.q.y) - (p.y-q.y)*(l.p.x*l.q.y - l.p.y*l.q.x))/((p.x-q.x)*(l.p.y-l.q.y) - (p.y - q.y)*(l.p.x - l.q.x)));
		Vector2d testv(p, test);
		return(test.x <= max(l.p.x, l.q.x) && test.x >= min(l.p.x, l.q.x) && test.y <= max(l.p.y, l.q.y) && test.y >= min(l.p.y, l.q.y) && (dot(testv, v) > 0));

	}
};
Point2d collapsePoint(Point &p, int ind){	
	switch (ind){
		case 0:
			return Point2d(p.y, p.z);
		case 1:
			return Point2d(p.x, p.z);
		default:
			return Point2d(p.x, p.y);
	}
}
vector<Point2d> collapsePoints(vector<Point> &p, int ind){
	vector<Point2d> c;
	for(int i = 0; i < p.size(); ++i){
		switch (ind){
			case 0:
				c.push_back(Point2d(p[i].y, p[i].z));
				break;
			case 1:
				c.push_back(Point2d(p[i].x, p[i].z));
				break;
			case 2:
				c.push_back(Point2d(p[i].x, p[i].y));
				break;
			default:
				break;
		}
	}
	return c;
}
struct Polygon : Plane{
	vector<Point> bounds;
	vector<Point2d> collapsedPoints;
	int collapse;
	Polygon(){
		collapse = 0;
		p = Point();
		n = Vector();
		surface = Surface();
	}
	Polygon (vector<Point> _bounds, Surface _surface){
		bounds = _bounds;

		Vector v1(bounds[0], bounds[1]);
		Vector v2(bounds[0], bounds[2]);

		n = cross(v1, v2);
		p = bounds[0];
		surface = _surface;

		if (abs(n.x) >= max(abs(n.y), abs(n.z))){
			collapse = 0;
		}else if(abs(n.y) >= max(abs(n.x), abs(n.z))){
			collapse = 1;
		}else{
			collapse = 2;
		}
		collapsedPoints = collapsePoints(bounds, collapse);
		collapsedPoints.push_back(collapsedPoints[0]);
	}

	bool isIn(Point &p){	
		if(p.x == INF) return false;
		Point2d testC = collapsePoint(p, collapse);
		bool in = false;

		for(int i = 0; i < collapsedPoints.size()-1; ++i){
			if((testC.y <= collapsedPoints[i].y) == (testC.y > collapsedPoints[i+1].y) && testC.x - collapsedPoints[i].x - (testC.y-collapsedPoints[i].y)*(collapsedPoints[i+1].x - collapsedPoints[i].x)/(collapsedPoints[i+1].y - collapsedPoints[i].y)<0.0){
				in = !in;
			}
		}
		return in;
	}

	virtual Point intersection(Ray &r){
		Point test = Plane::intersection(r);
		if (isIn(test)) return test;
		return Point(INF, INF, INF);
	}
};
ostream &operator<<(ostream &stream, Polygon &p){
	stream << p.p;
	for(int i = 0; i <p.bounds.size(); ++i){
		stream << p.bounds[i];
	}
	stream << p.n;
	stream << p.surface;
}
istream &operator>>(istream &stream, Polygon &p){
	vector<Point> b;
	string s;
	while(stream >> s){
		if(s.compare("Point") != 0) break;
		Point next;
		stream >> next;
		b.push_back(next);
	}
	assert(s.compare("Surface") == 0);
	Surface surface;
	stream >> surface;
	p = Polygon(b, surface);
	return stream;
}

struct Mesh: public Object{
	vector<Polygon> polys;
	int lastHit;
	Mesh(){
		surface = Surface();
		polys.clear();
		lastHit = 0;
	}
	Mesh(vector<Polygon> _polys, Surface _surface){
		polys = _polys;
		surface = _surface;
		lastHit = 0;
	}
	virtual Point intersection(Ray &r){
		Point minPoint(INF, INF, INF);
		double minDist = INF;
		for(int i = 0; i < polys.size(); ++i){
			Point test = polys[i].intersection(r);
			double testd = r.p.dist(test);
			if (testd < minDist && testd > ERR){
				minDist = testd;
				minPoint = test;
				lastHit = i;
			}
		}
		return minPoint;
	}
	virtual Vector normal(Ray &r){
		return polys[lastHit].normal(r);
		for(int i = 0; i < polys.size(); ++i){
			if(polys[i].isIn(r.p)){
				return polys[i].normal(r);
			}
		}
		assert(false);
		return Vector(1,0,0);
	}
};
istream &operator>>(istream &stream, Mesh &m){
	vector<Polygon> polys;
	string s;
	while(stream >> s){
		if(s.compare("Polygon") != 0) break;
		Polygon next;
		stream >> next;
		polys.push_back(next);
	}
	assert(s.compare("Surface") == 0);
	Surface surface;
	stream >> surface;
	m = Mesh(polys, surface);
	return stream;
}
struct Sphere : public Object{
	Point c;
	double r;
	Sphere(){
		c = Point();
		r = 1;
		surface = Surface();
	}
	Sphere(Point _c, double _r, Surface _surface){
		c = _c;
		r = _r;
		surface = _surface;
	}
	virtual Vector normal(Ray &r){
		Vector n(c, r.p);
		if (dot(n, r.v) > 0) return n;
		return -n;
	}
	virtual Point intersection(Ray &l){
		double t = (-.5 * sqrt(pow(-2.*l.v.x*c.x + 2.*l.v.x*l.p.x - 2.*l.v.y*c.y + 2.*l.v.y*l.p.y - 2.*l.v.z*c.z + 2.*l.v.z*l.p.z, 2.) - 4.*(pow(l.v.x, 2.) + pow(l.v.y, 2.) + pow(l.v.z, 2.))*(pow(c.x, 2.) - 2.*(c.x*l.p.x) + pow(c.y, 2.) - 2.*c.y*l.p.y + pow(c.z, 2.) - 2.*c.z*l.p.z - pow(r, 2.) + pow(l.p.x, 2.) + pow(l.p.y, 2.) + pow(l.p.z, 2.))) + l.v.x*c.x - l.v.x*l.p.x + l.v.y*c.y - l.v.y*l.p.y + l.v.z*c.z - l.v.z*l.p.z)/(pow(l.v.x, 2.) + pow(l.v.y, 2.) + pow(l.v.z, 2.));
		double t2 = (.5 * sqrt(pow(-2.*l.v.x*c.x + 2.*l.v.x*l.p.x - 2.*l.v.y*c.y + 2.*l.v.y*l.p.y - 2.*l.v.z*c.z + 2.*l.v.z*l.p.z, 2.) - 4.*(pow(l.v.x, 2.) + pow(l.v.y, 2.) + pow(l.v.z, 2.))*(pow(c.x, 2.) - 2.*(c.x*l.p.x) + pow(c.y, 2.) - 2.*c.y*l.p.y + pow(c.z, 2.) - 2.*c.z*l.p.z - pow(r, 2.) + pow(l.p.x, 2.) + pow(l.p.y, 2.) + pow(l.p.z, 2.))) + l.v.x*c.x - l.v.x*l.p.x + l.v.y*c.y - l.v.y*l.p.y + l.v.z*c.z - l.v.z*l.p.z)/(pow(l.v.x, 2.) + pow(l.v.y, 2.) + pow(l.v.z, 2.));
		if (t < t2){
			if (t > ERR){
				return Point(l.p.x + t*l.v.x, l.p.y + t*l.v.y, l.p.z + t*l.v.z);
			}else if(t2 > ERR){
				return Point(l.p.x + t2*l.v.x, l.p.y + t2*l.v.y, l.p.z + t2*l.v.z);	
			}else{
				return Point(INF, INF, INF);
			}
		}else if(t2 > ERR){
			return Point(l.p.x + t2*l.v.x, l.p.y + t2*l.v.y, l.p.z + t2*l.v.z);
		}else if(t > ERR){
			return Point(l.p.x + t*l.v.x, l.p.y + t*l.v.y, l.p.z + t*l.v.z);
		}else{
			return Point(INF, INF, INF);
		}
	}
};
istream &operator>>(istream &stream, Sphere &p){
	string s;
	stream >> s;
	assert(s.compare("Point") == 0);
	stream >> p.c;
	stream >> s;
	assert(s.compare("Radius") == 0);
	stream >> p.r;
	stream >> s;
	assert(s.compare("Surface") == 0);
	stream >> p.surface;
	return stream;
}
struct Light{
	Point p;
	Spectrum spectrum;
	double intensity;
	Light(Point _p, Spectrum _spectrum, double _intensity){
		p = _p;
		spectrum = _spectrum;
		intensity = _intensity;
	}
	Light(){

	}
};

istream &operator >>(istream &stream, Light &l){
	string s;
	stream >> s;
	assert(s.compare("Point") == 0);
	stream >> l.p;
	stream >> s;
	assert(s.compare("Spectrum") == 0);
	stream >> l.spectrum;
	stream >> s;
	assert(s.compare("Intensity")== 0);
	stream >> l.intensity;
	l.spectrum.normalize(l.intensity);
	return stream;
}

Vector randVector(double jiggle){
	double x = (double) rand()/(double) RAND_MAX;
	double y = (double) rand()/(double) RAND_MAX;
	double z = (double) rand()/(double) RAND_MAX;
	return Vector(x/jiggle, y/jiggle, z/jiggle);
}

struct Tracer{
	Point *focal;
	Vector *dir;
	Vector *horz;
	Vector *vert;
	CImg<char> *img;
	vector<Object *> objects;
	vector<Light> lights;
	Spectrum ambient;
	int width;
	int height;
	double pixsize;

	Tracer(){
		focal = new Point(0.0,20.0,-30.0);
		dir = new Vector(0.0, 0.0, 1.0);
		horz = new Vector(1.0, 0.0, 0.0);
		vert = new Vector(0.0, 1.0, 0.0);
		ambient = Spectrum(380, 750);
		ambient.normalize(AMBIENT);
		width = WIDTH;
		height = HEIGHT;
		pixsize = .0015*500/WIDTH;
	}
	void read(const char *const name){
		ifstream file;
		file.open(name);
		string type;
		while (file >> type){
			if(type.compare("Sphere") == 0){
				Sphere *s = new Sphere();
				file >> *s;
				objects.push_back(s);
			}else if(type.compare("Plane") == 0){
				Plane *p = new Plane();
				file >> *p;
				objects.push_back(p);
			}else if(type.compare("Polygon") == 0){
				Polygon *p = new Polygon();
				file >> *p;
				objects.push_back(p);
			}else if(type.compare("Mesh") == 0){
				Mesh *m = new Mesh();
				file >> *m;
				objects.push_back(m);
			}else if(type.compare("Quadric") == 0){
				Quadric *q = new Quadric();
				file >> *q;
				objects.push_back(q);
			}else if(type.compare("Ellipsoid") == 0){
				Ellipsoid *e = new Ellipsoid();
				file >> *e;
				objects.push_back(e);
			}else if(type.compare("Light") == 0){
				Light l;
				file >> l;
				lights.push_back(l);
			}else if(type.compare("Ambient") == 0){
				double a;
				file >> a;
				ambient.normalize(a);
			}else if(type.compare("Height") == 0){
				file >> height;
			}else if(type.compare("Width") == 0){
				file >> width;
				pixsize = .003*500/width;
			}else if(type.compare("Pixel") == 0){
				file >> pixsize;
			}else if(type.compare("Focal") == 0){
				file >> (*focal);
			}

		}
	}
	Spectrum traceRay(Ray &l, int depth, int currObject){
		Spectrum returnLight;

		int indSource = -1;
		double minDist = INF;
		Point minPoint(INF, INF, INF);
		for (int o = 0; o < objects.size(); ++o){
			Point intersect = objects[o]->intersection(l);
			double dist = intersect.dist(l.p);
			if (dist < minDist){
				minDist = dist;
				indSource = o;
				minPoint = intersect;
			}
		}
		if (indSource != -1){
			int nextObject = indSource;
			if(currObject == indSource){
				nextObject = -1;
			}
			Vector viewVector = l.v;
			viewVector.normalize();
			Vector backwardView = -viewVector;
			Ray backwardRay(minPoint, backwardView);
			Vector surfaceNormal = objects[indSource]->normal(backwardRay);
			surfaceNormal.normalize();


			if (depth < MAXDEPTH && objects[indSource]->surface.specular > 0.0){
				Vector backwardReflection = viewVector - 2.*dot(viewVector, surfaceNormal)*surfaceNormal;
				backwardReflection.normalize();
				double cosangle = dot(backwardReflection, surfaceNormal);
				if(false && cosangle >= 0.){
					for(int i = 0; i < 4; ++i){
						Vector jiggled = backwardReflection + randVector(objects[indSource]->surface.smoothness/cosangle);
						Ray reflectedRay(minPoint, jiggled);
						Spectrum reflectedSpectrum = traceRay(reflectedRay, depth+1, currObject);

						reflectedSpectrum = reflectedSpectrum * objects[indSource]->surface.specular;
						returnLight = returnLight + reflectedSpectrum*.25;
					}
				}else{
						Ray reflectedRay(minPoint, backwardReflection);
						Spectrum reflectedSpectrum = traceRay(reflectedRay, depth+1, currObject);

						reflectedSpectrum = reflectedSpectrum * objects[indSource]->surface.specular;
						returnLight = returnLight + reflectedSpectrum;
					
				}
			}

			if (depth < MAXDEPTH && objects[indSource]->surface.transmissive > 0.0){
				double currRefraction = 1.0003;
				if (currObject != -1){
					currRefraction = objects[currObject]->surface.refraction;
				}
				double nextRefraction = 1.0003;
				if(nextObject != -1){
					nextRefraction = objects[nextObject]->surface.refraction;
				}
				double refractionVal = currRefraction/nextRefraction;
				Vector backwardView = -1.*viewVector;
				double cosi = dot(surfaceNormal, backwardView);
				Vector refractedVector = refractionVal*viewVector + (refractionVal*cosi - sqrt(1 + pow(refractionVal, 2) * (pow(cosi, 2) - 1)))*surfaceNormal;
				Ray refractedRay(minPoint, refractedVector);
				Spectrum refractedSpectrum = traceRay(refractedRay, depth+1, nextObject);

				refractedSpectrum = refractedSpectrum*objects[indSource]->surface.transmissive;
				returnLight = returnLight + refractedSpectrum;
			}

			for(int i = 0; i < lights.size(); ++i){
				Light light = lights[i];
				Ray shadow(minPoint, light.p);
				int indShadow = -1;
				double minDistS = light.p.dist(minPoint);
				Point minPointS(INF, INF, INF);

				for (int o = 0; o < objects.size(); ++o){
					Point intersect = objects[o]->intersection(shadow);
					double dist = intersect.dist(minPoint);
					if (dist < minDistS){
						minDistS = dist;
						indShadow = o;
						minPointS = intersect;
						break;
					}
				}

				if(indShadow == -1){
					Vector lightVector = shadow.v;
					lightVector.normalize();
					double diffuseIntensity = max(0., dot(lightVector, surfaceNormal));
					Vector reflection = lightVector - 2.*dot(lightVector, surfaceNormal)*surfaceNormal;

					double specularIntensity = pow(max(0.,dot(reflection, viewVector)),objects[indSource]->surface.smoothness)*objects[indSource]->surface.specular; 
					returnLight = returnLight + light.spectrum*specularIntensity;

					diffuseIntensity *= objects[indSource]-> surface.diffuse;
					Spectrum diffuse = objects[indSource]->surface.spectrum*light.spectrum*diffuseIntensity;
					returnLight = returnLight + diffuse;

				}
			}
			Spectrum diffuse = objects[indSource]->surface.spectrum * ambient * objects[indSource]->surface.diffuse;
			returnLight = returnLight + diffuse;
		}
		return returnLight;
	}


	void drawScene(){
		img = new CImg<char>(width,height,1,3,0);
		for (int j = 0; j < height; ++j){
			for(int i = 0; i < width; ++i){
				Spectrum out;
				for(int x = 0; x < 1; x += 2){
					for(int y = 0; y < 1; y += 2){
						double xoff = (i-width/2 + x/2)*pixsize;
						double yoff = (height/2-j + y/2)*pixsize;
						Point px = (*focal)+(*dir)+xoff*(*horz) + yoff*(*vert);
						Vector d (*focal, px);
						Ray r(px, d);
						out = out + (traceRay(r, 1, -1));
					}
				}
				XYZ outColor = XYZ(out);
				RGB rgb = outColor.rgb();
				(*img)(i,j,0) = rgb.r;
				(*img)(i,j,1) = rgb.g;
				(*img)(i,j,2) = rgb.b;
			}
		}
	}
};

int main(int argc, char *argv[]) {
	red = colorSpec(610, 25);
	green = colorSpec(540, 30);
	blue = colorSpec(445, 25);
	Tracer t;
	t.read(argv[1]);
	t.drawScene();
	t.img->save(argv[2]);

	return 0;
}
