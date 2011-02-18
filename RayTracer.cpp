#include "CImg/CImg.h"
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
using namespace cimg_library;
using namespace std;
#define WIDTH 640
#define HEIGHT 400
#define PIXSIZE .05
#define INF 999999.9
#define FLEN 10
#define ERR 0.0001
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

struct Color{
	unsigned char r,g,b;
	Color(unsigned char _r, unsigned char _g, unsigned char _b){
		r = _r;
		g = _g;
		b = _b;
	}
	Color(){
		r = 0, g = 0, b = 0;
	}	
};

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
	double norm(){
		return sqrt(x*x+y*y+z*z);	
	}
};

struct Line{
	Point p;
	Vector v;
	Line(){
		p = Point();
		v = Vector();	
	}
	Line(Point &_p, Point &_q){
		p = _p;
		v = Vector(_q.x-p.x, _q.y-p.y, _q.z-p.z);
	}
	Line(Point &_p, Vector &_v){
		p = _p;
		v = _v;
	}
};

struct Sphere{
	Point c;
	double r;
	Color color;
	Sphere(){
		c = Point();
		r = 1;
		color = Color();
	}
	Sphere(Point _c, double _r, Color _color){
		c = _c;
		r = _r;
		color = _color;
	}
	Vector normal(Point &p){
		return Vector(c, p);	
	}
	Point intersection(Line &l){
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

Point operator+(const Point &p, const Vector &v){
	return Point(p.x + v.x, p.y + v.y, p.z + v.z);
}

Point operator*(const double &d, const Point &p){
	return Point(d*p.x, d*p.y, d*p.z);
}

Vector operator*(const double &d, const Vector &v){
	return Vector(d*v.x, d*v.y, d*v.z);
}

vector<Sphere> read(const char *const name){
	ifstream file;
	file.open(name);
	double x,y,z,rad;
	int r,g,b;
	vector<Sphere> obj;
	while (file >> x>>y>>z>>rad>>r>>g>>b){
		obj.push_back(Sphere(Point(x,y,z),rad,Color((char)r,(char)g,(char)b)));
	}
	return obj;
}

struct Tracer{
	Point *light;
	Point *focal;
	Vector *dir;
	Vector *horz;
	Vector *vert;
	CImg<unsigned char> *img;
	vector<Sphere> objects;

	Tracer(){
		light = new Point(1000.0, 1000.0, 0);
	        focal = new Point(0.0,0.0,-10.0);
		dir = new Vector(0.0, 0.0, 10.0);
        	horz = new Vector(1.0, 0.0, 0.0);
        	vert = new Vector(0.0, 1.0, 0.0);
        	img = new CImg<unsigned char>(WIDTH,HEIGHT,1,3,0);
	}
	vector<Sphere> read(const char *const name){
		ifstream file;
		file.open(name);
		double x,y,z,rad;
		int r,g,b;
		while (file >> x>>y>>z>>rad>>r>>g>>b){
			objects.push_back(Sphere(Point(x,y,z),rad,Color((char)r,(char)g,(char)b)));
		}
	}
	Color traceRay(Line &l){
		Color c(0,0,0);
		int indSource = -1;
		double minDist = INF;
		Point minPoint(INF, INF, INF);
		for (int o = 0; o < objects.size(); ++o){
			Point intersect = objects[o].intersection(l);
			double dist = intersect.dist(*focal);
			if (dist < minDist){
				minDist = dist;
				indSource = o;
				minPoint = intersect;
			}
		}
		if (indSource != -1){
			Line shadow(minPoint, *light);
			int indShadow = -1;
			double minDist = INF;
			Point minPoint(INF, INF, INF);

			for (int o = 0; o < objects.size(); ++o){
				Point intersect = objects[o].intersection(shadow);
				double dist = intersect.dist(*focal);
				if (dist < minDist){
					minDist = dist;
					indShadow = o;
					minPoint = intersect;
					break;
				}
			}
			if(indShadow != -1) {
				cout << minDist << endl;	
			}else{
				cout << "hey";
				c.r = objects[indSource].color.r;
				c.g = objects[indSource].color.g;
				c.b = objects[indSource].color.b;
			}
		}
		return c;
	}


	void drawScene(){
		for (int j = 0; j < HEIGHT; ++j){
			for(int i = 0; i < WIDTH; ++i){
				double xoff = (i-WIDTH/2)*PIXSIZE;
				double yoff = (HEIGHT/2-j)*PIXSIZE;
				Point px = (*focal)+(*dir)+xoff*(*horz) + yoff*(*vert);
				Line l(*focal, px);

				Color c = traceRay(l);

				(*img)(i,j,0) = c.r;
				(*img)(i,j,1) = c.g;
				(*img)(i,j,2) = c.b;

			}
		}
	}
};

int main() {
	Tracer t;
	t.read("in.txt");
	t.drawScene();
	t.img->save("out.bmp");
	return 0;
}
