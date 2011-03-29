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
#define AMBIENT .1
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
	/*
	RGB rgb(){
		RGB color(0,0,0);
		Spectrum r1 = Spectrum(380,475);
		Spectrum r2 = Spectrum(570, 750);
		Spectrum b = Spectrum(380, 495);
		Spectrum g = Spectrum(508, 631);

		double rmag = r1.magnitude() + r2.magnitude();
		double gmag = g.magnitude();
		double bmag = b.magnitude();

		Spectrum rInt1 = intersection(r1);
		Spectrum rInt2 = intersection(r2);
		Spectrum gInt = intersection(g);
		Spectrum bInt = intersection(b);

		color.r = (unsigned char)((rInt1.magnitude() + rInt2.magnitude()) * 255/rmag);
		color.g = (unsigned char)((gInt.magnitude()) * 255/gmag);
		color.b = (unsigned char)((bInt.magnitude()) * 255/bmag);
		return color;
	}
	*/
};

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

Spectrum red;
Spectrum green;
Spectrum blue;

RGB toRGB(Spectrum &s){
	double size = 370/NUMWAVES;
	double r = (red*s).integral();
	double g = (green*s).integral();
	double b = (blue*s).integral();
	double mag = max(r,max(g,b));
	r /= mag;
	g /= mag;
	b /= mag;
	return RGB((char) (r*255),(char) (g*255),(char) (b*255));
}

RGB operator*(const RGB &c, const double &d){
	double scale = min((255./max(c.r, max(c.g, c.b))), d);
	return RGB(c.r*scale, c.g*scale, c.b*scale);
}

struct Surface{
	Spectrum spectrum;
	double diffuse, specular, transmissive, refraction;
	Surface(Spectrum _spectrum, double _diffuse, double _specular, double _transmissive, double _refraction){
		spectrum = _spectrum;
		diffuse = _diffuse;
		specular = _specular;
		transmissive = _transmissive;
		refraction = _refraction;
	}
	Surface(){
	}
};
struct Photon{
	Spectrum spectrum;
	double intensity;
	Photon(double low, double high, double _intensity){
		spectrum = Spectrum(low, high); 
		intensity = _intensity;	
	}
	Photon(){
		spectrum = Spectrum();
		intensity = 0;
	}
	Photon(Spectrum &s, double &d){
		spectrum = s;
		intensity = d;
	}
	RGB rgb(){
		return toRGB(spectrum)*intensity;
	}
	Photon operator+(Photon &q){
		double pint = spectrum.integral();
		if(pint == 0.0) return q;
		double qint = q.spectrum.integral();
		if(qint == 0.0) return Photon(spectrum, intensity);
		Spectrum spec  = spectrum*(intensity/pint) + q.spectrum*(q.intensity/qint); 
		double ints  = intensity + q.intensity;
		return Photon( spec, ints );
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

double dot(Vector &n, Vector &m){
	return n.x*m.x + n.y*m.y + n.z*m.z;
}

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

struct Sphere{
	Point c;
	double r;
	Surface surface;
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
	Vector normal(Point &p){
		return Vector(c, p);	
	}
	Point intersection(Ray &l){
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
Sphere air;
Point operator+(const Point &p, const Vector &v){
	return Point(p.x + v.x, p.y + v.y, p.z + v.z);
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


struct Tracer{
	Point *light;
	Point *focal;
	Spectrum lightSpectrum;
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
		lightSpectrum = Spectrum(380,750);
	}
	void read(const char *const name){
		ifstream file;
		file.open(name);
		double x,y,z,rad, low, high, diff, spec, trans, difr;
		while (file >> x>>y>>z>>rad>>low>>high>>diff>>spec>>trans>>difr){
			objects.push_back(Sphere(Point(x,y,z),rad,Surface(Spectrum(low, high), diff, spec, trans, difr)));
		}
	}
	Photon traceRay(Ray &l, int depth, int currObject){
		Photon returnLight(380., 750., 0.);
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
			Ray shadow(minPoint, *light);
			int indShadow = -1;
			double minDistS = INF;
			Point minPointS(INF, INF, INF);

			for (int o = 0; o < objects.size(); ++o){
				Point intersect = objects[o].intersection(shadow);
				double dist = intersect.dist(*focal);
				if (dist < minDistS){
					minDistS = dist;
					indShadow = o;
					minPointS = intersect;
					break;
				}
			}
			int nextObject = indSource;
			if(currObject == indSource){
				nextObject = -1;
			}
			Vector viewVector = l.v;
			viewVector.normalize();

			Vector surfaceNormal = objects[indSource].normal(minPoint);
			surfaceNormal.normalize();
			if(nextObject == -1){
				surfaceNormal = -1 * surfaceNormal;
			}
			
			Vector lightVector = shadow.v;
			lightVector.normalize();
			
			if (depth < MAXDEPTH && objects[indSource].surface.specular > 0.0){
				Vector backwardReflection = viewVector - 2.*dot(viewVector, surfaceNormal)*surfaceNormal;
				Ray reflectedRay(minPoint, backwardReflection);
				Photon reflectedPhoton = traceRay(reflectedRay, depth+1, currObject);
				reflectedPhoton.intensity *= objects[indSource].surface.specular;
				returnLight = reflectedPhoton;
			}
			
			if (depth < MAXDEPTH && objects[indSource].surface.transmissive > 0.0){
				double currRefraction = 1.;
				if (currObject != -1){
					currRefraction = objects[currObject].surface.refraction;
				}
				double nextRefraction = 1.;
				if(nextObject != -1){
					nextRefraction = objects[nextObject].surface.refraction;
				}
				double refractionVal = currRefraction/nextRefraction;
				Vector backwardView = -1.*viewVector;
				double cosi = dot(surfaceNormal, backwardView);
				Vector refractedVector = refractionVal*viewVector + (refractionVal*cosi - sqrt(1 + pow(refractionVal, 2) * (pow(cosi, 2) - 1)))*surfaceNormal;
				Ray refractedRay(minPoint, refractedVector);
				Photon refractedPhoton = traceRay(refractedRay, depth+1, nextObject);
				if(nextObject != -1){
					refractedPhoton.intensity *= objects[indSource].surface.transmissive;
				}
				returnLight = returnLight + refractedPhoton;
			}
			
			double diffuseIntensity = max(dot(lightVector, surfaceNormal), AMBIENT)*objects[indSource].surface.diffuse;
			if(indShadow == -1){
				Vector reflection = lightVector - 2.*dot(lightVector, surfaceNormal)*surfaceNormal;

				double specularIntensity = max(0., pow(dot(reflection, viewVector),10))*objects[indSource].surface.specular; 
				Photon shiny(lightSpectrum,specularIntensity);
				returnLight = returnLight + shiny;
				
			}else{	
				diffuseIntensity =  AMBIENT*objects[indSource].surface.diffuse;
			}
			Photon diffuse =  Photon(objects[indSource].surface.spectrum, diffuseIntensity);
			returnLight = returnLight + diffuse;
		}
		return returnLight;
	}


	void drawScene(){
		for (int j = 0; j < HEIGHT; ++j){
			for(int i = 0; i < WIDTH; ++i){
				double xoff = (i-WIDTH/2)*PIXSIZE;
				double yoff = (HEIGHT/2-j)*PIXSIZE;
				Point px = (*focal)+(*dir)+xoff*(*horz) + yoff*(*vert);
				Ray l(*focal, px);

				Photon p = traceRay(l, 1, -1);
				
				RGB c = p.rgb();

				(*img)(i,j,0) = c.r;
				(*img)(i,j,1) = c.g;
				(*img)(i,j,2) = c.b;

			}
		}
	}
};

int main() {
	red = colorSpec(610, 25);
	green = colorSpec(540, 30);
	blue = colorSpec(445, 25);
	air = Sphere(Point(0,0,0), 0.0, Surface(Spectrum(), 0.0, 0.0, 0.0, 1.0));
	Tracer t;
	t.read("in.txt");
	t.drawScene();
	t.img->save("out.bmp");
	return 0;
}
