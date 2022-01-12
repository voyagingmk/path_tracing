#include <vector>
#include <iostream>
#include <type_traits>
#include <random>
#include <math.h>   // smallpt, a Path Tracer by Kevin Beason, 2008 
#include <stdlib.h> // Make : g++ -O3 -fopenmp smallpt.cpp -o smallpt 
#include <stdio.h>  //        Remove "-fopenmp" for g++ version < 4.2 
#include <omp.h>
#include <assimp/Importer.hpp>     // C++ importer interface
#include <assimp/Exporter.hpp>     // C++ importer interface
#include <assimp/scene.h>          // Output data structure
#include <assimp/postprocess.h>     // Post processing flags

#include "bvh/GRIDQBVH.h"


//#define DEBUG 1

#ifdef DEBUG
#define LOG printf
#else
#define LOG
#endif


#if 0

# define M_PI           3.14159265358979323846
# define M_1_PI			(1.0/3.14159265358979323846)

std::default_random_engine generator;
std::uniform_real_distribution<double> distr(0.0, 1.0);
double erand48(unsigned short * X) {
	return distr(generator);
}

struct Vec {        // Usage: time ./explicit 16 && xv image.ppm
	double x, y, z;                  // position, also color (r,g,b)
	Vec(double x_ = 0, double y_ = 0, double z_ = 0) { x = x_; y = y_; z = z_; }
	Vec operator+(const Vec &b) const { return Vec(x + b.x, y + b.y, z + b.z); }
	Vec operator-(const Vec &b) const { return Vec(x - b.x, y - b.y, z - b.z); }
	Vec operator*(double b) const { return Vec(x*b, y*b, z*b); }
	Vec mult(const Vec &b) const { return Vec(x*b.x, y*b.y, z*b.z); }
	Vec& norm() { return *this = *this * (1 / sqrt(x*x + y * y + z * z)); }
	double dot(const Vec &b) const { return x * b.x + y * b.y + z * b.z; } // cross:
	Vec operator%(Vec&b) { return Vec(y*b.z - z * b.y, z*b.x - x * b.z, x*b.y - y * b.x); }
};
struct Ray { Vec o, d; Ray(Vec o_, Vec d_) : o(o_), d(d_) {} };
enum Refl_t { DIFF, SPEC, REFR };  // material types, used in radiance()
struct Sphere {
	double rad;       // radius
	Vec p, e, c;      // position, emission, color
	Refl_t refl;      // reflection type (DIFFuse, SPECular, REFRactive)
	Sphere(double rad_, Vec p_, Vec e_, Vec c_, Refl_t refl_) :
		rad(rad_), p(p_), e(e_), c(c_), refl(refl_) {}
	double intersect(const Ray &r) const { // returns distance, 0 if nohit
		Vec op = p - r.o; // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
		double t, eps = 1e-4, b = op.dot(r.d), det = b * b - op.dot(op) + rad * rad;
		if (det < 0) return 0; else det = sqrt(det);
		return (t = b - det) > eps ? t : ((t = b + det) > eps ? t : 0);
	}
};
Sphere spheres[] = {//Scene: radius, position, emission, color, material
  Sphere(1e5, Vec(1e5 + 1,40.8,81.6), Vec(),Vec(.75,.25,.25),DIFF),//Left
  Sphere(1e5, Vec(-1e5 + 99,40.8,81.6),Vec(),Vec(.25,.25,.75),DIFF),//Rght
  Sphere(1e5, Vec(50,40.8, 1e5),     Vec(),Vec(.75,.75,.75),DIFF),//Back
  //Sphere(1e5, Vec(50,40.8,-1e5 + 170), Vec(),Vec(),           DIFF),//Frnt
  Sphere(1e5, Vec(50, 1e5, 81.6),    Vec(),Vec(.75,.75,.75),DIFF),//Botm
  Sphere(1e5, Vec(50,-1e5 + 81.6,81.6),Vec(),Vec(.75,.75,.75),DIFF),//Top
  Sphere(16.5,Vec(27,16.5,47),       Vec(),Vec(1,1,1)*.999, SPEC),//Mirr
  Sphere(16.5,Vec(73,16.5,78),       Vec(),Vec(1,1,1)*.999, REFR),//Glas
  Sphere(1.5, Vec(50,81.6 - 16.5,81.6),Vec(4,4,4) * 100,  Vec(), DIFF),//Lite
};
int numSpheres = sizeof(spheres) / sizeof(Sphere);
inline double clamp(double x) { return x < 0 ? 0 : x>1 ? 1 : x; }
inline int toInt(double x) { return int(pow(clamp(x), 1 / 2.2) * 255 + .5); }
inline bool intersect(const Ray &r, double &t, int &id) {
	double n = sizeof(spheres) / sizeof(Sphere), d, inf = t = 1e20;
	for (int i = int(n); i--;) if ((d = spheres[i].intersect(r)) && d < t) { t = d; id = i; }
	return t < inf;
}
Vec radiance(const Ray &r, int depth, unsigned short *Xi, int E = 1) {
	double t;                               // distance to intersection
	int id = 0;                               // id of intersected object
	if (!intersect(r, t, id)) return Vec(); // if miss, return black
	const Sphere &obj = spheres[id];        // the hit object
	Vec x = r.o + r.d*t, n = (x - obj.p).norm(), nl = n.dot(r.d) < 0 ? n : n * -1, f = obj.c;
	double p = f.x > f.y && f.x > f.z ? f.x : f.y > f.z ? f.y : f.z; // max refl
	if (++depth > 5 || !p) if (erand48(Xi) < p) f = f * (1 / p); else return obj.e*E;
	if (obj.refl == DIFF) {                  // Ideal DIFFUSE reflection
		double r1 = 2 * M_PI*erand48(Xi), r2 = erand48(Xi), r2s = sqrt(r2);
		Vec w = nl, u = ((fabs(w.x) > .1 ? Vec(0, 1) : Vec(1)) % w).norm(), v = w % u;
		Vec d = (u*cos(r1)*r2s + v * sin(r1)*r2s + w * sqrt(1 - r2)).norm();

		// Loop over any lights
		Vec e;
		for (int i = 0; i < numSpheres; i++) {
			const Sphere &s = spheres[i];
			if (s.e.x <= 0 && s.e.y <= 0 && s.e.z <= 0) continue; // skip non-lights

			Vec sw = s.p - x, su = ((fabs(sw.x) > .1 ? Vec(0, 1) : Vec(1)) % sw).norm(), sv = sw % su;
			double cos_a_max = sqrt(1 - s.rad*s.rad / (x - s.p).dot(x - s.p));
			double eps1 = erand48(Xi), eps2 = erand48(Xi);
			double cos_a = 1 - eps1 + eps1 * cos_a_max;
			double sin_a = sqrt(1 - cos_a * cos_a);
			double phi = 2 * M_PI*eps2;
			Vec l = su * cos(phi)*sin_a + sv * sin(phi)*sin_a + sw * cos_a;
			l.norm();
			if (intersect(Ray(x, l), t, id)) {
				if (id == i) {  // shadow ray
					double omega = 2 * M_PI * (1 - cos_a_max);
					e = e + f.mult(s.e * l.dot(nl) * omega) * M_1_PI;  // 1/pi for brdf
				}
			}
			else {
				printf("light ray not hit\n");
			}
		}

		return obj.e*E + e + f.mult(radiance(Ray(x, d), depth, Xi, 0));
	}
	else if (obj.refl == SPEC)              // Ideal SPECULAR reflection
		return obj.e + f.mult(radiance(Ray(x, r.d - n * 2 * n.dot(r.d)), depth, Xi));
	Ray reflRay(x, r.d - n * 2 * n.dot(r.d));     // Ideal dielectric REFRACTION
	bool into = n.dot(nl) > 0;                // Ray from outside going in?
	double nc = 1, nt = 1.5, nnt = into ? nc / nt : nt / nc, ddn = r.d.dot(nl), cos2t;
	if ((cos2t = 1 - nnt * nnt*(1 - ddn * ddn)) < 0)    // Total internal reflection
		return obj.e + f.mult(radiance(reflRay, depth, Xi));
	Vec tdir = (r.d*nnt - n * ((into ? 1 : -1)*(ddn*nnt + sqrt(cos2t)))).norm();
	double a = nt - nc, b = nt + nc, R0 = a * a / (b*b), c = 1 - (into ? -ddn : tdir.dot(n));
	double Re = R0 + (1 - R0)*c*c*c*c*c, Tr = 1 - Re, P = .25 + .5*Re, RP = Re / P, TP = Tr / (1 - P);
	return obj.e + f.mult(depth > 2 ? (erand48(Xi) < P ?   // Russian roulette
		radiance(reflRay, depth, Xi)*RP : radiance(Ray(x, tdir), depth, Xi)*TP) :
		radiance(reflRay, depth, Xi)*Re + radiance(Ray(x, tdir), depth, Xi)*Tr);
}
int main(int argc, char *argv[]) {
	int w = 512, h = 384, samps = 32; // # samples 
	Ray cam(Vec(50, 52, 295.6), Vec(0, -0.042612, -1).norm()); // cam pos, dir
	float fov = .5135;
	Vec cx = Vec(w* fov / h), cy = (cx%cam.d).norm()* fov, r, *c = new Vec[w*h];
#pragma omp parallel for schedule(dynamic, 1) private(r)       // OpenMP
	for (int y = 0; y < h; y++) {                       // Loop over image rows
		fprintf(stderr, "\rRendering (%d spp) %5.2f%%", samps * 4, 100. * y / (h - 1));
		for (unsigned short x = 0, Xi[3] = { 0,0,y * y * y }; x < w; x++)   // Loop cols
#if 1
		{
			int i = (h - y - 1) * w + x;
			r = Vec();
			for (int s = 0; s < samps; s++) {
				Vec d = cx * (x / float(w) - .5) + cy * (y / float(h) - .5) + cam.d;
				r = r + radiance(Ray(cam.o + d * 0, d.norm()), 0, Xi) * (1. / samps);
			}
			c[i] = c[i] + Vec(clamp(r.x), clamp(r.y), clamp(r.z));
		}
#else
			for (int sy = 0, i = (h - y - 1)*w + x; sy < 2; sy++)     // 2x2 subpixel rows
				for (int sx = 0; sx < 2; sx++, r = Vec()) {        // 2x2 subpixel cols
					for (int s = 0; s < samps; s++) {
						double r1 = 2 * erand48(Xi), dx = r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);
						double r2 = 2 * erand48(Xi), dy = r2 < 1 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);
						dx = 0.0;
						dy = 0.0;
						Vec d = cx * (((sx + .5 + dx) / 2 + x) / w - .5) +
							cy * (((sy + .5 + dy) / 2 + y) / h - .5) + cam.d;
						r = r + radiance(Ray(cam.o + d * 140, d.norm()), 0, Xi)*(1. / samps);
					} // Camera rays are pushed ^^^^^ forward to start in interior
					c[i] = c[i] + Vec(clamp(r.x), clamp(r.y), clamp(r.z))*.25;
				}
#endif
	}
	FILE *f;
	errno_t err = fopen_s(&f, "image2.ppm", "w");
	//FILE *f = fopen("image.ppm", "w");         // Write image to PPM file. 
	fprintf(f, "P3\n%d %d\n%d\n", w, h, 255);
	for (int i = 0; i < w*h; i++)
		fprintf(f, "%d %d %d ", toInt(c[i].x), toInt(c[i].y), toInt(c[i].z));
	return 0;
}

#else
using namespace accel;

Vector3 origin(278, 273, -1000);
Vector3 cam_direction(0, 0, 1);
Vector3 light_origin(300, 400, 300);
float radius = 25.0f;
s32 width = 640;
s32 height = 480;
const int samples = 10;
Vector3 cameraPosition;
Vector3 cameraForward;
Vector3 cameraUp;
Vector3 cameraRight;
f32 iw;
f32 ih;

Vector3 sx;
Vector3 sy;

RGB* rgb = NULL;

void initScene(const Vector3 position, const Vector3& direction)
{
	cameraPosition = position;
	cameraForward = normalize(direction);
	cameraUp = Vector3(0.0f, 1.0f, 0.0f);
	cameraRight = normalize(cross(cameraUp, cameraForward));
	f32 aspect = static_cast<f32>(width) / height;

	f32 fovy = ::tanf(0.5f * 60.0f / 180.0f * static_cast<f32>(M_PI));
	f32 dx = fovy * aspect;
	f32 dy = fovy;

	iw = 1.0f / (width - 1);
	ih = 1.0f / (height - 1);
	sx = dx * cameraRight;
	sy = dy * cameraUp;

	rgb = (RGB*)ACC_MALLOC(sizeof(RGB) * width * height);
}

void termScene()
{
	ACC_FREE(rgb);
}
# define M_PI           3.14159265358979323846
# define M_1_PI			(1.0/3.14159265358979323846)

std::default_random_engine generator;
const std::uniform_real_distribution<double> distr(0.0, 1.0);

std::random_device rdevice;
std::vector<std::default_random_engine> generators;
double erand48(unsigned short* X = nullptr) {
	return distr(generators[omp_get_thread_num()]);
}
enum Refl_t { DIFF, SPEC, REFR };  // material types, used in radiance()


class MyFace : public Face
{
public:
	MyFace() {}
	Vector3 getNormalByHitPoint(const Vector3& hitPoint) const {
		if (0&&hasNormal_) {
			float u = 0, v = 0, w = 1;
			Barycentric(hitPoint, u, v, w);
			Vector3 n = normal_[0] * u + normal_[1] * v + normal_[2] * w;
			n.norm();
			Vector3 n2 = getNormal();
			return n;
		}
		return getNormal();
	}
	void Barycentric(const Vector3& hitPoint, float& u, float& v, float& w) const
	{
		/*
		P_hit = p_1 + u * (p_2 - p_1) + v * (p_3 - p_1)
		u * (p_2 - p_1) + v * (p_3 - p_1) = P_hit - p_1
		B = p_2 - p_1
		C = p_3 - p_1
		A = P_hit - p_1
		=>
		u * B + v * C = A   £¨it is vector£¡£©
		convert to linear system =>
		(u * B + v * C)¡¤B = A¡¤B		£¨it is scalar£¡£©
		(u * B + v * C)¡¤C = A¡¤C

		u * B¡¤B  +  v * C¡¤B = A¡¤B
		u * B¡¤C  +  v * C¡¤C = A¡¤C

		¡¾ B¡¤B  C¡¤B ¡¿ ¡¾ u ¡¿   A¡¤B
		¡¾ B¡¤C  C¡¤C ¡¿ ¡¾ v ¡¿ = A¡¤C

		¡¾ d00  d01 ¡¿ ¡¾ u ¡¿   d20
		¡¾ d01  d11 ¡¿ ¡¾ v ¡¿ = d21

		det   = | d00  d01 |
			    | d01  d11 |

		det_u = | d20  d01 |
			    | d21  d11 |

		det_v = | d00  d20 |
			    | d01  d21 |
		*/
		Vector3 B = point_[1] - point_[0], C = point_[2] - point_[0], A = hitPoint - point_[0];
		float d00 = B.dot(B);
		float d01 = B.dot(C);
		float d11 = C.dot(C);
		float d20 = A.dot(B);
		float d21 = A.dot(C);
		float denom = d00 * d11 - d01 * d01;
		float denom_u = d20 * d11 - d01 * d21;
		float denom_v = d00 * d21 - d20 * d01;
		float inv_denom = 1 / denom;
		u = denom_u * inv_denom;
		v = denom_v * inv_denom;
		w = 1.0 - u - v;
	}
	Vector3 uv_[2];
	int id = 0;
	Vector3 e, c;      // emission, color
	Refl_t refl;      // reflection type (DIFFuse, SPECular, REFRactive)
};

//

Vector3 radiance(GRIDQBVH<MyFace>* bvh, accel::Ray& r, int depth, unsigned short* Xi, int E = 1) {
	                            // distance to intersection
	++depth;
	HitRecord hitRecord = bvh->intersect(r);
	if (!hitRecord.primitive_) {
		LOG("depth %d miss\n", depth);
		return Vector3(0, 0, 0);
	}
	const MyFace* face = reinterpret_cast<const MyFace*>(hitRecord.primitive_);                        // distance to intersection
	int id = face->id;
	double t = hitRecord.t_;
	if (id == 6 && depth == 1)
	{
		id = id + 0;

	}
	Vector3 x = r.origin_ + r.direction_ * t;
	Vector3 n = face->getNormalByHitPoint(x);
	Vector3 nl = n.dot(r.direction_) < 0 ? n : n * -1;
	Vector3 f = face->c;
	double p = f.x() > f.y() && f.x() > f.z() ? f.x() : f.y() > f.z() ? f.y() : f.z(); // max refl
	if (depth > 5 || !p) if (erand48(Xi) < p) f = f * (1 / p); else return face->e * E;
	if (face->refl == DIFF) {                  // Ideal DIFFUSE reflection
		double r1 = 2 * M_PI * erand48(Xi), r2 = erand48(Xi), r2s = sqrt(r2);
		Vector3 w = nl, u = ((fabs(w.x()) > .1 ? Vector3(0, 1, 0) : Vector3(1, 0, 0)).cross(w)).norm(), v = w.cross(u);
		Vector3 d = (u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2)).norm();

		Vector3 e(0, 0, 0);
		
		// Loop over any lights
		double srad = radius;
		Vector3 sp = light_origin;
		Vector3 sw = sp - x, su = ((fabs(sw.x_) > .1 ? Vector3(0, 1, 0) : Vector3(1,0,0)).cross(sw)).norm(), sv = sw.cross(su);
		//sw.norm();
		double cos_a_max = sqrt(1 - srad * srad / (x - sp).dot(x - sp));
		double eps1 = erand48(Xi), eps2 = erand48(Xi);
		double cos_a = 1 - eps1 + eps1 * cos_a_max;
		double sin_a = sqrt(1 - cos_a * cos_a);
		double phi = 2 * M_PI * eps2;
		Vector3 l = su * cos(phi) * sin_a + sv * sin(phi) * sin_a + sw * cos_a;
		l.norm();
		double omega = 2 * M_PI * (1 - cos_a_max);
		HitRecord hitRecord2 = bvh->intersect(accel::Ray(x, l));
		if (hitRecord2.primitive_) {
			const MyFace* face2 = reinterpret_cast<const MyFace*>(hitRecord2.primitive_);
			if (face2->id == 0) {
				e = e + f.mult(face2->e * l.dot(nl) * omega) * M_1_PI;  // 1/pi for brdf
				LOG("depth %d hit light e %.5f %.5f %.5f\n", depth, e.x_, e.y_, e.z_);
			}
			else {
				LOG("depth %d blocked %d\n", depth, face2->id);
			}
		}
		else {
			LOG("light ray not hit\n");
		}
		Vector3 bounce_radiance = radiance(bvh, accel::Ray(x, d), depth, Xi, 0);
		LOG("depth %d bounce_radiance %.5f %.5f %.5f\n", depth, bounce_radiance.x_, bounce_radiance.y_, bounce_radiance.z_);
		return face->e * E + e + f.mult(bounce_radiance);
	}
	else if (face->refl == SPEC) {           // Ideal SPECULAR reflection
		double n_dot_r = n.dot(r.direction_);
		Vector3 n2 = n * 2.0;
		Vector3 n2_n_dot_r = n2 * n_dot_r;
		Vector3 r2 = r.direction_ - n2_n_dot_r;
		Vector3 r3 = r.direction_ - n * 2.0 * n.dot(r.direction_);

		return face->e + f.mult(radiance(bvh, accel::Ray(x, r2), depth, Xi));
		return face->e + f.mult(radiance(bvh, accel::Ray(x, r.direction_ - n * 2.0 * n.dot(r.direction_)), depth, Xi));
	}
	accel::Ray reflRay(x, r.direction_ - n * 2 * n.dot(r.direction_));     // Ideal dielectric REFRACTION
	bool into = n.dot(nl) > 0;                // Ray from outside going in?
	double nc = 1, nt = 1.5, nnt = into ? nc / nt : nt / nc, ddn = r.direction_.dot(nl), cos2t;
	if ((cos2t = 1 - nnt * nnt * (1 - ddn * ddn)) < 0)    // Total internal reflection
		return face->e + f.mult(radiance(bvh, reflRay, depth, Xi));
	Vector3 tdir = (r.direction_ * nnt - n * ((into ? 1 : -1) * (ddn * nnt + sqrt(cos2t)))).norm();
	double a = nt - nc, b = nt + nc, R0 = a * a / (b * b), c = 1 - (into ? -ddn : tdir.dot(n));
	double Re = R0 + (1 - R0) * c * c * c * c * c, Tr = 1 - Re, P = .25 + .5 * Re, RP = Re / P, TP = Tr / (1 - P);
	return face->e + f.mult(depth > 2 ? (erand48(Xi) < P ?   // Russian roulette
		radiance(bvh, reflRay, depth, Xi) * RP : radiance(bvh, accel::Ray(x, tdir), depth, Xi) * TP) :
		radiance(bvh, reflRay, depth, Xi) * Re + radiance(bvh, accel::Ray(x, tdir), depth, Xi) * Tr);
}

inline double clamp(double x) { return x < 0 ? 0 : x>1 ? 1 : x; }
inline int toInt(double x) { return int(pow(clamp(x), 1 / 2.2) * 255 + .5); }

template<class T>
void test(
	s32& depth,
	f64& buildTime,
	f64& renderTime,
	f64& raysPerSecond,
	s32 numPrimitives,
	MyFace* faces,
	const char* filename)
{
	Timer<true> timer;
	T* bvh = new T;

	timer.start();
	bvh->build(numPrimitives, faces);
	timer.stop();
	buildTime = timer.getAverage();

	depth = bvh->getDepth();
	timer.reset();

	timer.start();
	for (s32 i = 0; i < height; ++i) {
		fprintf(stderr, "\rRendering  %5.2f%%",  100. * i / (height - 1));
#ifndef _DEBUG
#pragma omp parallel for schedule(dynamic, 1)
#endif
		for (s32 j = 0; j < width; ++j) {
			RGB& pixel = rgb[i * width + j];
			pixel.r_ = pixel.g_ = pixel.b_ = 0;
			for (int suby = 0; suby < 2; suby++)     // 2x2 subpixel rows
				for (int subx = 0; subx < 2; subx++) {        // 2x2 subpixel cols
					Vector3 r(0, 0, 0);
					for (int s = 0; s < samples; s++) {

						///Vector3 d = sx * (x / float(w) - .5) + sy * (y / float(h) - .5) + cam.d;
						//Vector3 vx = sx * (2.0f * j  / (width - 1) - 1.0f);
						//Vector3 vy = sy * (1.0f - 2.0f * i / (height - 1));
						//Vec d = cx * (((sx + .5 + dx) / 2 + x) / w - .5) +
						//	cy * (((sy + .5 + dy) / 2 + y) / h - .5) + cam.d;
						double r1 = 2 * erand48(), dx = r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);
						double r2 = 2 * erand48(), dy = r2 < 1 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);
						//Vector3 vx_ = sx * (j / (double)width - 0.5f);
						//Vector3 vy_ = sy * ((height - i - s32(1)) / (double)height - 0.5f);
						Vector3 d = sx * (((subx + .5 + dx) / 2 + j) / width - .5) +
							sy * (((suby + .5 + dy) / 2 + height - i - 1) / height - .5) + cameraForward;
						Vector3 r_ = radiance(bvh, Ray(cameraPosition, d), 0, nullptr);
						//Vector3 r_ = radiance(bvh, Ray(cameraPosition, cameraForward + vx_ + vy_), 0, nullptr);
						//Vector3 r_ = radiance(bvh, Ray(cameraPosition, normalize(cameraForward + vx + vy)), 0, nullptr);
						r = r + r_ * (1.0 / double(samples));
					}
					pixel.r_ += r.x() * .25;
					pixel.g_ += r.y() * .25;
					pixel.b_ += r.z() * .25;
				}
		}
	}
	timer.stop();
	delete bvh;
	renderTime = timer.getAverage();
	raysPerSecond = (height * width) / renderTime;
	printImage(filename, rgb, width, height);
	printf("%s is done.\n", filename);
}

template<class T>
void test(
	FILE* file,
	s32 numPrimitives,
	MyFace* faces,
	const char* name,
	const char* filename)
{
	s32 depth = 0;
	f64 buildTime = 0.0;
	f64 renderTime = 0.0;
	f64 raysPerSecond = 0.0;
	test<T>(depth, buildTime, renderTime, raysPerSecond, numPrimitives, faces, filename);
	fprintf(file, "%s depth: %d, build: %f, render: %f, rays: %f\r\n", name, depth, buildTime, renderTime, raysPerSecond);
}


Assimp::Importer importer;
const aiScene* loadAssimp(const char* filepath) {
	// Create an instance of the Importer class
	importer.SetPropertyInteger(AI_CONFIG_PP_RVC_FLAGS, aiComponent_TEXTURES);
	return importer.ReadFile(filepath,
		aiProcess_RemoveComponent |
		aiProcess_OptimizeMeshes |
		aiProcess_FindInvalidData |
		aiProcess_Triangulate |
		aiProcess_ValidateDataStructure |
		aiProcess_JoinIdenticalVertices |
		aiProcess_SortByPType);
}

int main()
{
	for (int i = 0, N = omp_get_max_threads(); i < N; ++i) {
		generators.emplace_back(std::default_random_engine(rdevice()));
	}
	//initScene(Vector3(300.0f, 300.0f, 300.0f), Vector3(-1.0f, -1.0f, -1.0f));
	initScene(origin, cam_direction);
	//const aiScene* scene = loadAssimp("./zhuzi_bihe.fbx");
	struct Config
	{
		Refl_t refl;
		const char* path;
		Vector3 emission;
		Vector3 albedo;
		Vector3 offset = Vector3(0,0,0);
		bool normalize = false;
		double scale = 1.0f;
		int mNumFaces = 0;
		MyFace* faces = nullptr;
	};
	Config config[] = {
		{DIFF, "./sphere.fbx", Vector3(100, 100, 100),Vector3(0, 0, 0), light_origin, true, radius},
		{DIFF, "./cbox/meshes/cbox_back.obj", Vector3(0, 0, 0),Vector3(0.75, 0.75, 0.75)},
		{DIFF, "./cbox/meshes/cbox_floor.obj", Vector3(0, 0, 0),Vector3(0.75, 0.75, 0.75)},
		{DIFF, "./cbox/meshes/cbox_greenwall.obj", Vector3(0, 0, 0),Vector3(0.105421, 0.37798, 0.076425)},
		{DIFF, "./cbox/meshes/cbox_redwall.obj", Vector3(0, 0.0, 0),Vector3(0.570068, 0.0430135, 0.0443706)},
		{DIFF, "./cbox/meshes/cbox_ceiling.obj", Vector3(0, 0, 0),Vector3(0.75, 0.75, 0.75)},
		{SPEC, "./cbox/meshes/cbox_smallbox.obj", Vector3(0, 0, 0),Vector3(0.99, 0.99, 0.99)},
		{DIFF, "./cbox/meshes/cbox_largebox.obj", Vector3(0, 0, 0.0),Vector3(0.4, 0.6, 0.6)},
		{SPEC, "./sphere.fbx",  Vector3(0, 0, 0), Vector3(0.99, 0.99, 0.99), Vector3(400, 150,150), true, 70.0f},

		//{DIFF, "./sphere.fbx", Vector3(1, 1, 1),Vector3(0, 0, 0)},
		//{DIFF, "./cbox/meshes/cbox_back.obj", Vector3(1, 1, 1),Vector3(0.0, 0.0, 0.0)},
		//{DIFF, "./cbox/meshes/cbox_floor.obj", Vector3(1, 1, 1),Vector3(0.0, 0.0, 0.0)},
		//{DIFF, "./cbox/meshes/cbox_greenwall.obj", Vector3(0, 0, 1),Vector3(0, 0, 0)},
		//{DIFF, "./cbox/meshes/cbox_redwall.obj", Vector3(0, 1, 0),Vector3(0, 0, 0)},
		//{DIFF, "./cbox/meshes/cbox_ceiling.obj", Vector3(1, 1, 1),Vector3(0, 0, 0)},
		//{SPEC, "./cbox/meshes/cbox_smallbox.obj", Vector3(0, 0.9, 0.9),Vector3(0, 0, 0)},
		//{DIFF, "./cbox/meshes/cbox_largebox.obj", Vector3(0.9, 0, 0.9),Vector3(0, 0, 0)},
	};
	int count = sizeof(config) / sizeof(Config);
	int facesTotal = 0;
	for (int idx = 0; idx < count; idx++)
	{
		const aiScene* scene = loadAssimp(config[idx].path);
		aiMesh* mesh = scene->mMeshes[0];
		int mNumFaces = mesh->mNumFaces;
		config[idx].mNumFaces = mNumFaces;
		MyFace* &faces = config[idx].faces;
		faces = new MyFace[mNumFaces];
		facesTotal += mNumFaces;

		for (int i = 0; i < mNumFaces; ++i) {
			aiFace& f = mesh->mFaces[i];
			s32 idx0 = f.mIndices[0];
			s32 idx1 = f.mIndices[1];
			s32 idx2 = f.mIndices[2];
			aiVector3D& v1 = mesh->mVertices[idx0];
			aiVector3D& v2 = mesh->mVertices[idx1];
			aiVector3D& v3 = mesh->mVertices[idx2];
			faces[i].id = idx;
			faces[i].e = config[idx].emission;
			faces[i].c = config[idx].albedo;
			faces[i].refl = config[idx].refl;
			faces[i].point_[0] = Vector3(v1.x, v1.y, v1.z);
			faces[i].point_[1] = Vector3(v2.x, v2.y, v2.z);
			faces[i].point_[2] = Vector3(v3.x, v3.y, v3.z);

			if (config[idx].normalize) {
				faces[i].point_[0].norm();
				faces[i].point_[1].norm();
				faces[i].point_[2].norm();
			}

			faces[i].point_[0] *= config[idx].scale;
			faces[i].point_[1] *= config[idx].scale;
			faces[i].point_[2] *= config[idx].scale;
			faces[i].point_[0] += config[idx].offset;
			faces[i].point_[1] += config[idx].offset;
			faces[i].point_[2] += config[idx].offset;
			int uv_count = mesh->GetNumUVChannels();
			if (mesh->GetNumUVChannels() > 0)
			{
				aiVector3D& uv1 = mesh->mTextureCoords[0][idx0];
				aiVector3D& uv2 = mesh->mTextureCoords[0][idx1];
				aiVector3D& uv3 = mesh->mTextureCoords[0][idx2];
				faces[i].uv_[0] = Vector3(uv1.x, uv1.y, 0);
				faces[i].uv_[1] = Vector3(uv2.x, uv2.y, 0);
				faces[i].uv_[2] = Vector3(uv3.x, uv3.y, 0);
			}

			if (mesh->mNormals)
			{
				aiVector3D& n1 = mesh->mNormals[idx0];
				aiVector3D& n2 = mesh->mNormals[idx1];
				aiVector3D& n3 = mesh->mNormals[idx2];
				faces[i].hasNormal_ = true;
				faces[i].normal_[0] = Vector3(n1.x, n1.y, n1.z);
				faces[i].normal_[1] = Vector3(n2.x, n2.y, n2.z);
				faces[i].normal_[2] = Vector3(n3.x, n3.y, n3.z);
			}
		} 

	}
	MyFace* faces = new MyFace[facesTotal];
	int n = 0;
	for (int i = 0; i < count; i++)
	{
		memcpy(faces + n, config[i].faces, sizeof(MyFace)* config[i].mNumFaces);
		n += config[i].mNumFaces;
	}

	FILE* file = NULL;
	fopen_s(&file, "statistics.txt", "wb");
	if (NULL == file) {
		return 0;
	}
	test<GRIDQBVH<MyFace> >(file, facesTotal, faces, "test", "test.ppm");
	termScene();
	return 0;
}


#endif