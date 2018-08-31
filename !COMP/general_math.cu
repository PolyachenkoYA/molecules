#include "general_math.cuh"

bool isNan(double3 v)
{
	return std::isnan(v.x) || std::isnan(v.y) || std::isnan(v.z);
}

__host__ __device__ char bool2int(bool b){ return b ? 1 : 0; }

__host__ __device__ bool is0(double3 v){ return (v.x == 0) && (v.y == 0) && (v.z == 0); }

__host__ __device__ double float_part(double x)
{
    return x - int(x);
}

__host__ __device__ double cos_sin(double x)
{
    return sqrtf(1-x*x);
}

__host__ __device__ int sgn(double x, double _eps)
{
    if(x >= _eps){
    	return 1;
    } else if(x <= -_eps){
        return -1;
    } else {
        return 0;
    }
}

__host__ __device__ double pow2(double a){ return a*a; }
__host__ __device__ double pow3(double a){ return a*a*a; }
__host__ __device__ double pow5(double a){ return a*a*a*a*a; }

double myRnd(void){ return rand()/double(RAND_MAX); }
double myRnd(double a, double b){
	if(a>b) swap(a,b);
	return myRnd()*(b-a)+a;
}
__host__ __device__ double3 vecByAngles(double phi, double tht )
{
	double ct = cos(tht);
	return make_double3(ct * cos(phi), ct * sin(phi), sin(tht));
}
double3 rndVec(double V)
{
	return vecByAngles(myRnd(-pi, pi), myRnd(-pi_d2, pi_d2))*V;
}
// this gauss is checked - it's really gauss
__host__ __device__ double gaussFnc(double x, double sgm, double x0)
{
	double b = (x-x0)/sgm;
	return exp(-b*b/2) / (sqrt(2*M_PI)*sgm);
}
double gaussRand(double sgm, double x0, double rng)
{
	rng *= sgm;
	double x;
	double y0 = 1.0/(sqrt(2*M_PI)*sgm); //y0 = gaussFnc(x0, sgm, x0); // max value
	double xl = x0-rng, xr = x0+rng;

	do{
		x = myRnd(xl, xr);
	}while(myRnd(0, y0) > gaussFnc(x, sgm, x0));

	return x;
}

__host__ __device__ bool almostEq(double a, double b, double _eps)
{
	return (b == 0 ? std::abs(a) : std::abs(a/b-1)) < _eps;
}
__host__ __device__ double epsDlt(double a, double b, double _eps)
{
	return b == 0 ? (a == 0 ? 0 : 1/_eps) : (std::abs(a) > std::abs(b) ? (a/b-1) : (b/a-1));
}

// -------------------------- shift ---------------------------------
__host__ __device__ double shiftX(double x, double R)
{
    if(x > R){ return x - 2*R; }
    else if(x < -R) { return x + 2*R; }
    return x;
}

__host__ __device__ double shiftX(double x, double R, double D)
{
    if(x > R){ return x - D; }
    else if(x < -R) { return x + D; }
    return x;
}

__host__ __device__ double3 shiftR(double3 r, double R)
{
    return shiftR(r, R, 2*R);
}

__host__ __device__ double3 shiftR(double3 r, double R, double D)
{
	r.x = shiftX(r.x, R, D);
	r.y = shiftX(r.y, R, D);
	r.z = shiftX(r.z, R, D);
    return r;
}
