#include "general_math.cuh"

bool isNan(double3 v)
{
	return std::isnan(v.x) || std::isnan(v.y) || std::isnan(v.z);
}

__host__ __device__ double powN(double x, int p)
{
	if(p < 0)
		return powN(1/x, -p);

	double res = 1;

	switch(p){
	case 0: return 1;

	case 1: return x;

	case 2: return pow2(x);

	case 3: return pow3(x);

	case 4: return pow2(pow2(x));

	case 5: return pow2(pow2(x))*x;

	case 6: return pow3(pow2(x));

	default:
		res = pow3(pow2(x)) * x; // p = 7
		for(int i = 7; i < p; ++i) res *= x;
		return res;
	}
}

__host__ __device__ bool isIntPow(int n, int p)
{
	return int(powN(int(pow(n, 1.0/p) + 0.5), p) + 0.5) == n;
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
__host__ __device__ double pow3(double a){ return pow2(a)*a; }
__host__ __device__ double pow5(double a){ return pow2(pow2(a))*a; }

double myRnd(void){ return rand()/double(RAND_MAX); }
double myRndEx(void){ return ((rand() % (RAND_MAX - 2)) + 1)/double(RAND_MAX); }
double myRnd(double a, double b){
	if(a > b) swap(a, b);
	return myRnd() * (b - a) + a;
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
__host__ __device__ double gaussFnc(double x, double sgm, double x0)
{
	double b = (x-x0)/sgm;
	return exp(-b*b/2) / (sqrt(2*M_PI)*sgm);
}
// this gauss is checked - it's really gauss
double gaussRand_my(double sgm, double x0, double rng)
{
	rng *= sgm;
	double x;
	double y0 = 1.0/(sqrt(2*M_PI)*sgm); //y0 = gaussFnc(x0, sgm, x0); // max value
	double xl = x0 - rng, xr = x0 + rng;

	do{
		x = myRnd(xl, xr);
	}while(myRnd(0, y0) > gaussFnc(x, sgm, x0));

	return x;
}
double gaussRand(double sgm, double x0)
{
	return x0 + sgm * sin(2 * M_PI * myRndEx()) * sqrt(-2 * log(myRndEx()));
}
double gauss3D(double T, double rng) // rng = 5 is ok
{
	// T = v0^2 = sgm^2
	rng *= sqrt(T);
	double C = sqrt(2 / M_PI) / pow(T, 1.5);
	double y0 = 2 / e_const * sqrt(2 / (M_PI * T));
	double v, v2;

	do{
		v = myRnd(0, rng);
		v2 = v * v;
	}while(myRnd(0, y0) > C * v2 * exp(-v2 / (2 * T)));

	return v;
}

__host__ __device__ bool almostEq(double a, double b, double _eps)
{
	//return (b == 0 ? std::abs(a) : std::abs(a/b-1)) < _eps;
	return epsDlt(a, b, _eps) < _eps;
}
__host__ __device__ double epsDlt(double a, double b, double _eps)
{
	return b == 0 ? (a == 0 ? 0 : 1/_eps) : (std::abs(a) > std::abs(b) ? (a/b-1) : (b/a-1));
}

// -------------------------- shift ---------------------------------
__host__ __device__ double shiftX(double x, double R)
{
    if(x > R){ return x - 2 * R; }
    else if(x < -R) { return x + 2 * R; }
    return x;
}

__host__ __device__ double shiftXtrue(double x, double R)
{
	return shiftXtrue(x, R, 2*R);
}

__host__ __device__ double shiftX(double x, double R, double D)
{
    if(x > R){ return x - D; }
    else if(x < -R) { return x + D; }
    return x;
}

__host__ __device__ double shiftXtrue(double x, double R, double D)
{
	return x = x - floor(x/D + 0.5)*D; // X in [-R;R]
    // return (x - floor(x/D + 0.5)*D); // X in [-R;R] // TODO - try this variant
}

__host__ __device__ double3 shiftR(double3 r, double R)
{
    return shiftR(r, R, 2 * R);
}

__host__ __device__ double3 shiftRtrue(double3 r, double R)
{
    return shiftRtrue(r, R, 2 * R);
}

__host__ __device__ double3 shiftR(double3 r, double R, double D)
{
	r.x = shiftX(r.x, R, D);
	r.y = shiftX(r.y, R, D);
	r.z = shiftX(r.z, R, D);
    return r;
}

__host__ __device__ double3 shiftRtrue(double3 r, double R, double D)
{
	r.x = shiftXtrue(r.x, R, D);
	r.y = shiftXtrue(r.y, R, D);
	r.z = shiftXtrue(r.z, R, D);
    return r;
}

__host__ __device__ double getForce(double r2)
{
	r2 = 1/r2;
	double ar6 = pow3(r2);
	return 48.0 * ar6 * r2 * (0.5 - ar6);
}

__host__ __device__ double getEp(double r2)
{
	r2 = 1/r2;
	double ar6 = pow3(r2);
	return 4 * ar6 * (ar6 - 1);
}
