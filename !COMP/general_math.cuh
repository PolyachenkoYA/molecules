#pragma once

#include "const.cuh"

// ----------------------------------- sgn -------------------------------------
__host__ __device__ int sgn(double x, double _eps = SYS_EPS);
__host__ __device__ double epsDlt(double a, double b, double _eps = SYS_EPS);
__host__ __device__ double float_part(double x);
__host__ __device__ bool almostEq(double a, double b, double _eps = SYS_EPS);
__host__ __device__ bool is0(double3 v);
__host__ __device__ char bool2int(bool b);
bool isNan(double3 v);

// --------------------------------------- rnd -----------------------------------
double myRnd(void);
double myRnd(double a, double b);
__host__ __device__ double3 vecByAngles(double phi, double tht = 0);
double3 rndVec(double V = 1);
__host__ __device__ double gaussFnc(double x, double sgm = 1, double x0 = 0);
double gaussRand(double sgm = 1, double x0 = 0, double rng = 5);

// -------------------------------------- common fncs -----------------------------
__host__ __device__ double pow2(double a);
__host__ __device__ double pow3(double a);
__host__ __device__ double pow5(double a);
__host__ __device__ double shiftX(double x, double R);
__host__ __device__ double shiftX(double x, double R, double D);
__host__ __device__ double3 shiftR(double3 r, double R);
__host__ __device__ double3 shiftR(double3 r, double R, double D);
__host__ __device__ double cos_sin(double x);


