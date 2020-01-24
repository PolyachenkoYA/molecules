#pragma once

#include <math.h>
#include <cmath>
#include <vector>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <omp.h>
#include <algorithm>
#include <random>
#include <numeric>

#include <iostream>
#include <stdlib.h>
#include <ctime>
#include <time.h>
#include <fstream>
#include <cstring>
#include <sys/stat.h>
#include <string>
#include <sstream>
#include <unistd.h>
#include <new>
#include <iomanip>
#include <cuda_profiler_api.h>
#include <utility>

//#include <math_helper.h>
#include "math_helper.h"

using namespace std;

// ----------------------------------- math -----------------------------------------------------------

#define pi 3.14159265359
#define pi_d2 1.57079631679
#define e_const 2.71828182846
#define SQRT_2 1.41421356237
#define SQRT_pi 1.77245385090
#define SQRT_2pi 2.50662827463
//#define CONST_R 8.31
//#define Na 6.02E23
#define SYS_EPS 1E-6
#define Tmp_triple 0.694
#define n_triple 0.84

// ------------------------------------- format -----------------------------------------------------

#define EmptyStr "     " // 5 spaces
#define TabStr "    " // 4 spaces
#define PrintPres 10

//#define CompLoadMod "CompLoadMod"
//#define PPLoadMod "PostprocLoadMod"

/*
R = 1 = 3.08567758E+19 m = 2.063E+8 au = 1kpc
G = 1 = 6.67408E-11 m^3/kg/t^2
V = 1 = 1E+5 m/s
t = 1 = 3.0856776E+14 sec = 9.778E+6 y = 10My
m = 1 = 4.6234E+39 kg = 2.32436E+9 SunM
F = 1 = 1.498334E+30 H
Ro = 1 = 1.573656E-19 kg/m^3
*/

// ----------------------------------- error handeling constants ---------------------------------------
//                         ------------------ CUDA errors ----------------------
#define CUDA_ERRORS_ID 10000
#define CUDA_MEMORY_ID (CUDA_ERRORS_ID + 1000)
#define CUDA_PARAMETERS_ID (CUDA_ERRORS_ID + 2000)

#define CUDA_WRONG_NBlockW (CUDA_PARAMETERS_ID + 1) // N must be N%BlockW==0

#define CUDA_ERROR_CHECK_CONST (CUDA_ERRORS_ID + 1)

#define CUDA_WRONG_DevMem_ALLOC (CUDA_MEMORY_ID + 1)
#define CUDA_WRONG_COPY_TO_DevMem (CUDA_MEMORY_ID + 2)
#define CUDA_WRONG_COPY_TO_HostMem (CUDA_MEMORY_ID + 3)

#define CopyValues 2
#define CopyStars 3

#define BlockW 64 // for double set 256 in case of share memory amount
/*
 * max = 512
 *    (for my GTX 960 formal max is 1024,
 *     but with 1024 there is too few __shared__ avalible
 *     so strange bugs may appear with 1024);
 * formaly best is 512, but it's only ~2% faster than 256 etc.
 * ---------------------------------------------------------------------
 * must be 2^n;
 */
#define MinThreadsPerBlock 3
#define NofDspacePr 4
const int MaxN = int(pow(2,31) - 0.5);

#define WrongCompMode CompModeID
#define CVRcompMode (CompModeID+1)
#define LPcompMode (CompModeID+2)
#define RKcompMode (CompModeID+3)
#define AdCVRcompMode (CompModeID+4)
#define MYcompMode (CompModeID+5)
#define VVRcompMode (CompModeID+6)

//                         ------------------ regular errors ----------------------

#define VelocitiesID 2000
#define DtID 3000
#define CompModeID 4000
#define MemoryID 5000
#define LoadModeID 6000
#define NID 100
#define FileID 200
#define KID 300
#define PostProcModID 400
#define QID 500
#define AID 600
#define LessOrEq0ID 700
#define IsInfID 800
#define IsNanID 900
#define JustErrID -100

#define YetUnsupportedInput (JustErrID-19)
#define nCrTooBig (JustErrID-18)
#define CantCenterCM (JustErrID-17)
#define AttemptToPrintOutsyncedData (JustErrID-16)
#define NoValidTforP (JustErrID-15)
#define TooDenseSystem (JustErrID-14)
#define WrongCondProcInp (JustErrID-13)
#define WrongChngParamInp (JustErrID-12)
#define GisNot1 (JustErrID-11)
#define ConstForCommonMagic (JustErrID-10)
#define TlessCalcedT (JustErrID-9)
#define WrongScriptFormat (JustErrID-8)
#define AbiggerR (JustErrID-7)
#define StopFileAlreadyExists (JustErrID-6)
#define SayIt (JustErrID-5)

#define CantOpenFile (FileID+1)
#define CantCreateFile (FileID+2)
#define CantOpenPrevFile (FileID+3)

#define Nis1 (NID+1)
#define NLessOrEq0 (NID+2)
#define NisTooBig (NID+3)
#define NisntCubeForCristal (NID+4)

#define Qis0_1 (QID+1)
#define QLess0 (QID+2)

#define MLessOrEq0 (LessOrEq0ID+1)
#define HLess0 (LessOrEq0ID+2)
#define TLessOrEq0 (LessOrEq0ID+3)
#define dtLessOrEq0 (LessOrEq0ID+4)
#define dumpDTLessOrEq0 (LessOrEq0ID+5)
#define Geq0 (LessOrEq0ID+6)
#define krLessOrEq0 (LessOrEq0ID+7)
#define khLessOrEq0 (LessOrEq0ID+8)
#define kaLess0 (LessOrEq0ID+9)
#define kvLess0 (LessOrEq0ID+10)
#define ALessOrEq0 (LessOrEq0ID+11)
#define dissipKless0 (LessOrEq0ID+12)
#define RcutLess0 (LessOrEq0ID+13)
#define TmpLess0 (LessOrEq0ID+14)
#define muLessEq0 (LessOrEq0ID+15)
#define TmpStabEps_LessEq0 (LessOrEq0ID+16)
#define TmpStabGap_Less0 (LessOrEq0ID+17)
#define nCr_Less0 (LessOrEq0ID+18)
#define NthreadsLessOrEq0 (LessOrEq0ID+19)
#define CalcedT0Less0 (LessOrEq0ID+20)
#define FnumLess0 (LessOrEq0ID+21)

#define WrongVelocities VelocitiesID
#define SimpleVelocities (VelocitiesID+1)
#define PlammerVelocities (VelocitiesID+2)

#define DtIs0 (DtID+1)
#define DtIsLess0 (DtID+3)

#define WrongPostProcMod PostProcModID
#define EPostProcMod (PostProcModID+1)
#define FromToBinPostProcMod (PostProcModID+2)
#define ConditionPostProcMod (PostProcModID+3)

#define QIsInf (IsInfID+1)
#define AIsInf (IsInfID+2)
#define krIsInf (IsInfID+3)
#define khIsInf (IsInfID+4)
#define kaIsInf (IsInfID+5)
#define kvIsInf (IsInfID+6)
#define MIsInf (IsInfID+7)
#define HIsInf (IsInfID+8)
#define TIsInf (IsInfID+9)
#define dumpDTIsInf (IsInfID+10)
#define dtIsInf (IsInfID+11)
#define GIsInf (IsInfID+13)
#define StarXIsInf (IsInfID+14)
#define StarMIsInf (IsInfID+15)
#define StarRIsInf (IsInfID+16)
#define StarVIsInf (IsInfID+17)
#define EkIsInf (IsInfID+18)
#define EpIsInf (IsInfID+19)
#define LIsInf (IsInfID+20)
#define nCrIsInf (IsInfID+21)

#define StarXIsNan (IsNanID+1)
#define StarMIsNan (IsNanID+2)
#define StarRIsNan (IsNanID+3)
#define StarVIsNan (IsNanID+4)
#define EkIsNan (IsNanID+5)
#define EpIsNan (IsNanID+6)
#define LIsNan (IsNanID+7)

#define NULLmalloc (MemoryID+1)

#define WrongLoadMode LoadModeID
#define CompLoadMode (LoadModeID + 1)
#define PostProcLoadMode (LoadModeID + 2)
