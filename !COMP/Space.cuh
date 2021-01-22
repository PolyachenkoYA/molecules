#pragma once

#ifndef SPACE_CUH_
#define SPACE_CUH_

#include "const.cuh"
#include "general_math.cuh"
#include "format.cuh"

__global__ void kernel_FindAllA(double3 *devX, double3 *devA, double *devPressure, double R, double r_cut, int Ntot);
__global__ void kernel_FindE(double3 *devX, double3 *devV, double3 *devE, double R, double r_cut, int Ntot);

class TSpace;

// --------------------------------------------- Main class ---------------------------------------------------

class TSpace{
public:
	// ---------- Parameters ----------
	//double Tnow,BigDt,CalcedT,Tmp,dt,StartDt,CalcedT0,eff;
	double localT, dumpDT, endT, glob_time, dt, eff, n, pressure;
	double R, dissipK, TmpStabEps, r_cut2, n_cr;
	double Tmp, CONST_R, CONST_k, CONST_Na, Tmp_curr, Hshadow;
	int TmpStabGap;
	bool thermostatOn, start_crystal_lattice;
	bool Ek_is_valid, Ep_is_valid, P_is_valid, VX_synced_for_print;

	// ------ Data for computation ------
	// nCudaB - number of active CUDA blocks
	double *devPressure, *hostPressure;
	double *hostH;
	double3 *devE, *hostE;
	double3 *devX, *hostX, *hostXreal;
	double3 *devV, *hostV;
	double3 *devA, *hostA;
	double3 *hostOldX;
	double3 **hostVK, **hostRK;
	double3 E;
	double3 L, Xcm;
	double *gCondFnc;
	int Nslice;

	// ---------- Format stuff ----------
	bool binOutF, useCuda;
	int Fnum, Ntot, nCudaB, Nfrm;
	double gt, rt, ct;  // timer display
	time_t global_start_time;
	int compMode, Nthreads;
	string sessionID, compModeStr;
	string ConvFName, GlxName, PrmFName, NewGlxName, SGname, HeadFName, StpFName, LogFName,
		   LFile, EFile, TFile, ParticleFileExt, ConditionFNameBase, FramesFilesFolder, PFile;
	vector<string> ParamFHead;
	vector<string> GalaxyFHead;

	// -------- Constructor -----------
	TSpace();
	~TSpace();

	// ----- pre&post-main stuff -----
	int ResizeStars(int _n);
	int CreateParticles(void);
	int AreClosePartcls(double3 c0, int Ncurr);
	int LoadParams(string FileName);
	int SafeLoadParams(string *goodName, string firstAttemptName = "");
	int SaveParams(string FileName);
	int LoadParticles(string FileName, bool bin, int mode);
	int SaveParticles(string FileName, bool bin, int mode);
	int CreateHeadF(string FileName);
	int sayLog(string s, bool printOnScreen = 0);
	int checkInput(int _n, double _dumpDT, double _dt);
	int applyInput(int _n, double _dumpDT, double _dt);
	int ConvertFromAgama(string _conName, string _glxName, string _prmName);
	int CheckPartcls(void);
	int postGetProc(int mode);

	// -------- main computations --------
	int main_compute(void);
	int doTimeStep(int cmpMode);
	int CPU_findAllA(void);
	int doLstep(void);
	int doEstep(void);
	int computeEk(void);
	int computeEp(void);
	int computePressure(void);
	int syncForPrint(void);
	int restoreForComp(void);
	int shiftAll(void);
	int doCondStep(double dr);
	int stabilizeTmp(void);
	int useThermostat(void);
	int stopMacroMovement(void);
	int saveDump(ofstream &FoutT, ofstream &FoutE, ofstream &FoutP, string &FramePath);
	double real_time_k(void);
	int centerCM(void);
	int findCM(void);
	int findAllA(void);

	// --------- numeric schemes --------
	int stepMY(void);
	int stepCVR(void);
	int stepVVR(void);
	int stepLP(void);
	int stepRK(void);
	int stepAdCVR(void);

	// ------ CUDA implementation -----
	int GPU_findAllA(void);
	int cudaDoEstep(void);
};

#endif /* SPACE_CUH_ */
