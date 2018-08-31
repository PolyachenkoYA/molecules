/*
 * Space.cuh
 *
 *  Created on: Jul 7, 2016
 *      Author: polyachyadt
 */

#pragma once

#ifndef SPACE_CUH_
#define SPACE_CUH_

#include "const.cuh"
#include "general_math.cuh"
#include "format.cuh"

__global__ void cudaFindAllA(double3 *devX, double3 *devA, double m, double R, double r_brd, int Ntot);
__global__ void cudaFindE(double3 *devX, double3 *devV, double2 *devE, double m, double R, double r_brd, int Ntot);

class TSpace;

// --------------------------------------------- Main class ---------------------------------------------------

class TSpace{
public:
	// ---------- Parameters ----------
	//double Tnow,BigDt,CalcedT,Tmp,dt,StartDt,CalcedT0,eff;
	double localT, dumpDT, totalT, totalT0, endT, dt, eff;
	double a, R, k12, k6, kv, dissipK, TmpStabEps, r_brd;
	double Tmp, mu, CONST_R, CONST_k, CONST_Na, Tmp_curr;
	int TmpStabGap;
	bool thermostatOn;

	// ------ Data for computation ------
	// nCudaB - number of active CUDA blocks
	double3 *devX, *hostX;
	double3 *devV, *hostV;
	double3 *devA, *hostA;
	double3 *hostF;
	double3 *hostOldX;
	double3 **hostVK, **hostRK;
	double2 *devE, *hostE;
	double *devM, *hostM;
	double2 E;
	double3 L;
	double *gCondFnc;
	int Nslice;

	// ---------- Format stuff ----------
	bool binOutF, useCuda, thermostatIsOn;
	int Fnum, Ntot, nCudaB, Nfrm;
	time_t rt, ct; // timer display
	int compMode;
	string sessionID;
	string ConvFName, GlxName, PrmFName, SGname, HeadFName, StpFName, LogFName, LFile, EFile, TFile, ParticleFileExt, ConditionFNameBase, FramesFilesFolder;
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
	int SafeLoadParams(string *goodName);
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
	double3 findA(double3 crd1, unsigned long int N);
	int findAllA(void);
	int doLstep(void);
	int doEstep(void);
	int computeEk(void);
	int computeEp(void);
	int shiftAll(void);
	int doCondStep(double dr);
	int stabilizeTmp(void);
	int useThermostat(void);
	int stopMacroMovement(void);

	// --------- numeric schemes --------
	int MYmove(void);
	int VRmove(void);
	int LPmove(void);
	int RKmove(void);
	int AdVRmove(void);
	int MyEGmove(void);

	// ------ CUDA numeric schemes -----
	int cudaMYmove(void);
	int cudaVRmove(void);
	int cudaLPmove(void);
	int cudaRKmove(void);
	int cudaAdVRmove(void);
	int cudaDoTimeStep(int cmpMode);
	int cudaDoEstep(void);
};

#endif /* SPACE_CUH_ */
