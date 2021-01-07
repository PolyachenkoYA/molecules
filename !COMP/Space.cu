
/**
 * PolyachYA Corporation.  All rights reserved.
 *
 * Please refer to the PolyachYA end user license agreement (EULA) associated
 * with this source code for terms and conditions that govern your use of
 * this software. Any use, reproduction, disclosure, or distribution of
 * this software and related documentation outside the terms of the EULA
 * is strictly prohibited.
 */

#include "Space.cuh"

//-----------------------------------------------------------------
// -------------------------- Space--------------------------------
//-----------------------------------------------------------------
TSpace::TSpace()
{
#ifndef _OPENMP
	#warning "OpenMP unused"
#endif
	global_start_time = time(0);
	srand(global_start_time);
	sessionID = string(ctime(&global_start_time));
	gt = global_start_time - omp_get_wtime();

	n_cr = glob_time = pressure = n = dissipK = TmpStabEps = r_cut2 = localT = dumpDT = endT = dt = eff = R = Tmp = L.x = L.y = L.z = E.x = E.y = E.z = 0;
	TmpStabGap = Nslice = Nfrm = Nthreads = Ntot = Fnum = compMode = binOutF = useCuda = Ek_is_valid = Ep_is_valid = P_is_valid = 0;
	hostXreal = devX = hostX = nullptr;
	devV = hostV = nullptr;
	hostOldX = nullptr;
	hostVK = hostRK = nullptr;
	devA = hostA = nullptr;
	devE = hostE = nullptr;
	gCondFnc = nullptr;
	devPressure = hostPressure = nullptr;
	hostH = nullptr;

    ParticleFileExt = "xyz";
	SGname = "particles." + ParticleFileExt;
    PrmFName = "param.dat";
    HeadFName = "head.txt";
    ConvFName = "convert.dat";
    StpFName = "stop";
    LFile = "L.txt";
    EFile = "E.txt";
    TFile = "Time.txt";
    PFile = "Pressure.txt";
    ConditionFNameBase = "condition";
    FramesFilesFolder = "frames";
    compModeStr = "";

    ParamFHead = {"Ntot", "n", "n_cr", "r_cut", "endT", "dumpDT", "dt", "Tmp",
    		      "dissipK", "TmpStabEps", "TmpStabGap", "compMode", "binOutF", "Nthreads"};
    GalaxyFHead = {"x", "y", "z", "Vx", "Vy", "Vz", "N"};
}

TSpace::~TSpace(){
	free(hostPressure);
	free(hostH);
	free(hostX);
	free(hostXreal);
	free(hostV);
	free(hostOldX);
	free(hostA);
	free(hostE);
	if(gCondFnc) free(gCondFnc);

	if(compMode == RKcompMode){
		int i;
		if(hostVK) for(i = 0; i < 4; ++i)
			free(hostVK[i]);
		if(hostRK) for(i = 0; i < 4; ++i)
			free(hostRK[i]);
		free(hostVK);
		free(hostRK);
	}

	if(useCuda){
		CUDA_CHECK(cudaFree(devPressure));
		CUDA_CHECK(cudaFree(devX));
		CUDA_CHECK(cudaFree(devV));
		CUDA_CHECK(cudaFree(devA));
		CUDA_CHECK(cudaFree(devE));
	}
}

int TSpace::sayLog(string s, bool printOnScreen){
	return SayLog(LogFName, s, printOnScreen);
}

int TSpace::AreClosePartcls(double3 c0, int Ncurr)
{
	int i;
	//double3 r;
	double D = 2*R;
	c0 = -c0;
	for(i = 0; i < Ncurr; ++i){
		if(dot(shiftR(hostX[i] + c0, R, D)) < 1) // r(x,x0) < a
			return 1;
	}
	return 0;
}

int TSpace::CreateParticles(void)
{
	sayLog("   CreateParticles{\n");
	int i;

	sayLog("      particles generation...\n");
	if(start_crystal_lattice){ // n > n_cr -> cristal
		int i2, i3;

		int k = int(pow(Ntot/4, 1.0/3) + 0.5); // it's ok because we check if N is a perfect cube before
		double l = pow(4/n, 1.0/3);
		int n_shift = k/2;
		for(i = 0; i < k; ++i) for(i2 = 0; i2 < k; ++i2) for(i3 = 0; i3 < k; ++i3) { // fcc lattice
			hostX[(i3 + (i2 + i * k) * k) * 4] =
					make_double3(i - n_shift, i2 - n_shift, i3 - n_shift) * l;
			hostX[(i3 + (i2 + i * k) * k) * 4 + 1] =
					make_double3(i - n_shift + 0.5, i2 - n_shift + 0.5, i3 - n_shift) * l;
			hostX[(i3 + (i2 + i * k) * k) * 4 + 2] =
					make_double3(i - n_shift + 0.5, i2 - n_shift, i3 - n_shift + 0.5) * l;
			hostX[(i3 + (i2 + i * k) * k) * 4 + 3] =
					make_double3(i - n_shift, i2 - n_shift + 0.5, i3 - n_shift + 0.5) * l;
		}

		/*     cubic lattice
		int k = int(pow(Ntot, 1.0/3) + 0.5); // it's ok because we check if N is a perfect cube before
		double l = pow(n, -1.0/3);
		int n_shift = k/2;
		for(i = 0; i < k; ++i) for(i2 = 0; i2 < k; ++i2) for(i3 = 0; i3 < k; ++i3) {
			hostX[i3 + (i2 + i * k) * k] = make_double3(i - n_shift, i2 - n_shift, i3 - n_shift) * l * (1 + SYS_EPS);
		}*/
	} else { // n < n_cr -> random gas
		for(i = 0; i < Ntot; ++i){
			// too close particles can cause huge x''
			do{
				hostX[i] = make_double3(myRnd(-R, R), myRnd(-R, R), myRnd(-R, R)) * (1 - SYS_EPS);
			}while(AreClosePartcls(hostX[i], i));
		}
	}
	for(i = 0; i < Ntot; ++i)
		hostXreal[i] = hostX[i];

	for(i = 0; i < Ntot; ++i){
		hostV[i] = rndVec(gauss3D(Tmp));
	}

	sayLog("      ... particles generation DONE\n");

	i = stopMacroMovement();
	if(i) return i;
	//i = centerCM();
	//if(i) return i;

	sayLog(string("      Center of mass stabilized\n")+
		   string("   }\n"));

	return 0;
}

int TSpace::findCM(void)
{
	int i;

	Xcm = make_double3(0,0,0);
	for(i = 0; i < Ntot; ++i){
		Xcm += hostX[i];
	}
	Xcm /= Ntot;

	return 0;
}

int TSpace::centerCM(void)
{
	int i, _k = 0;
	double r_av2 = pow(n, -2.0 / 3);
	do{
		findCM();

		for(i = 0; i < Ntot; ++i){
			hostX[i] = shiftR(hostX[i] - Xcm, R);
		}

		++_k;
		if(_k > 100)
			return CantCenterCM;
	}while(dot(Xcm) > r_av2 * 0.0001);

	return 0;
}

int TSpace::stopMacroMovement(void)
{
	int i;

	double3 Vc = make_double3(0,0,0);
	for(i = 0; i < Ntot; ++i){
		Vc += hostV[i];
	}
	Vc /= Ntot;

	for(i = 0; i < Ntot; ++i){
		hostV[i] -= Vc;
	}

	return 0;
}

int TSpace::CheckPartcls(void)
{
	int i,j;
	//double3 r;
	for(i = 0; i < Ntot; ++i) for(j = i+1; j < Ntot; ++j){
		if(dot(hostX[i] - hostX[j]) < 1)
			return 1;
	}

	return 0;
}

int TSpace::ResizeStars(int _n)
{
	if(_n != Ntot){

		int szf3 = sizeof(double3)*_n;
		//int szf2 = sizeof(double2)*_n;
		int szf1 = sizeof(double)*_n;

		hostPressure = (double*)malloc(szf1);
		if(!hostPressure) return NULLmalloc;

		hostH = (double*)malloc(szf1);
		if(!hostH) return NULLmalloc;

		hostX = (double3*)malloc(szf3);
		if(!hostX) return NULLmalloc;

		hostXreal = (double3*)malloc(szf3);
		if(!hostXreal) return NULLmalloc;

		hostV = (double3*)malloc(szf3);
		if(!hostV) return NULLmalloc;

		hostOldX = (double3*)malloc(szf3);
		if(!hostOldX) return NULLmalloc;

		if(compMode == RKcompMode){
			hostVK = (double3**)malloc(sizeof(double3*)*4);
			if(!hostVK) return NULLmalloc;
			hostRK = (double3**)malloc(sizeof(double3*)*4);
			if(!hostRK) return NULLmalloc;

			for(int i = 0; i < 4; ++i){
				hostVK[i] = (double3*)malloc(szf3);
				if(!hostVK[i]) return NULLmalloc;
				hostRK[i] = (double3*)malloc(szf3);
				if(!hostRK[i]) return NULLmalloc;
			}
		}

		hostA = (double3*)malloc(szf3);
		if(!hostA) return NULLmalloc;

		hostE = (double3*)malloc(szf3);
		if(!hostE) return NULLmalloc;

		if(useCuda){
			CUDA_CHECK(cudaMalloc((void**)&devPressure,szf1));
			CUDA_CHECK(cudaMalloc((void**)&devX,szf3));
			CUDA_CHECK(cudaMalloc((void**)&devV,szf3));
			CUDA_CHECK(cudaMalloc((void**)&devA,szf3));
			CUDA_CHECK(cudaMalloc((void**)&devE,szf3));
		}

		Ntot = _n;
	}

	return 0;
}

int TSpace::CreateHeadF(string FileName)
{
	sayLog("   CreateHeadF("+FileName+"){\n");
    ofstream Fout(FileName.c_str());
    if(!Fout) return CantCreateFile;
    Fout << Fnum + 1 << "\n"
         << R << "\n"
         << binOutF << "\n"
         << Ntot << "\n"
         << rt << "\n";
    Fout.close();
    sayLog("   }created\n");
    return 0;
}

int TSpace::SaveParams(string FileName)
{
	sayLog("   SaveParams(" + FileName + "){\n");
    int i,spForV = 17;
    ofstream Fout(FileName.c_str());
    if(!Fout){return CantOpenFile;}
    Fout.precision(PrintPres);

    for(i = 0; i < ParamFHead.size(); ++i){ Fout << setw(spForV) << ParamFHead[i]; }
    Fout << '\n';
    sayLog("      text part was written\n");

    Fout << setw(spForV) << Ntot;

    vector<double> doubleOutData;
    doubleOutData = {n, n_cr, sqrt(r_cut2), endT, 1/dumpDT, 1/dt, Tmp, dissipK, TmpStabEps};
    for(i = 0; i < doubleOutData.size(); ++i) Fout << setw(spForV) << doubleOutData[i];

    vector<int> intOutData;
    intOutData = {TmpStabGap, (compMode - CompModeID) * (useCuda ? -1 : 1), binOutF, Nthreads};
    for(i = 0; i < intOutData.size(); ++i) Fout << setw(spForV) << intOutData[i];

    Fout << '\n';
    Fout.close();

    sayLog(string("      numeric part was written\n")+
    	   string("   }\n"));

    return 0;
}

int TSpace::syncForPrint(void)
{
    double dt_d2 = dt/2;
    int i;

    switch(compMode){
    	case LPcompMode:
    		for(i = 0; i < Ntot; ++i){
    			hostV[i] -= hostA[i] * dt_d2;
    		}
    		//sayLog("         velocities were modified for saving\n");
    		break;
    }

    VX_synced_for_print = 1;
    return 0;
}

int TSpace::restoreForComp(void)
{
    double dt_d2 = dt/2;
    int i;

    switch(compMode){
    	case LPcompMode:
    		for(i = 0; i < Ntot; ++i){
    			hostV[i] += hostA[i] * dt_d2;
    		}
    		//sayLog("         velocities were restored\n");
    		break;
    }

    VX_synced_for_print = 0;
    return 0;
}

int TSpace::SaveParticles(string FileName, bool bin, int mode)
{
	sayLog(string("      SaveStars (") + FileName + string("){\n"));
	int i, spForV = 18;
    ofstream Fout;

    if(bin){
        Fout.open(FileName.c_str(),ios::binary | ios::trunc);
        if(!Fout){ return CantCreateFile; }
        Fout.write((char*)&(Ntot),sizeof(int));
        Fout.write((char*)&(glob_time),sizeof(double));
    } else {
        Fout.open(FileName.c_str(),ios::trunc);
        if(!Fout){ return CantCreateFile; }
        Fout.precision(PrintPres);
        Fout << Ntot
        	 << "\nLattice=\"1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0\" Properties=pos:R:3:velo:R:3 Time=" << toString(glob_time) << "\n";
    }
    sayLog("         Header was written\n");

    for(i = 0; i < Ntot; ++i)
    {
    	if(dot(hostX[i]) == INFINITY) return StarXIsInf;
    	if(dot(hostV[i]) == INFINITY) return StarVIsInf;
    	if((hostX[i].x == NAN) || (hostX[i].y == NAN) || (hostX[i].z == NAN) || (hostX[i].x == -NAN) || (hostX[i].y == -NAN) || (hostX[i].z == -NAN)) return StarXIsNan;
    	if((hostV[i].x == NAN) || (hostV[i].y == NAN) || (hostV[i].z == NAN) || (hostV[i].x == -NAN) || (hostV[i].y == -NAN) || (hostV[i].z == -NAN)) return StarVIsNan;
        if(bin){
            Fout.write((char*)&(hostXreal[i].x),sizeof(double));
            Fout.write((char*)&(hostXreal[i].y),sizeof(double));
            Fout.write((char*)&(hostXreal[i].z),sizeof(double));
            Fout.write((char*)&(hostV[i].x),sizeof(double));
            Fout.write((char*)&(hostV[i].y),sizeof(double));
            Fout.write((char*)&(hostV[i].z),sizeof(double));
        } else {
        	Fout << scientific;
        	Fout << setw(spForV) << hostXreal[i].x
        		 << setw(spForV) << hostXreal[i].y
        		 << setw(spForV) << hostXreal[i].z
        		 << setw(spForV) << hostV[i].x
        		 << setw(spForV) << hostV[i].y
        		 << setw(spForV) << hostV[i].z
        		 << "\n";
        }
    }
    Fout.close();
    sayLog("         Stars were saved(" + FileName + ")\n");
    return 0;
}

int TSpace::checkInput(int _n, double _dumpDT, double _dt)
{
    if(compMode < CompModeID){ compMode += CompModeID; }

    if((compMode != CVRcompMode) &&
       (compMode != LPcompMode) &&
       (compMode != RKcompMode) &&
       (compMode != AdCVRcompMode) &&
       (compMode != VVRcompMode) &&
       (compMode != MYcompMode)){
    	sayLog("\ncompMode = " + toString(compMode) + "\n", 1);
    	return WrongCompMode;
    }
    if(((compMode == RKcompMode) || (compMode == AdCVRcompMode))){
    	sayLog("\nRK4 and AdCVR schemes are not supported\n", 1);
    	return YetUnsupportedInput;
    }
    if(_n == 1){ sayLog("\nN == 1\n", 1); }
    if(_n <= 0){ return NLessOrEq0; }
    if(_n > MaxN){ return NisTooBig; }
    if(useCuda && (_n % BlockW != 0)){ return CUDA_WRONG_NBlockW; }
    if(_dt <= 0){ return dtLessOrEq0; }
    if(_dt == INFINITY){ return dtIsInf; }
    if(_dumpDT <= 0){ return dumpDTLessOrEq0; }
    if(_dumpDT == INFINITY){ return dumpDTIsInf; }
    if(endT <= 0){ return TLessOrEq0; }
    if(endT == INFINITY){ return TIsInf; }
    //if(dissipK < 0){ return dissipKless0; }
    if(std::abs(dissipK) < SYS_EPS){ sayLog("\n|dissipK| < " + toString(SYS_EPS) + "\n", 1); }
    if(r_cut2 < 0){ return RcutLess0; }
    if(Tmp < 0){ return TmpLess0; } else
    if(Tmp < SYS_EPS){ sayLog("\ninitial temperature < " + toString(SYS_EPS) + "\n", 1); }
    if(TmpStabEps <= 0){ return TmpStabEps_LessEq0; }
    if(TmpStabEps < SYS_EPS){ sayLog("\nTmpStabEps < " + toString(SYS_EPS) + "\n", 1); }
    if(TmpStabGap < 0){ return TmpStabGap_Less0; }
    if(TmpStabGap > 1000000){ sayLog("\nTmpStabGap > 1000000\n", 1); }
	if((n > n_triple * (1 + SYS_EPS)) && !start_crystal_lattice){ sayLog("\nn > " + toString(n_triple) + " but crystal structure is't used. May be impossible to build the system\n", 1); }
    if(start_crystal_lattice){
    	sayLog("\n(n,T) = (" + toString(n) + ";" + toString(Tmp) + "); (n, T)_tp = (" + toString(n_cr) + ";" + toString(Tmp_triple) + ") so cristal structure will be used for initila state\n", 1);
    	if(!isIntPow(_n*2, 3)){
    		return NisntCubeForCristal;
    	}
    	//if(n > 1 + SYS_EPS){
    	//	sayLog("\nn >= 1\n", 1);
    	//	return TooDenseSystem;
    	//}
    }
    if(Nthreads < 0){ return NthreadsLessOrEq0; }
    if(Nthreads == 0){
    	sayLog("Nthreads = " + toString(omp_get_max_threads()) + " will be used");
    }
    if(Nthreads > omp_get_max_threads()){
    	sayLog("Nthreads = " + toString(Nthreads) + " > maxNthreads = " + toString(omp_get_max_threads()) + ". It is likely to be not an optimal choice\n", 1);
    }
    if(Nthreads < omp_get_max_threads()){
    	sayLog("Nthreads = " + toString(Nthreads) + " < maxNthreads = " + toString(omp_get_max_threads()) + ". It might be not an optimal choice\n", 1);
    }
    if(binOutF == 1){
    	sayLog("\nbin files aren't supported for now\n", 1);
    	return 12345;
    }

    sayLog("         Parameters were checked\n");
    return 0;
}

int TSpace::applyInput(int _n, double _dumpDT, double _dt)
{
	sayLog("      Applying input{\n");

	useCuda = compMode < 0;
	compMode = std::abs(compMode);
	start_crystal_lattice = (n > n_cr * (1 + SYS_EPS)) || (Tmp < 0.7);
	int err = checkInput(_n, _dumpDT, _dt);
	if(err) return err;
	err = ResizeStars(_n);
	if(err) return err;

	nCudaB = (Ntot + BlockW - 1) / BlockW;
	dumpDT = 1.0 / _dumpDT;
	dt = 1.0 / _dt;
	r_cut2 = r_cut2 * r_cut2;

	R = pow(Ntot/n, 1.0/3) / 2;
	if(R*R < r_cut2){ sayLog("\nR < Rcut\n", 1); }

#ifdef _OPENMP
	if(Nthreads > 0) // if == 0, max avail will be used
		omp_set_num_threads(Nthreads);
#endif

	sayLog("      }\n");
	return 0;
}

int TSpace::SafeLoadParams(string *goodName, string firstAttemptName)
{
	if(firstAttemptName.empty()){
		firstAttemptName = GlxName;
	}
    string b_s = "./" + firstAttemptName + "_" + PrmFName;
    bool b_b = access(b_s.c_str(),0);
	if(b_b){
		sayLog(string("   '" + b_s + "' not found\n") +
		       string("   trying to find './" + PrmFName + "'\n"));
		cout << "'" << b_s << "' not found\n"
			 << "trying to find './" << PrmFName << "' ... ";
		b_s = "./" + PrmFName;
	}
	*goodName = b_s;
	int err = LoadParams(b_s);
    if(b_b) cout << (err ? "failed\n" : "success\n");
    if(err) return err;
    return 0;
}

int TSpace::LoadParams(string FileName)
{
    int _N, err;
    double _dumpDT, _dt;
    string buf_s;
    sayLog("   LoadParams (" + FileName + ") {\n");

    ifstream Fin(FileName.c_str());
    if(!Fin){ return CantOpenFile; }
    std::getline(Fin, buf_s);
    Fin >> _N >> n >> n_cr >> r_cut2
    	>> endT >> _dumpDT >> _dt
    	// endT - total time to compute from 0 to endT
    	// dumpDT - time between saving to file
    	// dt - time step for computation
    	>> Tmp >> dissipK >> TmpStabEps
    	>> TmpStabGap >> compMode >> binOutF >> Nthreads;
    Fin.close();

    sayLog("      Parameters were loaded\n");

    err = applyInput(_N, _dumpDT, _dt);
    if(err) return err;

    sayLog("   }\n   ");
    sayLog(string(useCuda ? "GPU" : "CPU") + " in use\n", 1);
    return 0;
}

int TSpace::LoadParticles(string FileName, bool bin, int mode)
{
    int i,N,err;
    string s_buf;
    //double d_buf;
    ifstream Fin;
    sayLog("   LoadParticles (" + FileName + ") {\n");
    if(bin){
        Fin.open(FileName.c_str(),ios::binary);
        if(!Fin){ return CantOpenFile; }
        Fin.read((char*)&(N), sizeof(int));
        Fin.read((char*)&(glob_time), sizeof(double));
    } else {
        Fin.open(FileName.c_str());
        if(!Fin){ return CantOpenFile; }
        Fin >> N;
        getline(Fin, s_buf);
        getline(Fin, s_buf);
    }
    if(N <= 0) return NLessOrEq0;
    err = ResizeStars(N);
    if(err) return err;
    sayLog("      Head data was loaded\n");
    for(i = 0; i < Ntot; ++i)
    {
        if(bin){
            Fin.read((char*)&(hostXreal[i].x), sizeof(double));
            Fin.read((char*)&(hostXreal[i].y), sizeof(double));
            Fin.read((char*)&(hostXreal[i].z), sizeof(double));
            Fin.read((char*)&(hostV[i].x), sizeof(double));
            Fin.read((char*)&(hostV[i].y), sizeof(double));
            Fin.read((char*)&(hostV[i].z), sizeof(double));
        } else {
        	Fin >> hostXreal[i].x >> hostXreal[i].y >> hostXreal[i].z >> hostV[i].x >> hostV[i].y >> hostV[i].z;
        }
        hostX[i] = shiftRtrue(hostXreal[i], R);
    }
    Fin.close();
    sayLog("      Stars were loaded\n");

    postGetProc(mode);

    sayLog("   }\n");
    return 0;
}

int TSpace::postGetProc(int mode)
{
	int i;

	if(useCuda){
		CUDA_CHECK(cudaMemcpy(devX,hostX,Ntot*sizeof(double3),cudaMemcpyHostToDevice));
		sayLog("      devX and devM were filled with host data\n");
	}

	if(mode == CompLoadMode){
		i = useCuda ? GPU_findAllA() : CPU_findAllA();
		if(i) return i;
		sayLog("      Acceleration were computed\n");

		double dt2_d2;
		switch(compMode){
		case CVRcompMode:
			dt2_d2 = dt*dt*0.5;
			for(i = 0; i < Ntot; ++i) hostOldX[i] = shiftR(hostX[i] - hostV[i]*dt + hostA[i]*dt2_d2, R);
			sayLog("      oldCrd for CVRmode were computed\n");
			break;
		case LPcompMode:
			dt2_d2 = dt*0.5;
			for(i = 0; i < Ntot; ++i) hostV[i] += hostA[i]*dt2_d2;
			sayLog("      v(1/2dt) for LPmode were computed\n");
			break;
		}
	}

	return 0;
}

int TSpace::stabilizeTmp(void)
{
	sayLog("   stabTmp{\n");

	int b_i, i;
	long int print_step, k_step = 0;
	double Tmp_av, std_disp, Tmp_av0;
	double t = 0;
    double start_time_glob, start_time_2;
    double *Tmp_inst = new double[TmpStabGap];

    b_i = postGetProc(CompLoadMode);
    if(b_i) return b_i;
    Tmp_av = Tmp;
    print_step = Ntot > 10000 ? 1 : 100000000 / (Ntot * Ntot);

    //time(&start_time_glob);
    start_time_glob = omp_get_wtime();
    do
    {
    	b_i = doTimeStep(compMode);
    	if(b_i) return b_i;

    	syncForPrint();
    	computeEk();
    	restoreForComp();

    	Tmp_av = (Tmp_av * (TmpStabGap - 1) + Tmp_curr) / TmpStabGap;
    	Tmp_inst[k_step % TmpStabGap] = Tmp_curr;
		if(k_step == TmpStabGap){
			Tmp_av0 = Tmp_av;
			//time(&start_time_2);
			start_time_2 = omp_get_wtime();
		}
    	if(k_step >= TmpStabGap){
    		std_disp = 0;
    		for(i = 0; i < TmpStabGap; ++i){
    			std_disp += pow2(Tmp_av - Tmp_inst[i]);
    		}
    		std_disp = sqrt(std_disp / TmpStabGap) / Tmp_av;

    	}

    	++k_step;
    	t += dt;

        if(k_step % print_step == 0)
        {
        	if(k_step < TmpStabGap){
            	time_progress(gt, start_time_glob, omp_get_wtime(), ((double)k_step / TmpStabGap),
            			"Tmp stabilizing: initial time gap\nT_target = " + toString(Tmp) + "; T_curr_av = " + toString(Tmp_av) + "; T_curr = " + toString(Tmp_curr),
            			1);

        	} else {
            	time_progress(gt, start_time_2, omp_get_wtime(),
            			log(epsDlt(Tmp_av0, Tmp) / epsDlt(Tmp_av, Tmp)) / log(epsDlt(Tmp_av0, Tmp) / TmpStabEps),
            			"Tmp stabilizing: waiting for equilibrium\nT_target = " + toString(Tmp) + "; T_curr_av = " + toString(Tmp_av) + "; T_curr = " + toString(Tmp_curr),
            			1);
        	}
        }
    }while(!(almostEq(Tmp_av, Tmp, TmpStabEps) && almostEq(Tmp_curr, Tmp, TmpStabEps) && (std_disp < TmpStabEps) && (k_step > TmpStabGap)));

    //for(i = 0; i < this->Ntot; ++i){
    //	this->hostV[i] *= sqrt(Tmp / Tmp_curr);
    //}

    delete [] Tmp_inst;
    sayLog("Tmp_target = " + toString(Tmp) + "; Tmp_av = " + toString(Tmp_av) + "; Tmp = " + toString(Tmp_curr) + "                        \n", 1);
    sayLog("   }\n");
    Tmp = Tmp_curr;
    return 0;
}

int TSpace::saveDump(ofstream &FoutT, ofstream &FoutE, ofstream &FoutP, string &FramePath)
{
	int b_i;
	syncForPrint();

	sayLog("      Saving " + toString(++Fnum) + " file...\n");
    FoutT << glob_time << endl;

    b_i = SaveParticles(FramePath + toString(Fnum) + string(".") + ParticleFileExt, binOutF, CompLoadMode);
    if(b_i) return b_i;
    sayLog("      ... saved\n");

    b_i = useCuda ? cudaDoEstep() : doEstep();
    if(b_i) return b_i;

    FoutE << E.x << " " << E.y << " " << (E.x + E.y) << " " << (E.x + E.y + E.z) << endl;
    sayLog("      Energy computed & saved\n");

	b_i = computePressure();
	if(b_i) return b_i;
    FoutP << pressure << endl;
    sayLog("      Pressure computed & saved\n");

    restoreForComp();

    return 0;
}

int TSpace::main_compute(void)
{
	sayLog("   main_compute{\n");
	if(!access(StpFName.c_str(),0)) return StopFileAlreadyExists;

    int b_i = 0;
    double rStime;
    double kt;
//    bool stpCalc = 0;
    string buf_s, BaseName = "./" + GlxName + "/";
    string FramePath = BaseName + FramesFilesFolder + "/";

    //ofstream FoutT((BaseName+TFile).c_str(),ios::app);
    // This was to make it possible to continue computation, but here it's not possible for other reasons.
    // so no need to keep this bothering feature
    ofstream FoutT((BaseName + TFile).c_str());
    if(!FoutT) return CantCreateFile;
    FoutT.precision(PrintPres);

    ofstream FoutP((BaseName + PFile).c_str());
    if(!FoutP) return CantCreateFile;
    FoutP.precision(PrintPres);

    ofstream FoutE((BaseName + EFile).c_str());
    if(!FoutE) return CantCreateFile;
    FoutE.precision(PrintPres);

    --Fnum;
	b_i = saveDump(FoutT, FoutE, FoutP, FramePath);
	if(b_i) return b_i;

    cout << "Model name : '" << GlxName << "'\n";

    kt = real_time_k();

    sayLog(string("      " + compModeStr + " computation mode\n") +
    	   string("      thermostat is " + toString(thermostatOn ? "used" : "not used") + "\n") +
    	   string("      Starting the computation{\n"));
    rStime = omp_get_wtime();
    do
    {
    	b_i = doTimeStep(compMode);
    	if(b_i) return b_i;
    	localT += dt;
    	glob_time += dt;

        if(localT > dumpDT*(1 - SYS_EPS))
        {
        	localT = 0;

        	b_i = saveDump(FoutT, FoutE, FoutP, FramePath);
        	if(b_i) return b_i;

            time_progress(gt, rStime, omp_get_wtime(), glob_time / endT / kt, "computing");
        }
    }while((glob_time < endT - dt*(1 - SYS_EPS)) && access(StpFName.c_str(),0));
    rt = omp_get_wtime() - rStime;
    FoutT.close();
    FoutP.close();
    FoutE.close();

    eff = glob_time/dt * Ntot*Ntot / rt*3600;
    //if(!useCuda) eff/=2;
    buf_s = glob_time < endT - dt * (1 - SYS_EPS) ? "terminated" : "done";
    buf_s = "Computation was " + buf_s;
    cout << buf_s << "\n";
    sayLog(string("      " + buf_s + "\n") +
    	   string("      }\n") +
    	   string("   }\n") +
    	   string("   Efficiency{\n") +
    	   string("      time computed = " + toString(glob_time) + "; dt = " + toString(dt) + "; usedT/dt = " + toString(glob_time / dt)+"\n") +
           string("      N = " + toString(Ntot) + "; real_t [h] = " + toString(rt / 3600.0) + "\n") +
           string("      efficiency (e = endT/dt*N^2/real_t) = " + toString(eff) + "\n") +
           string("   }\n"));


    return 0;
}

double TSpace::real_time_k(void)
{
    /*
     * kt is for accounting computeEp effect on time
     * */
    double kt = dt/dumpDT;
    switch(compMode){
		case MYcompMode:
			compModeStr = "MY";
			break;
    	case CVRcompMode:
    		compModeStr = "CVR";
    		break;
    	case LPcompMode:
    		compModeStr = "LP";
    		break;
    	case RKcompMode:
    		kt /= 4;
    		compModeStr = "RK";
    		break;
    	case AdCVRcompMode:
    		kt /= 4;
    		compModeStr = "AdCVR";
    		break;
    }
    return 1 + kt * 0.35; // TODO find k
}

int TSpace::shiftAll(void)
{
	double D = 2*R;
    for(int i = 0; i < Ntot; ++i) hostX[i] = shiftR(hostX[i], R, D);
    return 0;
}

int TSpace::doTimeStep(int cmpMode)
{
    switch (cmpMode){
        case CVRcompMode: stepCVR(); break;
        case LPcompMode: stepLP(); break;
        case RKcompMode: stepRK(); break;
        case AdCVRcompMode: stepAdCVR(); break;
    	case MYcompMode: stepMY(); break;
    	case VVRcompMode: stepVVR(); break;
        default: return WrongCompMode;
    }

    // shiftAll();

    // this is necessary here - without it there is a slow drift of CM
    stopMacroMovement();

    Ek_is_valid = 0;
    Ep_is_valid = 0;
    P_is_valid = 0;
    VX_synced_for_print = 0;

    return 0;
}

int TSpace::CPU_findAllA(void)
{
    int i,j;
    double3 bv, m_x_curr;
    double r2, _f;
    double D = 2 * R;

   for(i = 0; i < Ntot; ++i){
	   hostA[i].x = hostA[i].y = hostA[i].z = 0;
	   hostPressure[i] = 0;
   }

    #ifdef _OPENMP
    	#pragma omp parallel for private(bv, r2, j, m_x_curr, _f)
    #endif
    for(i = 0; i < Ntot; ++i){
    	m_x_curr = -hostX[i];
        for(j = 0; j < i; ++j){
        	 bv = shiftR(hostX[j] + m_x_curr, R, D);
        	 r2 = dot(bv, bv);
        	 if(r2 < r_cut2){
        		 hostA[i] += bv * getForce(r2);
        	 }
        }
        for(j = i + 1; j < Ntot; ++j){
        	 bv = shiftR(hostX[j] + m_x_curr, R, D);
        	 r2 = dot(bv, bv);
        	 if(r2 < r_cut2){
        		 _f = getForce(r2);
        		 hostPressure[i] += r2 * _f;
        		 hostA[i] += bv * _f;
        	 }
        }
        hostPressure[i] = -hostPressure[i];
    }

    return 0;
}

int TSpace::useThermostat(void)
{
    if(thermostatOn){
    	double b_d, sT = sqrt(Tmp);
    	//double p_coll = dt * dissipK;
    	int i;

    	b_d = sqrt(2 * Tmp * dissipK / dt);
    	// there shouldn't be *2 in the end. For some reason Temperature converges to Tmp/2 without it.
    	for(i = 0; i < Ntot; ++i){ // this should be auto-paralleled
    		hostA[i] += (make_double3(gaussRand(b_d), gaussRand(b_d), gaussRand(b_d)) - hostV[i] * dissipK); // Langevin
    		//if(myRnd() < p_coll){
    		//	hostV[i] = make_double3(gaussRand(sT), gaussRand(sT), gaussRand(sT)); // Andersen
    		//}
    	}
    }

    return 0;
}

int TSpace::stepVVR(void)
{
    int i;
    double3 dx;
    double D = 2 * R;

#ifdef _OPENMP
	#pragma omp parallel for private(dx)
#endif
    for(i = 0; i < Ntot; ++i){
    	hostV[i] += hostA[i] * (dt * 0.5);
    	dx = hostV[i] * dt;
    	hostX[i] = shiftR(hostX[i] + dx, R, D);
    	hostXreal[i] += dx;
    }

    findAllA();

    for(i = 0; i < Ntot; ++i)
    	hostV[i] += hostA[i] * (dt * 0.5);

    return 0;
}

int TSpace::stepMY(void)
// This one may seem similar to LP one. But it's not. In LP we have v(n+1/2) and x(n) and here we have v(n) and x(n)
{
    int i;
    double D = 2 * R;
    double3 dx;

#ifdef _OPENMP
	#pragma omp parallel for private(dx)
#endif
    for(i = 0; i < Ntot; ++i){
    	dx = (hostV[i] + hostA[i] * (dt * 0.5)) * dt;
    	hostX[i] = shiftR(hostX[i] + dx, R, D);
    	hostXreal[i] += dx;
    	hostOldX[i] = hostA[i]; // here it is actually old forces but not old coordinates
    }

    findAllA();

    for(i = 0; i < Ntot; ++i)
    	hostV[i] += (hostA[i] + hostOldX[i]) * (0.5 * dt);

    return 0;
}

int TSpace::stepCVR(void)
{
    int i;
    double3 bv, dx;
    double D = 2*R, dt2 = dt*dt;

    findAllA();
#ifdef _OPENMP
	#pragma omp parallel for private(bv, dx)
#endif
    for(i = 0; i < Ntot; ++i){
        bv = hostX[i];
        dx = hostX[i] - hostOldX[i] + hostA[i]*dt2;
        hostX[i] = shiftR(hostX[i] + dx, R, D);
        // here we use 3 points - X[i-1], X[i], X[i+1]. Periodic boundaries are not checked to work correctly
        hostXreal[i] += dx;
        hostV[i] = shiftR(hostX[i] - hostOldX[i], R, D)/(2*dt);
        hostOldX[i] = bv;
    }
    return 0;
}

int TSpace::stepLP(void)
{
    int i;
    double D = 2*R;
    double3 dx;

#ifdef _OPENMP
	#pragma omp parallel for private(dx)
#endif
    // V is 1/2 time-step ahead of X
    for(i = 0; i < Ntot; ++i){
    	dx = hostV[i] * dt;
    	hostX[i] = shiftR(hostX[i] + dx, R, D); // finish previous step by modifying X
    	hostXreal[i] += dx;
    }

    findAllA();                                        // compute new A

    for(i = 0; i < Ntot; ++i)
    	hostV[i] += hostA[i]*dt; // compute new V

    return 0;
}

int TSpace::stepRK(void)
{
    int i;
    double dt_d2 = dt*0.5;

    findAllA();
    for(i=0;i<Ntot;++i){ // Here all loops should auto-paralelize themselves
    	hostOldX[i] = hostX[i];

    	hostVK[0][i] = hostA[i];
    	hostRK[0][i] = hostV[i];
    	hostX[i] = hostOldX[i] + hostRK[0][i] * dt_d2;
    }

    findAllA();
    for(i=0;i<Ntot;++i){
    	hostVK[1][i] = hostA[i];
    	hostRK[1][i] = hostV[i] + hostVK[0][i] * dt_d2;
    	hostX[i] = hostOldX[i] + hostRK[1][i] * dt_d2;
    }

    findAllA();
    for(i=0;i<Ntot;++i){
    	hostVK[2][i] = hostA[i];
    	hostRK[2][i] = hostV[i] + hostVK[1][i] * dt_d2;
    	hostX[i] = hostOldX[i] + hostRK[2][i] * dt_d2;
    }

    findAllA();
	for(i = 0; i < Ntot; ++i){
    	hostVK[3][i] = hostA[i];
    	hostRK[3][i] = hostV[i] + hostVK[2][i] * dt;
		hostA[i] = (hostVK[0][i] + 2 * hostVK[1][i] + 2 * hostVK[2][i] + hostVK[3][i]) / 6;
		hostX[i] = hostOldX[i] + (hostRK[0][i] + 2 * hostRK[1][i] + 2 * hostRK[2][i] + hostRK[3][i]) * (dt/6);
	}

    for(i = 0; i < Ntot; ++i){
    	hostV[i] += hostA[i] * dt;
    }

    return 0;
}

int TSpace::stepAdCVR(void)
{
    int i;
    double k1 = 0.1786178958448091;
    double k2 = -0.2123418310626054;
    double k3 = -0.0662645826698185;

    for(i=0;i<Ntot;++i) hostX[i] += hostV[i]*(k1*dt);
    findAllA();
    for(i=0;i<Ntot;++i){
    	hostV[i] += hostA[i]*((0.5-k2)*dt);
    	hostX[i] += hostV[i]*(dt*k3);
    }
    findAllA();
    for(i=0;i<Ntot;++i){
    	hostV[i] += hostA[i]*(k2*dt);
    	hostX[i] += hostV[i]*((1-2*(k1+k3))*dt);
    }
    findAllA();
    for(i=0;i<Ntot;++i){
    	hostV[i] += hostA[i]*(k2*dt);
    	hostX[i] += hostV[i]*(k3*dt);
    }
    findAllA();
    for(i=0;i<Ntot;++i){
    	hostV[i] += hostA[i]*((0.5-k2)*dt);
    	hostX[i] += hostV[i]*(k1*dt);
    }

    return 0;
}

/*
*/
// TSpace
//------------------------------------------------------------------------------

int TSpace::doCondStep(double dr)
{
	int i, j, k;
	double D = 2*R;

	for(i = 0; i < Nslice; ++i){
		gCondFnc[i] = 0;
	}

	for(i = 0; i < Ntot; ++i){
		for(j = 0; j < i; ++j){
			k = floor( length( shiftR(hostX[j] - hostX[i], R, D) )/dr );
			if(k < Nslice)
				++gCondFnc[k];
		}
		// every particle except i0
		for(i = i+1; i < Ntot; ++i){
			k = floor( length( shiftR(hostX[j] - hostX[i], R, D) )/dr );
			if(k < Nslice)
				++gCondFnc[k];
		}
	}

	double kg = 4 * pi * dr*dr * dr * (Ntot/pow3(2*R));
	for(i = 0; i < Nslice; ++i){
		gCondFnc[i] /= (kg * (i+0.5)*(i+0.5));
	}

	return 0;
}

int TSpace::computeEk(void)
{
	E.x = 0;
	for(int i = 0; i < Ntot; ++i){
		E.x += dot(hostV[i]);
	}
    E.x *= 0.5;
	if(std::abs(E.x) == INFINITY) return EkIsInf;
	if((E.x == NAN) || (E.x == -NAN)) return EkIsNan;

	Tmp_curr = 2.0/3 * E.x/Ntot;

	Ek_is_valid = 1;
	return 0;
}

int TSpace::computePressure(void)
// Forces and Ek are expected to be relevant
{
	if(!Ek_is_valid)
		return NoValidTforP;

	int i;

	pressure = 0;
	for(i = 0; i < Ntot; ++i){
		pressure += hostPressure[i];
	}
	pressure /= (3 * pow3(2 * R));
	pressure += n * Tmp_curr;

	P_is_valid = 1;
	return 0;
}

int TSpace::computeEp(void)
{
	unsigned int i, j;
	double3 r, m_r_curr, m_v_curr, dv;
	double r2, E_loc = 0, H_loc = 0;
	double D = 2*R;

// this openmp is checked - it has no influence on accuracy
#ifdef _OPENMP
	#pragma omp parallel for private(r2, j, r, m_r_curr, m_v_curr, dv) schedule(guided) reduction(+: E_loc, H_loc)
#endif
	for(i = 0; i < Ntot; ++i){
		m_r_curr = -hostX[i];
		m_v_curr = -hostV[i];
		for(j = 0; j < i; ++j){
			r = shiftR(hostX[j] + m_r_curr, R, D);
			r2 = dot(r,r);
			if(r2 < r_cut2){
				r2 = 1 / r2;
				dv = hostV[j] + m_v_curr;
				H_loc += getd2H(r2, r, dv, m_v_curr);
				E_loc += getEp(r2);
			}
		}
		for(j = i + 1; j < Ntot; ++j){
			r = shiftR(hostX[j] + m_r_curr, R, D);
			r2 = dot(r);
			if(r2 < r_cut2){
				r2 = 1 / r2;
				dv = hostV[j] + m_v_curr;
				H_loc += getd2H(r2, r, dv, m_v_curr);
			}
		}
	}
	// E_loc and H_loc are already summed thanks to reduction(+:E_loc) in the OMP initiation
	E.y = E_loc - getEp(1 / r_cut2) * Ntot * (Ntot - 1) / 2; // substract Ecut

	E.z = 0;
	for(i = 0; i < Ntot; ++i)
		E.z += dot(hostA[i], hostA[i]);
	E.z = (H_loc * 2 - E.z / 2) * (dt*dt/12);

	if(std::abs(E.y) == INFINITY) return EpIsInf;
	if((E.y == NAN) || (E.y == -NAN)) return EpIsNan;

	Ep_is_valid = 1;
	return 0;
}

int TSpace::doEstep(void)
{
	if(!VX_synced_for_print)
		return AttemptToPrintOutsyncedData;

	int err_handl = computeEk();
	if(err_handl) return err_handl;
	err_handl = computeEp();
	if(err_handl) return err_handl;

	return 0;
}

int TSpace::cudaDoEstep(void)
{
	if(!VX_synced_for_print)
		return AttemptToPrintOutsyncedData;

	int i;
	CUDA_CHECK(cudaMemcpy(devV, hostV, Ntot * sizeof(double3), cudaMemcpyHostToDevice));
	kernel_FindE<<< nCudaB,BlockW >>>(devX, devV, devE, R, r_cut2, Ntot);
	CUDA_CHECK(cudaDeviceSynchronize());
	CUDA_CHECK(cudaMemcpy(hostE, devE, nCudaB * sizeof(double3), cudaMemcpyDeviceToHost));
	E.x = E.y = E.z = 0;
	for(i = 0; i < nCudaB; ++i) E += hostE[i];

	E.y /= 2; // because we summed all pairs 2 twice in CUDA

	double f2 = 0;
	for(i = 0; i < Ntot; ++i) f2 += dot(hostA[i], hostA[i]);
	E.z = (E.z * 2 - f2 / 2) * (dt * dt / 12);

	Tmp_curr = 2.0/3 * E.x / Ntot;
	Ek_is_valid = 1;
	Ep_is_valid = 1;

	if(E.x == INFINITY) return EkIsInf;
	if(E.x == NAN) return EkIsNan;
	if(std::abs(E.y) == INFINITY) return EpIsInf;
	if((E.y == NAN) || (E.y == -NAN)) return EpIsNan;

    return 0;
}

__global__ void kernel_FindE(double3 *devX, double3 *devV, double3 *devE, double R, double r_cut2, int Ntot)
{
	__shared__ double3 c1[BlockW], c2[BlockW], v1[BlockW], v2[BlockW];
	__shared__ double3 e1[BlockW];
	int gind = BlockW * blockIdx.x + threadIdx.x;

	if(gind >= Ntot) return;

	int tind, tile, i, tx = threadIdx.x;
	double3 r, dv;
	double r2, D = 2 * R;
	c1[tx] = -devX[gind];
	v1[tx] = -devV[gind];
	e1[tx].x = dot(v1[tx], v1[tx]) * 0.5;
	e1[tx].y = 0;
	e1[tx].z = 0;

	for(tile = 0; tile * BlockW + tx < Ntot; ++tile){
		tind = tile * BlockW + tx;
		c2[tx] = devX[tind];
		v2[tx] = devV[tind];
		__syncthreads();
		for(i = 0; i < BlockW; ++i){
			if(tile * BlockW + i != gind){
				r = shiftR(c2[i] + c1[tx], R, D);
				if(r2 < r_cut2){
					dv = v2[i] + v1[tx];
					r2 = 1 / dot(r, r);
					e1[tx].y += getEp(r2);
					e1[tx].z += getd2H(r2, r, dv, v1[tx]);
				}
			}
		}
		__syncthreads();
	}

	i = BlockW / 2;
	while(i > 0){
		if(tx < i){
			tind = tx + i;
			e1[tx].x += e1[tind].x;
			e1[tx].y += e1[tind].y;
			e1[tx].z += e1[tind].z;
		}
		__syncthreads();
		i /= 2;
		/*
		*	For correct work Ntot%BlockW must be ==0, because if it's not, then
		*	/=2 won't properly wrap indexes
		*/
	}

	if(tx==0){
		devE[blockIdx.x].x = e1[0].x;
		devE[blockIdx.x].y = e1[0].y;
		devE[blockIdx.x].z = e1[0].z;
	}
}
/*
 */

int TSpace::GPU_findAllA(void)
{
	CUDA_CHECK(cudaMemcpy(devX, hostX, Ntot*sizeof(double3), cudaMemcpyHostToDevice));
	kernel_FindAllA<<< nCudaB,BlockW >>>(devX, devA, devPressure, R, r_cut2, Ntot);
	CUDA_CHECK(cudaDeviceSynchronize());
	CUDA_CHECK(cudaMemcpy(hostA, devA, Ntot*sizeof(double3), cudaMemcpyDeviceToHost));
	CUDA_CHECK(cudaMemcpy(hostPressure, devPressure, Ntot*sizeof(double1), cudaMemcpyDeviceToHost));

	return 0;
}

int TSpace::findAllA(void)
{
	int res = useCuda ? GPU_findAllA() : CPU_findAllA();
	useThermostat();
	return res;
}

__global__ void kernel_FindAllA(double3 *devX, double3 *devA, double *devPressure, double R, double r_cut2, int Ntot)
{
	__shared__ double3 c1[BlockW], c2[BlockW], a1[BlockW];
	__shared__ double p[BlockW];
	int gind = BlockW * blockIdx.x + threadIdx.x; // global index

	if(gind >= Ntot) return;
	/*
	*	if(Ntot%BlockW != 0) then some threads in one block will never reach
	*	__syncthreads() command because of this "if" statement, and it's bad. So Ntot%BlockW must be == 0
	*/

	// set start values
	int tind, c2_tind, tile, i, tx = threadIdx.x;
	double3 r;
	double _f, r2, D = 2 * R;

	// copy "left column" of stars from global - the ones I'll sum TO
	c1[tx] = -devX[gind];
	a1[tx].x = a1[tx].y = a1[tx].z = 0;
	p[tx] = 0;

	for(tile = 0; tile*BlockW + tx < Ntot; ++tile){
		tind = tile*BlockW + tx;
		// copy part of "up row" from global - the 2nd part of current tile - the ones I'll sum WITH
		c2[tx] = devX[tind];
		// m2[tx] = devM[tind];
		__syncthreads();

		for(i = 0; i < BlockW; ++i){
			c2_tind = tile*BlockW + i;
			if(c2_tind != gind){
				r = shiftR(c2[i] + c1[tx], R, D);
				r2 = dot(r, r);
				if(r2 < r_cut2){
					_f = getForce(r2);
					if(c2_tind > gind)
						p[tx] += r2 * _f;
					a1[tx] += r * _f;
				}
			}
		}

		__syncthreads();
	}

	devA[gind] = a1[tx];
	devPressure[gind] = -p[tx];
}
