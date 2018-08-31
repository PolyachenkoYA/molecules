
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
	srand(time(0));

	dissipK = TmpStabEps = r_brd = localT = dumpDT = totalT = endT = dt = eff = a = R = k12 = k6 = Tmp = mu = L.x = L.y = L.z = E.x = E.y = 0;
	TmpStabGap = Nslice = Nfrm = Ntot = Fnum = compMode = binOutF = useCuda = 0;
	devX = hostX = nullptr;
	devV = hostV = nullptr;
	hostF = nullptr;
	hostOldX = nullptr;
	hostVK = hostRK = nullptr;
	devA = hostA = nullptr;
	devE = hostE = nullptr;
	devM = hostM = nullptr;
	gCondFnc = nullptr;

	time_t now = time(0);
	sessionID = string(ctime(&now));

    ParticleFileExt = "xyz";
	SGname = "particles." + ParticleFileExt;
    PrmFName = "param.dat";
    HeadFName = "head.txt";
    ConvFName = "convert.dat";
    StpFName = "stop";
    LFile = "L.txt";
    EFile = "E.txt";
    TFile = "Time.txt";
    ConditionFNameBase = "condition";
    FramesFilesFolder = "frames";

    ParamFHead = {"Ntot","k12","k6","R","r_max","endT","totalT","dumpDT","dt","Tmp","mu","kv","dissipK","TmpStabEps","TmpStabGap","compMode","binOutF"};
    GalaxyFHead = {"x","y","z","Vx","Vy","Vz","m","N","Fnum","totalT"};
}

TSpace::~TSpace(){
	free(hostX);
	free(hostV);
	free(hostF);
	free(hostOldX);
	free(hostA);
	free(hostE);
	free(hostM);
	if(gCondFnc) free(gCondFnc);

	int i;
	if(hostVK) for(i = 0; i < 4; ++i)
		free(hostVK[i]);
	if(hostRK) for(i = 0; i < 4; ++i)
		free(hostRK[i]);
	free(hostVK);
	free(hostRK);

	if(useCuda){
		CUDA_CHECK(cudaFree(devX));
		CUDA_CHECK(cudaFree(devV));
		CUDA_CHECK(cudaFree(devA));
		CUDA_CHECK(cudaFree(devE));
		CUDA_CHECK(cudaFree(devM));
	}
}

int TSpace::sayLog(string s, bool printOnScreen){
	return SayLog(LogFName, s, printOnScreen);
}

int TSpace::AreClosePartcls(double3 c0, int Ncurr)
{
	int i;
	double3 r;
	double D = 2*R;
	for(i = 0; i < Ncurr; ++i){
		r = shiftR(hostX[i] - c0, R, D);
		if(dot(r,r) < 1) // r(x,x0) < a
			return 1;
	}
	return 0;
}

int TSpace::CreateParticles(void)
{
	sayLog("   CreateParticles{\n");
	int i;
	double v0, Vmn, Vmx;

	v0 = sqrt(3*Tmp/mu/(1 + kv*kv/3));
	Vmn = (1 - kv)*v0;
	Vmx = (1 + kv)*v0;
	// 3kT = m * <V^2>

	sayLog("      particles generation...\n");
	for(i = 0; i < Ntot; ++i){
		// too close particles can cause huge x''
		do{
			hostX[i] = make_double3(myRnd(-R, R), myRnd(-R, R), myRnd(-R, R));
		}while(AreClosePartcls(hostX[i], i));

		hostV[i] = rndVec(myRnd(Vmn, Vmx));

		hostM[i] = mu;
	}

	sayLog("      ... particles generation DONE\n");

	stopMacroMovement();

	sayLog(string("      Center of mass was stabilized\n")+
		   string("   }\n"));

	return 0;
}

int TSpace::stopMacroMovement(void)
{
	int i;

	double3 Vc = make_double3(0,0,0);
	double Mc = 0;
	for(i = 0; i < Ntot; ++i){
		Mc += hostM[i];
		Vc += hostV[i]*hostM[i];
	}
	Vc /= Mc;

	for(i = 0; i < Ntot; ++i){
		hostV[i] -= Vc;
	}

	return 0;
}

int TSpace::CheckPartcls(void)
{
	int i,j;
	double3 r;
	for(i = 0; i < Ntot; ++i) for(j = i+1; j < Ntot; ++j){
		r = hostX[i] - hostX[j];
		if(dot(r) < 1){
			return 1;
		}
	}

	return 0;
}

int TSpace::ResizeStars(int _n)
{
	if(_n!=Ntot){
		int szf3 = sizeof(double3)*_n;
		int szf2 = sizeof(double2)*_n;
		int szf = sizeof(double)*_n;

		hostX = (double3*)malloc(szf3);
		if(!hostX) return NULLmalloc;

		hostV = (double3*)malloc(szf3);
		if(!hostV) return NULLmalloc;

		hostF = (double3*)malloc(szf3);
		if(!hostF) return NULLmalloc;

		hostOldX = (double3*)malloc(szf3);
		if(!hostOldX) return NULLmalloc;

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

		hostA = (double3*)malloc(szf3);
		if(!hostA) return NULLmalloc;

		hostE = (double2*)malloc(szf2);
		if(!hostE) return NULLmalloc;

		hostM = (double*)malloc(szf);
		if(!hostM) return NULLmalloc;

		if(useCuda){
			CUDA_CHECK(cudaMalloc((void**)&devX,szf3));
			CUDA_CHECK(cudaMalloc((void**)&devV,szf3));
			CUDA_CHECK(cudaMalloc((void**)&devA,szf3));
			CUDA_CHECK(cudaMalloc((void**)&devE,szf2));
			CUDA_CHECK(cudaMalloc((void**)&devM,szf));
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
    Fout << Fnum+1 << "\n"
         << R << "\n"
         << binOutF << "\n"
         << Ntot << "\n";
    Fout.close();
    sayLog("   }created\n");
    return 0;
}

int TSpace::SaveParams(string FileName)
{
	sayLog("   SaveParams(" + FileName + "){\n");
    int i,spForV = 14;
    ofstream Fout(FileName.c_str());
    if(!Fout){return CantOpenFile;}
    Fout.precision(PrintPres);

    for(i=0;i<ParamFHead.size();++i){ Fout << setw(spForV) << ParamFHead[i]; }
    Fout << '\n';
    sayLog("      text part was written\n");

    vector<double> outData;
    Fout << setw(spForV) << Ntot;

    outData = {k12, k6, R, sqrt(r_brd), endT, totalT, 1/dumpDT, 1/dt, Tmp, mu, kv, dissipK, TmpStabEps};
    for(i = 0; i < outData.size(); ++i) Fout << setw(spForV) << outData[i];
    Fout << setw(spForV) << TmpStabGap;
    Fout << setw(spForV) << (compMode-CompModeID)*(useCuda ? -1 : 1);
    Fout << setw(spForV) << binOutF;
    Fout << '\n';
    Fout.close();

    sayLog(string("      numeric part was written\n")+
    	   string("   }\n"));

    return 0;
}

int TSpace::SaveParticles(string FileName, bool bin, int mode)
{
	sayLog(string("      SaveStars (") + FileName + string("){\n"));
	int i, spForV = 14;
    ofstream Fout;

    if(bin){/*
        Fout.open(FileName.c_str(),ios::binary | ios::trunc);
        if(!Fout){ return CantCreateFile; }
        int buf = Ntot;
        Fout.write((char*)&(buf),sizeof(int));
        Fout.write((char*)&(Fnum),sizeof(int));
        Fout.write((char*)&(totalT),sizeof(double));
        */
    } else {
        Fout.open(FileName.c_str(),ios::trunc);
        if(!Fout){ return CantCreateFile; }
        Fout.precision(PrintPres);
        Fout << Ntot
        	 << "\nLattice=\"1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0\" Properties=pos:R:3:velo:R:3:mass:R:1 Time=" << toString(totalT) << "\n";
        // for(i = 0; i < GalaxyFHead.size(); ++i){ Fout << setw(spForV) << GalaxyFHead[i]; }
        // Fout << "\n" << setw(spForV*8) << Ntot << setw(spForV)  << Fnum << setw(spForV)  << totalT << "\n";
    }
    sayLog("         Header was written\n");
    if(mode == CompLoadMode) switch(compMode){
    	case LPcompMode:
    		for(i = 0; i < Ntot; ++i)	hostV[i] -= hostA[i]*(dt/2);
    		sayLog("         velocities were modified for saving\n");
    		break;
    }
    //Fout << scientific;
    for(i = 0; i < Ntot; ++i)
    {
    	if(dot(hostX[i]) == INFINITY) return StarXIsInf;
    	if(hostM[i] == INFINITY) return StarMIsInf;
    	if(dot(hostV[i]) == INFINITY) return StarVIsInf;
    	if((hostX[i].x == NAN) || (hostX[i].y == NAN) || (hostX[i].z == NAN) || (hostX[i].x == -NAN) || (hostX[i].y == -NAN) || (hostX[i].z == -NAN)) return StarXIsNan;
    	if(hostM[i] == NAN) return StarMIsNan;
    	if((hostV[i].x == NAN) || (hostV[i].y == NAN) || (hostV[i].z == NAN) || (hostV[i].x == -NAN) || (hostV[i].y == -NAN) || (hostV[i].z == -NAN)) return StarVIsNan;
        if(bin){
            Fout.write((char*)&(hostX[i].x),sizeof(double));
            Fout.write((char*)&(hostX[i].y),sizeof(double));
            Fout.write((char*)&(hostX[i].z),sizeof(double));
            Fout.write((char*)&(hostV[i].x),sizeof(double));
            Fout.write((char*)&(hostV[i].y),sizeof(double));
            Fout.write((char*)&(hostV[i].z),sizeof(double));
            Fout.write((char*)&(hostM[i]),sizeof(double));
        } else {
        	Fout << setw(spForV) << hostX[i].x
        		 << setw(spForV) << hostX[i].y
        		 << setw(spForV) << hostX[i].z
        		 << setw(spForV) << hostV[i].x
        		 << setw(spForV) << hostV[i].y
        		 << setw(spForV) << hostV[i].z
        		 << setw(spForV) << hostM[i]
        		 << "\n";
        }
    }
    Fout.close();
    sayLog("         Stars were saved\n");
    if(mode == CompLoadMode) switch(compMode){
    	case LPcompMode:
    		for(i = 0; i < Ntot; ++i)	hostV[i] += hostA[i]*(dt/2);
    		sayLog("         velocities were restored\n");
    		break;
    }
    sayLog("      }(" + FileName + ")\n");
    return 0;
}

int TSpace::checkInput(int _n, double _dumpDT, double _dt)
{
    if(compMode < CompModeID){ compMode += CompModeID; }

    if((compMode != VRcompMode) &&
       (compMode != LPcompMode) &&
       (compMode != RKcompMode) &&
       (compMode != AdVRcompMode) &&
       (compMode != MYcompMode))
      { return WrongCompMode; }
    if(_n == 1){ sayLog("\nN == 1\n",1); }
    if(_n <= 0){ return NLessOrEq0; }
    if(_n > MaxN){ return NisTooBig; }
    if(useCuda && (_n % BlockW != 0)){ return CUDA_WRONG_NBlockW; }
    if(a <= 0){ return ALessOrEq0; }
    if(a > R){ sayLog("\na > R\n",1); }
    if(a == INFINITY){ return AIsInf; }
    if(_dt <= 0){ return dtLessOrEq0; }
    if(_dt == INFINITY){ return dtIsInf; }
    if(_dumpDT <= 0){ return dumpDTLessOrEq0; }
    if(_dumpDT == INFINITY){ return dumpDTIsInf; }
    if((k12 == 0) || (k6 == 0)){ return Geq0; }
    if((k12 < 0) || (k6 < 0)){ sayLog("\nk12 or k6 < 0\n",1); }
    if((k12 == INFINITY) || (k6 == INFINITY)){ return GIsInf; }
    if(endT <= 0){ return TLessOrEq0; }
    if(endT + (1/_dumpDT)*(1 + SYS_EPS) < totalT){ return TlessCalcedT; }
    if(endT == INFINITY){ return TIsInf; }
    //if(dissipK < 0){ return dissipKless0; }
    if(std::abs(dissipK) < SYS_EPS){ sayLog("\n|dissipK| < " + toString(SYS_EPS) + "\n",1); }
    if(r_brd < 0){ return RbrdLess0; }
    if(Tmp < 0){ return TmpLess0; } else
    if(Tmp < SYS_EPS){ sayLog("\ninitial temperature < " + toString(SYS_EPS) + "\n",1); }
    if(mu <= 0){ return muLessEq0; }
    if(kv < 0){ return kvLess0; }
    if(TmpStabEps <= 0){ return TmpStabEps_LessEq0; }
    if(TmpStabEps < SYS_EPS){ sayLog("\nTmpStabEps < " + toString(SYS_EPS) + "\n",1); }
    if(TmpStabGap < 0){ return TmpStabGap_Less0; }
    if(TmpStabGap > 100000){ sayLog("\nTmpStabGap > 100000\n",1); }
    if(Ntot/pow3(2*R) > 0.7 + SYS_EPS){
    	sayLog("\nn > 0.7\n",1);
    	if(pow3(2*R) <= Ntot){
    		sayLog("\nNtot/(2R)^3 = n >= 1\n",1);
    		return TooDenseSystem;
    	}
    }

    if(binOutF == 1){
    	sayLog("\nbin files aren't supported for now\n",1);
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
	int err = checkInput(_n, _dumpDT, _dt);
	if(err) return err;
	err = ResizeStars(_n);
	if(err) return err;

	nCudaB = (Ntot+BlockW-1)/BlockW;
	dumpDT = 1/_dumpDT;
	dt = 1/_dt;
	r_brd = r_brd*r_brd;

	sayLog("      }\n");
	return 0;
}

int TSpace::SafeLoadParams(string *goodName)
{
    string b_s = "./" + GlxName + "_" + PrmFName;
    bool b_b = access(b_s.c_str(),0);
	if(b_b){
		sayLog(string("   '" + b_s + "' not found\n") +
		       string("   trying to find './"+PrmFName+"'\n"));
		cout << "'" << b_s << "' not found\n"
			 << "trying to find './" << PrmFName << "' ... ";
		b_s = "./"+PrmFName;
	}
	*goodName = b_s;
	int err = LoadParams(b_s);
    if(b_b) cout << (err ? "failed\n" : "success\n");
    if(err) return err;
    return 0;
}

int TSpace::LoadParams(string FileName)
{
    int _n, err;
    double _dumpDT, _dt;
    string buf_s;
    sayLog("   LoadParams (" + FileName + ") {\n");

    ifstream Fin(FileName.c_str());
    if(!Fin){ return CantOpenFile; }
    std::getline(Fin, buf_s);
    Fin >> _n >> k12 >> k6
    	>> R >> r_brd// >> a
    	>> endT >> totalT >> _dumpDT >> _dt
    	// endT - total time to compute from 0 to endT
    	// totalT - total time already computed
    	// dumpDT - time between saving to file
    	// dt - time step for computation
    	>> Tmp >> mu >> kv >> dissipK >> TmpStabEps
    	>> TmpStabGap >> compMode >> binOutF;
    a = 1;
    Fin.close();

    sayLog("      Parameters were loaded\n");

    err = applyInput(_n, _dumpDT, _dt);
    if(err) return err;

    sayLog("   }\n   " + string(useCuda ? "GPU" : "CPU") + " in use\n");
    return 0;
}

int TSpace::LoadParticles(string FileName, bool bin, int mode)
{
    int i,N,err;
    string buf;
    ifstream Fin;
    sayLog("   LoadParticles (" + FileName + ") {\n");
    if(bin){
        Fin.open(FileName.c_str(),ios::binary);
        if(!Fin){ return CantOpenFile; }
        Fin.read((char*)&(N),sizeof(int));
        Fin.read((char*)&(Fnum),sizeof(int));
        Fin.read((char*)&(totalT0),sizeof(double));
    } else {
        Fin.open(FileName.c_str());
        if(!Fin){ return CantOpenFile; }
        //Fin >> N >> Fnum >> totalT0;
        Fin >> N;
        Fnum = 0;
        totalT0 = 0;
        getline(Fin,buf);
        getline(Fin,buf);
    }
    if(N <= 0) return NLessOrEq0;
    if(Fnum > 3000000000) return FnumLess0;
    if(totalT0 < 0) return CalcedT0Less0;
    totalT = totalT0;
    err=ResizeStars(N);
    if(err) return err;
    sayLog("      Head data was loaded\n");
    for(i = 0; i < Ntot; ++i)
    {
        if(bin){
            Fin.read((char*)&(hostX[i].x),sizeof(double));
            Fin.read((char*)&(hostX[i].y),sizeof(double));
            Fin.read((char*)&(hostX[i].z),sizeof(double));
            Fin.read((char*)&(hostV[i].x),sizeof(double));
            Fin.read((char*)&(hostV[i].y),sizeof(double));
            Fin.read((char*)&(hostV[i].z),sizeof(double));
            Fin.read((char*)&(hostM[i]),sizeof(double));
        } else {
            Fin >> hostX[i].x >> hostX[i].y >> hostX[i].z >> hostV[i].x >> hostV[i].y >> hostV[i].z >> hostM[i];
        }
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
		CUDA_CHECK(cudaMemcpy(devM,hostM,Ntot*sizeof(double),cudaMemcpyHostToDevice));
		sayLog("      devX and devM were filled with host data\n");
	}

	if(mode == CompLoadMode){
		if(useCuda){
			cudaFindAllA<<< nCudaB,BlockW >>>(devX, devA, hostM[0], R, r_brd, Ntot);
    		CUDA_CHECK(cudaThreadSynchronize());
    		CUDA_CHECK(cudaMemcpy(hostA,devA,Ntot*sizeof(double3),cudaMemcpyDeviceToHost));
		} else {
			findAllA();
		}
		sayLog("      Acceleration were computed\n");
		switch(compMode){
		case VRcompMode:
			for(i = 0; i < Ntot; ++i) hostOldX[i] = shiftR(hostX[i] - hostV[i]*dt + hostA[i]*(dt*dt/2), R);
			sayLog("      oldCrd for VRmode were computed\n");
			break;
		case LPcompMode:
			for(i = 0; i < Ntot; ++i) hostV[i] += hostA[i]*(dt/2);
			sayLog("      v(1/2dt) for LPmode were computed\n");
			break;
		}
	}

	return 0;
}

int TSpace::stabilizeTmp(void)
{
	sayLog("   stabTmp{\n");

	int b_i;
	long int print_step, k_step = 0;
	double TmpAv;
	double t = 0;
    time_t rStime;

    b_i = postGetProc(CompLoadMode);
    if(b_i) return b_i;
    TmpAv = Tmp;
    print_step = Ntot > 10000 ? 1 : 100000000/(Ntot*Ntot);

    time(&rStime);
    do
    {
    	b_i = useCuda ? cudaDoTimeStep(compMode) : doTimeStep(compMode);
    	if(b_i) return b_i;
    	++k_step;
    	t += dt;
    	computeEk();
    	TmpAv = (TmpAv*(TmpStabGap-1) + Tmp_curr)/TmpStabGap;

        if(k_step%print_step == 0)
        {
        	time_progress(rStime, time(0), 1.0/(1 + log(epsDlt(TmpAv, Tmp)/TmpStabEps)/dissipK / t * pow2(time(0) - rStime)), "Tmp stabilizing");
        }
    }while(!(almostEq(TmpAv, Tmp, TmpStabEps) && almostEq(Tmp_curr, Tmp, TmpStabEps*2) && (k_step > TmpStabGap)));
    Tmp = Tmp_curr;

    sayLog("Tmp_res = " + toString(Tmp_curr) + "; Tmp_av = " + toString(TmpAv) + "                        \n", 1);
    sayLog("   }\n");
    return 0;
}

int TSpace::main_compute(void)
{
	sayLog("   main_compute{\n");
	if(!access(StpFName.c_str(),0)) return StopFileAlreadyExists;

    int b_i = 0;
    time_t rStime;
    double kt;
    bool stpCalc = 0;
    string buf_s, BaseName = "./" + GlxName + "/";
    string FramePath = BaseName + FramesFilesFolder + "/";

    //ofstream FoutT((BaseName+TFile).c_str(),ios::app);
    // This was to make possible to continue computation, but here it's not possible for other reasons.
    // so no need to keep this bothering reature
    ofstream FoutT((BaseName + TFile).c_str());
    if(!FoutT) return CantCreateFile;
    FoutT.precision(PrintPres);

    if(totalT == 0) FoutT << totalT0 << endl;
    SaveParticles(FramePath + toString(Fnum) + string(".") + ParticleFileExt, binOutF, CompLoadMode);
    cout << "Galaxy name : '" << GlxName << "'\n"
    	 << "Done start part = " << 100*totalT0/endT << "%\n";
    /*
     * kt is for accounting that after calc will finish then some more time will be needed for post_proc E
     * */

    kt = dt/dumpDT;
    switch(compMode){
		case MYcompMode:
			buf_s = "MY";
			break;
    	case VRcompMode:
    		buf_s = "VR";
    		break;
    	case LPcompMode:
    		buf_s = "LP";
    		break;
    	case RKcompMode:
    		kt /= 4;
    		buf_s = "RK";
    		break;
    	case AdVRcompMode:
    		kt /= 4;
    		buf_s = "AdVR";
    		break;
    }
    kt = 1 + kt * 0.788; // TODO find k - it's not 0.788 here

    sayLog(string("      "+buf_s+" computation mode\n")+
    	   string("      Starting the computation{\n"));
    time(&rStime);
    do
    {
    	b_i = useCuda ? cudaDoTimeStep(compMode) : doTimeStep(compMode);
    	//b_i = doTimeStep(compMode);
    	if(b_i) return b_i;
    	localT += dt;
    	totalT += dt;

        if(localT > dumpDT*(1 - SYS_EPS))
        {
        	sayLog("      Saving " + toString(++Fnum) + " file...\n");
            localT = 0;
            FoutT << totalT << endl;

            time_progress(rStime, time(0), totalT / endT, "computing");

            b_i = SaveParticles(FramePath + toString(Fnum) + string(".") + ParticleFileExt, binOutF, CompLoadMode);
            if(b_i) return b_i;
            sayLog("      ... saved\n");
            stpCalc = !((totalT < endT-dt*(1 - SYS_EPS)) && access(StpFName.c_str(),0));
        }
    }while(!stpCalc);
    FoutT.close();

    double usedT = totalT - totalT0;
    eff = usedT/dt * Ntot*Ntot / rt*3600;
    //if(!useCuda) eff/=2;
    buf_s = totalT < endT - dt * (1 - SYS_EPS) ? "terminated" : "done";
    buf_s = "Computation was " + buf_s;
    cout << buf_s << "\n";
    sayLog(string("      "+buf_s+"\n")+
    	   string("      }\n")+
    	   string("   }\n")+
    	   string("   Efficiency{\n")+
    	   string("      time computed = "+toString(usedT)+"; dt = "+toString(dt)+"; usedT/dt = "+toString(usedT/dt)+"\n")+
           string("      N = " + toString(Ntot) + "; real_t [h] = " + toString(1.0/3600*rt) + "\n")+
           string("      efficiency (e = endT/dt*N^2/t_real) = " + toString(eff) + "\n")+
           string("   }\n"));


    return 0;
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
    	case MYcompMode: MYmove(); break;
        case VRcompMode: VRmove(); break;
        case LPcompMode: LPmove(); break;
        case RKcompMode: RKmove(); break;
        case AdVRcompMode: AdVRmove(); break;
        default: return WrongCompMode;
    }

    //shiftAll();

    // this is necessary here - without it there is a slow drift of CM
    stopMacroMovement();

    return 0;
}

double3 TSpace::findA(double3 crd1, unsigned long N)
{
    unsigned long i;
    double3 bv, sumf;
    double r2, ar2, ar6;

    sumf.x = sumf.y = sumf.z = 0;
    for(i = 0; i < N; ++i)
    {
        bv = hostX[i] - crd1;
        r2 = dot(bv);
        ar2 = a*a/r2;
        ar6 = ar2 * ar2 * ar2;
        sumf += bv * 6*ar6/r2*(2*ar6 - 1); // F = -r/abs(r) * phi'
    }
    for(i=N+1;i<Ntot;++i)
    {
        bv = hostX[i] - crd1;
        r2 = dot(bv);
        ar2 = a*a/r2;
        ar6 = ar2 * ar2 * ar2;
        sumf += bv * 6*ar6/r2*(2*ar6 - 1); // F = -r/abs(r) * phi'
    }

    return sumf;
}


int TSpace::findAllA(void)
{
    unsigned long i,j;
    double3 bv;
    double ar2, ar6;
    double D = 2*R;
    double m_u = 24/hostM[0];
/*
    for(i = 0; i < Ntot; ++i){ hostF[i].x = hostF[i].y = hostF[i].z = 0; }
    for(i = 0; i < Ntot; ++i){
        for(j = i+1; j < Ntot; ++j){
        	 bv = shiftR(hostX[j] - hostX[i], R, D);
        	 ar2 = 1/dot(bv);
        	 ar6 = ar2*ar2*ar2;
        	 bv *= 24*ar6*ar2*(1 - 2*ar6);
        	 hostF[i] += bv;
        	 hostF[j] -= bv;
        }
    }
    for(i = 0; i < Ntot; ++i){ hostA[i]  = hostF[i]/hostM[i]; }
*/
   for(i = 0; i < Ntot; ++i){ hostA[i].x = hostA[i].y = hostA[i].z = 0; }
   // this openmp is checked 1 time - it has no effect on accuracy
    #ifdef _OPENMP
    	#pragma omp parallel for private(bv, ar2, ar6, j)
	#else
		#warning "OpenMP unused"
    #endif
    for(i = 0; i < Ntot; ++i){
    	//m_u = 24/hostM[i];
        for(j = 0; j < i; ++j){
        	 bv = shiftR(hostX[j] - hostX[i], R, D);
        	 ar2 = dot(bv);
        	 if(ar2 < r_brd){
        		 ar2 = 1/ar2;
        		 ar6 = ar2*ar2*ar2;
        		 bv *= m_u*ar6*ar2*(1 - 2*ar6);
        		 hostA[i] += bv;
        	 }
        }
        for(j = i+1; j < Ntot; ++j){
        	 bv = shiftR(hostX[j] - hostX[i], R, D);
        	 ar2 = dot(bv);
        	 if(ar2 < r_brd){
        		 ar2 = 1/ar2;
        		 ar6 = ar2*ar2*ar2;
        		 bv *= m_u*ar6*ar2*(1 - 2*ar6);
        		 hostA[i] += bv;
        	 }
        }
    }

    return 0;
}

int TSpace::useThermostat(void)
{
	double b_d;
	int i;

    if(thermostatOn){
    	b_d = sqrt(2*Tmp*dissipK / (hostM[0] * dt));
#ifdef _OPENMP
	#pragma omp parallel for
#else
	#warning "OpenMP unused"
#endif

    	for(i = 0; i < Ntot; ++i){
    		hostA[i] += (hostV[i] * (-dissipK)) + make_double3(gaussRand() * b_d, gaussRand() * b_d, gaussRand() * b_d);
    	}
    }

    return 0;
}

int TSpace::MYmove(void)
{
    int i;
    double D = 2*R;

    findAllA();
    useThermostat();
#ifdef _OPENMP
	#pragma omp parallel for
#else
	#warning "OpenMP unused"
#endif
    for(i = 0; i < Ntot; ++i){
    	hostV[i] += hostA[i] * dt;
    	hostX[i] = shiftR(hostX[i] + hostV[i] * dt, R, D);
    }

    return 0;
}

int TSpace::VRmove(void)
{
    int i;
    double3 bv;
    double D = 2*R;

    findAllA();
    useThermostat();
#ifdef _OPENMP
	#pragma omp parallel for private(bv)
#else
	#warning "OpenMP unused"
#endif
    for(i = 0; i < Ntot; ++i){
        bv = hostX[i];
        hostX[i] = shiftR(hostX[i]*2 - hostOldX[i] + hostA[i]*(dt*dt), R, D);
        hostV[i] = shiftR(hostX[i] - hostOldX[i], R, D)/(2*dt);
        hostOldX[i] = bv;
    }
    return 0;
}

int TSpace::LPmove(void)
{
    int i;
    double D = 2*R;

    // V is 1/2 time-step ahead of X
    for(i = 0; i < Ntot; ++i) hostX[i] = shiftR(hostX[i] + hostV[i] * dt, R, D); // finish previous step by modifying X
    findAllA();                                        // compute new A
    useThermostat();
    for(i = 0; i < Ntot; ++i) hostV[i] += hostA[i]*dt; // compute new V

    return 0;
}

int TSpace::RKmove(void)
{
    int i;

    findAllA();
#ifdef _OPENMP
	#pragma omp parallel for
#else
	#warning "OpenMP unused"
#endif
    for(i=0;i<Ntot;++i){
    	hostOldX[i] = hostX[i];

    	hostVK[0][i] = hostA[i];
    	hostRK[0][i] = hostV[i];
    	hostX[i] = hostOldX[i] + hostRK[0][i]*(dt/2);
    }

    findAllA();
#ifdef _OPENMP
	#pragma omp parallel for
#else
	#warning "OpenMP unused"
#endif
    for(i=0;i<Ntot;++i){
    	hostVK[1][i] = hostA[i];
    	hostRK[1][i] = hostV[i] + hostVK[0][i]*(dt/2);
    	hostX[i] = hostOldX[i] + hostRK[1][i]*(dt/2);
    }

    findAllA();
#ifdef _OPENMP
	#pragma omp parallel for
#else
	#warning "OpenMP unused"
#endif
    for(i=0;i<Ntot;++i){
    	hostVK[2][i] = hostA[i];
    	hostRK[2][i] = hostV[i] + hostVK[1][i]*(dt/2);
    	hostX[i] = hostOldX[i] + hostRK[2][i]*dt;
    }

    findAllA();

#ifdef _OPENMP
	#pragma omp parallel for
#else
	#warning "OpenMP unused"
#endif
	for(i = 0; i < Ntot; ++i){
    	hostVK[3][i] = hostA[i];
    	hostRK[3][i] = hostV[i] + hostVK[2][i]*dt;
		hostA[i] = (hostVK[0][i] + 2*hostVK[1][i] + 2*hostVK[2][i] + hostVK[3][i])/6;
		hostX[i] = hostOldX[i] + (hostRK[0][i] + 2*hostRK[1][i] + 2*hostRK[2][i] + hostRK[3][i])*(dt/6);
	}
	useThermostat();

    for(i = 0; i < Ntot; ++i){
    	hostV[i] += hostA[i]*dt;
    }

    return 0;
}

int TSpace::AdVRmove(void)
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

	for(i = 0; i < Nslice; ++i){
		gCondFnc[i] = 0;
	}

	for(i = 0; i < Ntot; ++i){
		for(j = 0; j < i; ++j){
			k = floor( length( shiftR(hostX[j] - hostX[i], R) )/dr );
			if(k < Nslice)
				++gCondFnc[k];
		}
		// every particle except i0
		for(i = i+1; i < Ntot; ++i){
			k = floor( length( shiftR(hostX[j] - hostX[i], R) )/dr );
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
		E.x += hostM[i]*dot(hostV[i])/2;
	}
	if(std::abs(E.x) == INFINITY) return EkIsInf;
	if((E.x == NAN) || (E.x == -NAN)) return EkIsNan;

	Tmp_curr = 2.0/3 * E.x/Ntot;

	return 0;
}

int TSpace::computeEp(void)
{
	int i, j;
	double3 r;
	double ar2, ar6;
	double D = 2*R;

	E.y = 0;
	for(i = 0; i < Ntot; ++i){
		hostE[i].y = 0;
	}

// this openmp is checked - it has no influence on accuracy
#ifdef _OPENMP
	#pragma omp parallel for private(ar2, ar6, j, r)
#else
	#warning "OpenMP unused"
#endif
	for(i = 0; i < Ntot; ++i) for(j = i+1; j < Ntot; ++j){
		r = shiftR(hostX[j] - hostX[i], R, D);
		ar2 = dot(r,r);
		if(ar2 < r_brd){
			ar2 = 1/ar2;
			ar6 = ar2 * ar2 * ar2; // (a/r)^6
			hostE[i].y += 4*ar6*(ar6 - 1);
		}
	}
	for(i = 0; i < Ntot; ++i) E.y += hostE[i].y;
	if(std::abs(E.y) == INFINITY) return EpIsInf;
	if((E.y == NAN) || (E.y == -NAN)) return EpIsNan;

	return 0;
}

int TSpace::doEstep(void)
{
	int err_handl = computeEk();
	if(err_handl) return err_handl;
	err_handl = computeEp();
	if(err_handl) return err_handl;

    return 0;
}

int TSpace::cudaDoEstep(void)
{
	int i;
	CUDA_CHECK(cudaMemcpy(devV,hostV,Ntot*sizeof(double3),cudaMemcpyHostToDevice));
	cudaFindE<<< nCudaB,BlockW >>>(devX,devV,devE,hostM[0],R,r_brd,Ntot);
	CUDA_CHECK(cudaThreadSynchronize());
	CUDA_CHECK(cudaMemcpy(hostE,devE,nCudaB*sizeof(double2),cudaMemcpyDeviceToHost));
	E.x=E.y=0;
	for(i = 0; i < nCudaB; ++i) E += hostE[i];
	E.y/=2;

	if(E.x == INFINITY) return EkIsInf;
	if(E.x == NAN) return EkIsNan;
	if(std::abs(E.y) == INFINITY) return EpIsInf;
	if((E.y == NAN) || (E.y == -NAN)) return EpIsNan;
    return 0;
}

__global__ void cudaFindE(double3 *devX, double3 *devV, double2 *devE, double m, double R, double r_brd, int Ntot)
{
	__shared__ double3 c1[BlockW],c2[BlockW],v1[BlockW];
	__shared__ double2 e1[BlockW];
	//__shared__ double m1[BlockW],m2[BlockW];
	int gind = BlockW*blockIdx.x + threadIdx.x;

	if(gind >= Ntot) return;

	int tind,tile,i,tx=threadIdx.x;
	double3 r;
	double ar2, ar6, m_u = m/2, D = 2*R;
	c1[tx].x = devX[gind].x;
	c1[tx].y = devX[gind].y;
	c1[tx].z = devX[gind].z;
	v1[tx].x = devV[gind].x;
	v1[tx].y = devV[gind].y;
	v1[tx].z = devV[gind].z;
	//m1[tx] = devM[gind];
	//e1[tx].x = m1[tx]*dot(v1[tx])/2;
	e1[tx].x = m_u*dot(v1[tx]);
	e1[tx].y = 0;

	for(tile=0; tile*BlockW + tx < Ntot; ++tile){
		tind = tile*BlockW+tx;
		c2[tx].x = devX[tind].x;
		c2[tx].y = devX[tind].y;
		c2[tx].z = devX[tind].z;
		//m2[tx] = devM[tind];
		__syncthreads();
		for(i=0; i<BlockW; ++i){
			if(tile*BlockW + i != gind){
				r = shiftR(c2[i] - c1[tx], R, D);
				ar2 = dot(r);
				if(ar2 < r_brd){
					ar2 = 1/ar2;
					ar6 = ar2 * ar2 * ar2; // (a/r)^6
					e1[tx].y += 4*ar6*(ar6 - 1);
				}
			}
		}
		__syncthreads();
	}

	i=BlockW/2;
	while(i>0){
		if(tx<i){
			tind = tx+i;
			e1[tx].x += e1[tind].x;
			e1[tx].y += e1[tind].y;
		}
		__syncthreads();
		i/=2;
		/*
		*	For correct work Ntot%BlockW must be ==0, because if it's not, then
		*	/=2 won't properly wrap indexes
		*/
	}

	if(tx==0){
		devE[blockIdx.x].x = e1[0].x;
		devE[blockIdx.x].y = e1[0].y;
	}
}
/*
 */
int TSpace::cudaDoTimeStep(int cmpMode){
	int err;
	switch(cmpMode){
		case MYcompMode: err = cudaMYmove(); break;
		case VRcompMode: err = cudaVRmove(); break;
		case LPcompMode: err = cudaLPmove(); break;
		case RKcompMode: err = cudaRKmove(); break;
		case AdVRcompMode: err = cudaAdVRmove(); break;
		default: err = WrongCompMode;
	}

    //shiftAll();
    stopMacroMovement();

	return err;
}

int TSpace::cudaMYmove(void)
{
	CUDA_CHECK(cudaMemcpy(devX, hostX, Ntot*sizeof(double3), cudaMemcpyHostToDevice));
	cudaFindAllA<<< nCudaB,BlockW >>>(devX, devA, hostM[0], R, r_brd, Ntot);
	CUDA_CHECK(cudaThreadSynchronize());
	CUDA_CHECK(cudaMemcpy(hostA, devA, Ntot*sizeof(double3), cudaMemcpyDeviceToHost));
	useThermostat();

    int i;

#ifdef _OPENMP
	#pragma omp parallel for
#else
	#warning "OpenMP unused"
#endif
    for(i = 0; i < Ntot; ++i){
    	hostV[i] += hostA[i] * dt;
    	hostX[i] = shiftR(hostX[i] + hostV[i] * dt, R);
    }

    return 0;
}

int TSpace::cudaVRmove(void)
{
	CUDA_CHECK(cudaMemcpy(devX, hostX, Ntot*sizeof(double3), cudaMemcpyHostToDevice));
	cudaFindAllA<<< nCudaB,BlockW >>>(devX, devA, hostM[0], R, r_brd, Ntot);
	CUDA_CHECK(cudaThreadSynchronize());
	CUDA_CHECK(cudaMemcpy(hostA, devA, Ntot*sizeof(double3), cudaMemcpyDeviceToHost));
	useThermostat();

    int i;
    double3 bv;

#ifdef _OPENMP
	#pragma omp parallel for private(bv)
#else
	#warning "OpenMP unused"
#endif
    for(i = 0; i < Ntot; ++i){
        bv = hostX[i];
        hostX[i] = shiftR(hostX[i]*2 - hostOldX[i] + hostA[i]*(dt*dt), R);
        hostV[i] = shiftR(hostX[i] - hostOldX[i], R)/(2*dt);
        hostOldX[i] = bv;
    }

    return 0;
}

int TSpace::cudaLPmove(void)
{
    int i;
    double D = 2*R;

    // V is 1/2 time-step ahead of X

    for(i = 0; i < Ntot; ++i) hostX[i] = shiftR(hostX[i] + hostV[i] * dt, R, D); // finish previous step by modifying X
	CUDA_CHECK(cudaMemcpy(devX, hostX, Ntot*sizeof(double3), cudaMemcpyHostToDevice));
	cudaFindAllA<<< nCudaB,BlockW >>>(devX, devA, hostM[0], R, r_brd, Ntot); // compute new A
	CUDA_CHECK(cudaThreadSynchronize());
	CUDA_CHECK(cudaMemcpy(hostA, devA, Ntot*sizeof(double3), cudaMemcpyDeviceToHost));
	useThermostat();

    for(i = 0; i < Ntot; ++i) hostV[i] += hostA[i]*dt; // compute new V

    return 0;
}

int TSpace::cudaRKmove(void)
{
    int i;

	CUDA_CHECK(cudaMemcpy(devX, hostX, Ntot*sizeof(double3), cudaMemcpyHostToDevice));
	cudaFindAllA<<< nCudaB,BlockW >>>(devX, devA, hostM[0], R, r_brd, Ntot); // compute new A
	CUDA_CHECK(cudaThreadSynchronize());
	CUDA_CHECK(cudaMemcpy(hostA, devA, Ntot*sizeof(double3), cudaMemcpyDeviceToHost));
#ifdef _OPENMP
	#pragma omp parallel for
#else
	#warning "OpenMP unused"
#endif
    for(i = 0; i < Ntot; ++i){
    	hostOldX[i] = hostX[i];
    	hostVK[0][i] = hostA[i];
    	hostRK[0][i] = hostV[i];
    	hostX[i] = hostOldX[i] + hostRK[0][i]*(dt/2);
    }

	CUDA_CHECK(cudaMemcpy(devX, hostX, Ntot*sizeof(double3), cudaMemcpyHostToDevice));
	cudaFindAllA<<< nCudaB,BlockW >>>(devX, devA, hostM[0], R, r_brd, Ntot); // compute new A
	CUDA_CHECK(cudaThreadSynchronize());
	CUDA_CHECK(cudaMemcpy(hostA, devA, Ntot*sizeof(double3), cudaMemcpyDeviceToHost));
#ifdef _OPENMP
	#pragma omp parallel for
#else
	#warning "OpenMP unused"
#endif
    for(i = 0; i < Ntot; ++i){
    	hostVK[1][i] = hostA[i];
    	hostRK[1][i] = hostV[i] + hostVK[0][i]*(dt/2);
    	hostX[i] = hostOldX[i] + hostRK[1][i]*(dt/2);
    }

	CUDA_CHECK(cudaMemcpy(devX, hostX, Ntot*sizeof(double3), cudaMemcpyHostToDevice));
	cudaFindAllA<<< nCudaB,BlockW >>>(devX, devA, hostM[0], R, r_brd, Ntot); // compute new A
	CUDA_CHECK(cudaThreadSynchronize());
	CUDA_CHECK(cudaMemcpy(hostA, devA, Ntot*sizeof(double3), cudaMemcpyDeviceToHost));
#ifdef _OPENMP
	#pragma omp parallel for
#else
	#warning "OpenMP unused"
#endif
    for(i = 0; i < Ntot; ++i){
    	hostVK[2][i] = hostA[i];
    	hostRK[2][i] = hostV[i] + hostVK[1][i]*(dt/2);
    	hostX[i] = hostOldX[i] + hostRK[2][i]*dt;
    }

	CUDA_CHECK(cudaMemcpy(devX, hostX, Ntot*sizeof(double3), cudaMemcpyHostToDevice));
	cudaFindAllA<<< nCudaB,BlockW >>>(devX, devA, hostM[0], R, r_brd, Ntot); // compute new A
	CUDA_CHECK(cudaThreadSynchronize());
	CUDA_CHECK(cudaMemcpy(hostA, devA, Ntot*sizeof(double3), cudaMemcpyDeviceToHost));

#ifdef _OPENMP
	#pragma omp parallel for
#else
	#warning "OpenMP unused"
#endif
	for(i = 0; i < Ntot; ++i){
    	hostVK[3][i] = hostA[i];
    	hostRK[3][i] = hostV[i] + hostVK[2][i]*dt;
		hostA[i] = (hostVK[0][i] + 2*hostVK[1][i] + 2*hostVK[2][i] + hostVK[3][i])/6;
		hostX[i] = hostOldX[i] + (hostRK[0][i] + 2*hostRK[1][i] + 2*hostRK[2][i] + hostRK[3][i])*(dt/6);
	}
	useThermostat();

    for(i = 0; i < Ntot; ++i){
    	hostV[i] += hostA[i]*dt;
    }

    return 0;
}

int TSpace::cudaAdVRmove(void)
{
    int i;
    double k1 = 0.1786178958448091;
    double k2 = -0.2123418310626054;
    double k3 = -0.0662645826698185;

    for(i=0;i<Ntot;++i) hostX[i] += hostV[i]*(k1*dt);
	CUDA_CHECK(cudaMemcpy(devX,hostX,Ntot*sizeof(double3),cudaMemcpyHostToDevice));
	//cudaFindAllA<<< nCudaB,BlockW >>>(devX, devA, hostM[0], R, Ntot);
	CUDA_CHECK(cudaThreadSynchronize());
	CUDA_CHECK(cudaMemcpy(hostA,devA,Ntot*sizeof(double3),cudaMemcpyDeviceToHost));
    for(i=0;i<Ntot;++i){
    	hostV[i] += hostA[i]*((0.5-k2)*dt);
    	hostX[i] += hostV[i]*(dt*k3);
    }
	CUDA_CHECK(cudaMemcpy(devX,hostX,Ntot*sizeof(double3),cudaMemcpyHostToDevice));
	//cudaFindAllA<<< nCudaB,BlockW >>>(devX, devA, hostM[0], R, Ntot);
	CUDA_CHECK(cudaThreadSynchronize());
	CUDA_CHECK(cudaMemcpy(hostA,devA,Ntot*sizeof(double3),cudaMemcpyDeviceToHost));
    for(i=0;i<Ntot;++i){
    	hostV[i] += hostA[i]*(k2*dt);
    	hostX[i] += hostV[i]*((1-2*(k1+k3))*dt);
    }
	CUDA_CHECK(cudaMemcpy(devX,hostX,Ntot*sizeof(double3),cudaMemcpyHostToDevice));
	//cudaFindAllA<<< nCudaB,BlockW >>>(devX, devA, hostM[0], R, Ntot);
	CUDA_CHECK(cudaThreadSynchronize());
	CUDA_CHECK(cudaMemcpy(hostA,devA,Ntot*sizeof(double3),cudaMemcpyDeviceToHost));
    for(i=0;i<Ntot;++i){
    	hostV[i] += hostA[i]*(k2*dt);
    	hostX[i] += hostV[i]*(k3*dt);
    }
	CUDA_CHECK(cudaMemcpy(devX,hostX,Ntot*sizeof(double3),cudaMemcpyHostToDevice));
	//cudaFindAllA<<< nCudaB,BlockW >>>(devX, devA, hostM[0], R, Ntot);
	CUDA_CHECK(cudaThreadSynchronize());
	CUDA_CHECK(cudaMemcpy(hostA,devA,Ntot*sizeof(double3),cudaMemcpyDeviceToHost));
    for(i=0;i<Ntot;++i){
    	hostV[i] += hostA[i]*((0.5-k2)*dt);
    	hostX[i] += hostV[i]*(k1*dt);
    }

    return 0;
}


__global__ void cudaFindAllA(double3 *devX, double3 *devA, double m, double R, double r_brd, int Ntot)
{
	__shared__ double3 c1[BlockW],c2[BlockW],a1[BlockW];
	//__shared__ double m2[BlockW];
	int gind = BlockW*blockIdx.x + threadIdx.x; // global index

	if(gind >= Ntot) return;
	/*
	*	if(Ntot%BlockW != 0) then some threads in one block will never reach
	*	__syncthreads() command because of this "if" statement, and it's bad. So Ntot%BlockW must be == 0
	*/

	// set start values
	int tind, tile, i, tx = threadIdx.x;
	double3 r;
	double ar2, ar6, D = 2*R;
//	double m_u = 24.0/devM[tind];
	double m_u = 24.0/m;

	// copy "left column" of stars from global - the ones I'll sum TO
	c1[tx] = devX[gind];
	a1[tx].x = a1[tx].y = a1[tx].z = 0;

	for(tile = 0; tile*BlockW + tx < Ntot; ++tile){
		tind = tile*BlockW + tx;
		// copy part of "up row" from global - the 2nd part of current tile - the ones I'll sum WITH
		c2[tx] = devX[tind];
		// m2[tx] = devM[tind];
		__syncthreads();

		for(i = 0; i < BlockW; ++i){
			if(tile*BlockW + i != gind){
				r = shiftR(c2[i] - c1[tx], R, D); // this line takes ~60% of computational time of the whole N^2 block
				ar2 = dot(r);
				if(ar2 < r_brd){
					ar2 = 1/ar2;
					ar6 = ar2 * ar2 * ar2;
					r *= m_u*ar6*ar2*(1 - 2*ar6);
					a1[tx] += r;
				}
			}
		}

		__syncthreads();
	}

	devA[gind] = a1[tx];
}
