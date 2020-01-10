#include "format.cuh"

template<typename T> void stp(T str)
{
	cout << str;
	cin.get();
}

void print_sys_info(ostream &out)
{
	out << "sys info:\n"
		<< "eps = " << SYS_EPS << "\n"
		<< "CUDA_BlockW = " << BlockW << "\n";
}

void time_progress(double gt, double real_start_t, double curr_t, double done_part, string proc_name, int extra_strN)
{
	double real_t = curr_t - real_start_t;
	double left_t = done_part > 0 ? real_t * (1 / done_part - 1) : -1;
	int i, b_i = done_part > 0 ? int(left_t) : -1;
	time_t curr_ctime = int(gt + curr_t);
	time_t done_ctime = int(gt + curr_t + left_t);

	cout << proc_name << "\n"
		 << 100*done_part << " %          \n"
		 << "time spent " << int(real_t)/3600 << ":" << (int(real_t)%3600)/60 << ":" << int(real_t)%60 << "          \n"
		 << "time left " << b_i/3600 << ":" << (b_i%3600)/60 << ":" << b_i%60 << "          \n"
		 << "last save: " << string(ctime(&curr_ctime))
		 << "finish   : " << string(ctime(&done_ctime))
		 << "\r";      // goto begin & erase
	for(i = 0; i < extra_strN + 6; ++i){
		cout << "\033[A"; // up & erase
	}
	/*
		 * At first sight it seems nothing will be printed because I print and erase all right away.
		 * But I suppose that in some case \smth affect output-buffer but doesn't cause any screen.reprint()
		 * so everything (empty strings from previous cout) is printed in the beginning of next time print when actual chars are printed
		 * So in fact user sees everything erased just before new information to be printed
		 */
}

// -------------------------- string++ ---------------------------------
vector<double> d3ToV(double3 v){
	double vp[3] = {v.x, v.y, v.z};
	vector<double> vv;
	vv.assign(vp, vp+3);
	return vv;
}

string d3ToStr(double3 d)
{
	return "(" + toString(d.x) + ";" + toString(d.y) + ";" + toString(d.z) + ")";
}

template<typename T>
string vectorToStr(vector<T> v, string sp)
{
        string s="";
        for(int i = 0; i < v.size()-1; ++i) s+= (toString(v[i])+sp);
        return s+toString(v[v.size()-1]);
}
template<typename T>
void printVector(vector<T> v,string sp1, string sp2, string sp3)
{
        cout << sp1;
        for(int i = 0; i < v.size()-1; ++i) cout << v[i] << sp2;
        cout << v[v.size()-1] << sp3;
}
template<typename T>
T sumVector(vector<T> v)
{
        T s=0;
        for(int i = 0; i < v.size(); ++i) s+=v[i];
        return s;
}

// -------------------------- errors handling ---------------------------------
template<typename Ttype>
int SayLog(string logFname, Ttype s, bool printOnScreen){
	return sayError(logFname, SayIt, s, printOnScreen);
}

string getErrStr(int n, string s)
{
	string s_ret;
    switch(n){
		case SayIt: s_ret = ""; break;

		case CantOpenFile: s_ret = "Can't open file\n"; break;
		case CantOpenPrevFile: s_ret = "Can't open .prev file\n"; break;
		case CantCreateFile: s_ret = "Can't create file\n"; break;
		case WrongScriptFormat: s_ret = ""; break;
		case Nis1: s_ret = "N == 1"; break;
		case NLessOrEq0: s_ret = "N <= 0"; break;
		case NisTooBig: s_ret = "N is too big: maxN = " + toString(MaxN); break;
		case ALessOrEq0: s_ret = "a (base R) <= 0"; break;
		case WrongVelocities: s_ret = "Wrong Velocity mod"; break;
		case DtIs0: s_ret = "dt == 0"; break;
		case dtIsInf: s_ret = "dt is inf"; break;
		case DtIsLess0: s_ret = "dt < 0"; break;
		case WrongCompMode: s_ret = "Wrong computation mod"; break;
		case krLessOrEq0: s_ret = "kr <= 0"; break;
		case khLessOrEq0: s_ret = "kh <= 0"; break;
		case kaLess0: s_ret = "ka < 0"; break;
		case kvLess0: s_ret = "kv < 0"; break;
		case nCr_Less0: s_ret = "n_cr < 0"; break;
		case QLess0: s_ret = "Q < 0"; break;
		case StopFileAlreadyExists: s_ret = "Stop file already existed then computation has started"; break;
		case ConstForCommonMagic: s_ret = "Nothing to worry about, keep doing your work, everything is OK! Sorry for the strange bug! :-)"; break;
		case WrongPostProcMod: s_ret = "Wrong post processing mod, you entered"; break;
		case FnumLess0: s_ret = "Fnum (number of current out galaxy file) < 0"; break;
		case CalcedT0Less0: s_ret = "totalT0 (model time which has computed already) < 0"; break;
		case MLessOrEq0: s_ret = "M <= 0"; break;
		case TLessOrEq0: s_ret = "endT <= 0"; break;
		case TlessCalcedT: s_ret = "endT < totalT"; break;
		case dumpDTLessOrEq0: s_ret = "dumpDT <= 0"; break;
		case Geq0: s_ret = "G == 0"; break;
		case QIsInf: s_ret = "Q is inf"; break;
		case AIsInf: s_ret = "a is inf"; break;
		case krIsInf: s_ret = "kr is inf"; break;
		case khIsInf: s_ret = "kh is inf"; break;
		case kaIsInf: s_ret = "ka is inf"; break;
		case kvIsInf: s_ret = "kv is inf"; break;
		case MIsInf: s_ret = "M is inf"; break;
		case HIsInf: s_ret = "H is inf"; break;
		case TIsInf: s_ret = "endT is inf"; break;
		case EpIsInf: s_ret = "Ep is inf"; break;
		case EkIsInf: s_ret = "Ek is inf"; break;
		case LIsInf: s_ret = "L is inf"; break;
		case nCrIsInf: s_ret = "n_cr is inf"; break;
		case dumpDTIsInf: s_ret = "dumpDT is inf"; break;
		case GIsInf: s_ret = "G is inf"; break;
		case GisNot1: s_ret = "G != 1 (for efficiency)"; break;
		case StarXIsInf: s_ret = "X of some particle is inf"; break;
		case StarMIsInf: s_ret = "M of some particle is inf"; break;
		case StarRIsInf: s_ret = "R of some particle is inf"; break;
		case StarVIsInf: s_ret = "V of some particle is inf"; break;
		case NULLmalloc: s_ret = "Bad host memory allocation"; break;
		case WrongChngParamInp: s_ret = "chng_params : No parameter with name"; break;
		case WrongCondProcInp: s_ret = "Format for condition post_proc is:\n ./post_proc    mode(=3)     name_of_system    dr\n"; break;
		case dissipKless0: s_ret = "dissipK < 0"; break;
		case RcutLess0: s_ret = "r_cut < 0"; break;
		case TmpLess0: s_ret = "initial temperature < 0"; break;
		case muLessEq0: s_ret = "mu <= 0"; break;
		case TmpStabEps_LessEq0: s_ret = "TmpStabEps <= 0"; break;
		case TmpStabGap_Less0: s_ret = "TmpStabGap < 0"; break;
		case TooDenseSystem: s_ret = "N/V > 1 - not allowed conditions\n"; break;
		case NoValidTforP: s_ret = "Attempt to compurte Pressure with no valid T"; break;
		case AttemptToPrintOutsyncedData: s_ret = "Attempt to save out-synced data"; break;
		case CantCenterCM: s_ret = "Can't center CM"; break;
		case NisntCubeForCristal: s_ret = "System density may result cristal state but N isn't valid for generating one. N must be N * 2 = m^3 for FCC lattice."; break;
		case YetUnsupportedInput: s_ret = "Input params you set aren't fully supported currently"; break;
		case NthreadsLessOrEq0: s_ret = "Nthreads <= 0"; break;

		case CUDA_ERROR_CHECK_CONST: s_ret = ""; break;
		case CUDA_WRONG_DevMem_ALLOC: s_ret = "cuda: wrong device memory allocation"; break;
		case CUDA_WRONG_COPY_TO_DevMem: s_ret = "cuda: wrong copy to device memory"; break;
		case CUDA_WRONG_COPY_TO_HostMem: s_ret = "cuda: wrong copy to host memory"; break;
		case CUDA_WRONG_NBlockW: s_ret = "cuda: wrong N: N%" + toString(BlockW) + " != 0"; break;
		case cudaErrorUnknown: s_ret = "cudaErrorUnknown:"; break;

        default: s_ret = "\nUnknown error #" + toString(n);
    }
    //s_ret = s_ret.empty() ? s : (s_ret + "\n(" + s + ")\n");

    return s_ret.empty() ? s : (s_ret + "\n(" + s + ")\n");
}

template<typename Ttype>
int sayError(string logFname, int n, Ttype s, bool printOnScreen)
{
	//return 0;
    if(!n){ return 0; }
    ofstream Fout(logFname.c_str(),ios::app);
    if(!Fout){ return CantCreateFile; }

    string _s = "error #" + toString(n) + "\n";
    if(n != SayIt){
    	Fout << "\n" << _s;
    	if(printOnScreen) cerr << _s;
    }
    _s = getErrStr(n, toString(s));
	Fout << _s;
	if(printOnScreen) cerr << _s;
    if(n != SayIt){ Fout << "\n"; }

    Fout.close();

    return n;
}

// CUDA_CHECK(code_instruction) - for use
void MYcudaCheck(cudaError_t code, const char *file, int line, bool abort)
{
	if (code != cudaSuccess)
	{
		sayError("CUDA_CHECK_errors.log", CUDA_ERROR_CHECK_CONST,
				"\ncheckCudaErrors:\n" +
      			string("error: ") + string(cudaGetErrorString(code)) + "\n" +
      			"file = \"" + string(file) + "\"\n" +
      			"line = " + toString(line) + "\n\n");
		if(abort) exit(code);
	}
}

void DoCommonMagic(void)
{
	string s = "str";
    fromString<int>("1");
    fromString<int>((char*)" ");
    toString(1);
    toString((int)1);
    toString((bool)1);
    sayError("aaa", 1,"abc");
    sayError(s, 1,"abc");
    stp("start");
}

template <typename Ttype> string toString(Ttype val)
{
    std::ostringstream oss;
    oss<< val;
    return oss.str();
}

template<typename Ttype> Ttype fromString(const string& s)
{
  std::istringstream iss(s);
  Ttype res;
  iss >> res;
  return res;
}

string toLower(string s)
{
	char d = 'a'-'A';
	for(int i = 0; i<s.size(); ++i) if((s[i]>='A') && (s[i]<='Z')) s[i]+=d;
	return s;
}
string toUpper(string s)
{
	char d = 'A'-'a';
	for(int i = 0; i<s.size(); ++i) if((s[i]>='a') && (s[i]<='z')) s[i]+=d;
	return s;
}

template<typename Ttype> int addToLogFile(Ttype s)
{
    ofstream Fout("log.log",ios::app);
    if(!Fout){ return CantCreateFile; }

    Fout << s;

    Fout.close();
    return 0;
}
