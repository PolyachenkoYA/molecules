#include "format.cuh"

template<typename T> void stp(T str)
{
	cout << str;
	cin.get();
}

void time_progress(time_t real_start_t, time_t curr_t, double done_part, string proc_name)
{
	time_t real_t = curr_t - real_start_t;
	double left_t = real_t * (1/done_part - 1);
	int b_i = int(left_t);
	time_t b_t = curr_t + b_i;

	cout << proc_name << " " << 100*done_part << " %          \n"
		 << "time used " << real_t/3600 << ":" << (real_t%3600)/60 << ":" << real_t%60 << "          \n"
		 << "time left " << b_i/3600 << ":" << (b_i%3600)/60 << ":" << b_i%60 << "          \n"
		 << "last save: " << string(ctime(&curr_t))
		 << "finish   : " << string(ctime(&b_t))
		 << "\r"      // goto begin & erase
		 << "\033[A"  // up & erase
		 << "\033[A"  // up & erase
		 << "\033[A"  // up & erase
		 << "\033[A"  // up & erase
		 << "\033[A"; // up & erase
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
    switch(n){
    	case SayIt: return s;

        case CantOpenFile: return "Can't open file\n" + s + "\n";
        case CantOpenPrevFile: return "Can't open .prev file\n" + s + "\n";
        case CantCreateFile: return "Can't create file\n" + s + "\n";
        case WrongScriptFormat: return s + "\n";
        case Nis1: return "N == 1 (" + s + ")\n";
        case NLessOrEq0: return "N <= 0 (" + s + ")\n";
        case NisTooBig: return "N is too big: maxN = " + toString(MaxN) + " (" + s + ")\n";
        case ALessOrEq0: return "a (base R) <= 0(" + s + ")\n";
        case WrongVelocities: return "Wrong Velocity mod (" + s + ")\n";
        case DtIs0: return "dt == 0 (" + s + ")\n";
        case dtIsInf: return "dt is inf (" + s + ")\n";
        case DtIsLess0: return "dt < 0 (" + s + ")\n";
        case WrongCompMode: return "Wrong computation mod (" + s + ")\n";
        case krLessOrEq0: return "kr <= 0 (" + s + ")\n";
        case khLessOrEq0: return "kh <= 0 (" + s + ")\n";
        case kaLess0: return "ka < 0 (" + s + ")\n";
        case kvLess0: return "kv < 0 (" + s + ")\n";
        case QLess0: return "Q < 0 (" + s + ")\n";
        case StopFileAlreadyExists: return "Stop file already existed then computation has started (" + s + ")\n";
        case ConstForCommonMagic: return "Nothing to worry about, keep doing your work, everything is OK! Sorry for the strange bug! :-) (" + s + ")\n";
        case WrongPostProcMod: return "Wrong post processing mod, you entered " + s + " mode\n";
        case FnumLess0: return "Fnum (number of current out galaxy file) < 0 (" + s + ")\n";
        case CalcedT0Less0: return "totalT0 (model time which has computed already) < 0 (" + s + ")\n";
        case MLessOrEq0: return "M <= 0 (" + s + ")\n";
        case TLessOrEq0: return "endT <= 0 (" + s + ")\n";
        case TlessCalcedT: return "endT < totalT (" + s + ")\n";
        case dumpDTLessOrEq0: return "dumpDT <= 0 (" + s + ")\n";
        case Geq0: return "G == 0 (" + s + ")\n";
        case QIsInf: return "Q is inf (" + s + ")\n";
        case AIsInf: return "a is inf (" + s + ")\n";
        case krIsInf: return "kr is inf (" + s + ")\n";
        case khIsInf: return "kh is inf (" + s + ")\n";
        case kaIsInf: return "ka is inf (" + s + ")\n";
        case kvIsInf: return "kv is inf (" + s + ")\n";
        case MIsInf: return "M is inf (" + s + ")\n";
        case HIsInf: return "H is inf (" + s + ")\n";
        case TIsInf: return "endT is inf (" + s + ")\n";
        case dumpDTIsInf: return "dumpDT is inf (" + s + ")\n";
        case GIsInf: return "G is inf (" + s + ")\n";
        case GisNot1: return "G != 1 (for efficiency) (" + s + ")\n";
        case StarXIsInf: return "X of some star is inf (" + s + ")\n";
        case StarMIsInf: return "M of some star is inf (" + s + ")\n";
        case StarRIsInf: return "R of some star is inf (" + s + ")\n";
        case StarVIsInf: return "V of some star is inf (" + s + ")\n";
        case NULLmalloc: return "Bad host memory allocation (" + s + ")\n";
        case WrongChngParamInp: return "No parameter with name " + s + " (chng_params)\n";
        case WrongCondProcInp: return "Format for condition post_proc is:\n ./post_proc    mode(=3)     name_of_system    dr\n";
        case dissipKless0: return "dissipK < 0 (" + s + ")\n";
        case RbrdLess0: return "r_brd < 0 (" + s + ")\n";
        case TmpLess0: return "initial temperature < 0 (" + s + ")\n";
        case muLessEq0: return "mu <= 0 (" + s + ")\n";
        case TmpStabEps_LessEq0: return "TmpStabEps <= 0 (" + s + ")\n";
        case TmpStabGap_Less0: return "TmpStabGap < 0 (" + s + ")\n";
        case TooDenseSystem: return "N/V > 1 - not allowed conditions\n";

        case CUDA_ERROR_CHECK_CONST: return s;
        case CUDA_WRONG_DevMem_ALLOC: return "cuda: wrong device memory allocation (" + s + ")\n";
        case CUDA_WRONG_COPY_TO_DevMem: return "cuda: wrong copy to device memory (" + s + ")\n";
        case CUDA_WRONG_COPY_TO_HostMem: return "cuda: wrong copy to host memory (" + s + ")\n";
        case CUDA_WRONG_NBlockW: return "cuda: wrong N: N%" + toString(BlockW) + " != 0 (" + s + ")\n";
        case cudaErrorUnknown: return "cudaErrorUnknown:" + s;

        default: return s + "\nUnknown error #" + toString(n) + "\n";
    }
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
