#pragma once

#include "const.cuh"

template<typename T> void stp(T str); // hand-made breakpoint

// ----------------------------------------- usefull print stuff ------------------------------------------------
template<typename T> void printVector(vector<T> v, string sp1="\n", string sp2=" ", string sp3="\n");
template<typename T> string vectorToStr(vector<T> v, string sp=" ");
template<typename T> T sumVector(vector<T> v);
vector<double> d3ToV(double3 v);
string d3ToStr(double3 d);
template<typename T> string toString(T val);
template<typename T> T fromString(const string& s);
string toLower(string s0);
string toUpper(string s0);

void time_progress(time_t real_start_t, time_t curr_t, double done_part, string proc_name, int extra_strN = 0);

// ------------------------------------------ error handling --------------------------------------------------

string getErrStr(int n, string s);
template<typename Ttype> int sayError(string logFname, int n, Ttype s, bool printOnScreen = 0);
template<typename Ttype> int SayLog(string logFname, Ttype s, bool printOnScreen = 0);
template<typename Ttype> int addToLogFile(Ttype s);
void DoCommonMagic(void);

#define CUDA_CHECK(ans) { MYcudaCheck((ans), __FILE__, __LINE__); }
void MYcudaCheck(cudaError_t code, const char *file, int line, bool abort=1);
