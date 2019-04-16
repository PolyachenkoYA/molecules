/*
 ============================================================================
 Name        : gen.cu
 Author      : PolyachYA
 Version     :
 Copyright   : Your copyright notice
 Description : CUDA compute reciprocals
 ============================================================================
 */

#include "../../Space.cuh"

int main(int argc, char* argv[])
{
	TSpace Space;
	if(argc == 2){
		//int er;
		Space.GlxName = string(argv[1]);
		Space.LogFName = Space.GlxName + "_gen.log";
		Space.sayLog("generation session " + Space.sessionID);
		string b_s;

		if(sayError(Space.LogFName, Space.SafeLoadParams(&b_s),"LoadParams")) return 1;
		Space.thermostatOn = (Space.dissipK > 0) && (Space.TmpStabGap > 0);
		Space.dissipK = std::abs(Space.dissipK);

		if(sayError(Space.LogFName, Space.CreateParticles(),"SetupStars")) return 1;

		if(Space.thermostatOn){
			if(sayError(Space.LogFName, Space.stabilizeTmp(),"stabilizeTmp")) return 1;
		}

        b_s = "./"+Space.GlxName+"_"+Space.SGname;
		if(sayError(Space.LogFName, Space.SaveParticles(b_s,Space.binOutF,PostProcLoadMode),b_s)) return 1;
		cout << "New system was saved in '" << b_s << "'\n";
		//b_s = "./"+Space.GlxName+"_new_"+Space.PrmFName;
		b_s = "./"+Space.GlxName+"_"+Space.PrmFName;
		if(sayError(Space.LogFName, Space.SaveParams(b_s),b_s)) return 1;
		cout << "Parameters were resaved in '" << b_s << "'\n";

		Space.sayLog("all is ok | gen " + Space.sessionID);
	} else {
    	cout << "Format is:\n"
    		 <<	argv[0] << " model_name\n";
    	print_sys_info(cout);
    	return 1;
	}
    return 0;
}
