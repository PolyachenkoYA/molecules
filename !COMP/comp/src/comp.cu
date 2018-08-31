/*
 ============================================================================
 Name        : calc.cu
 Author      : yura
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
        Space.GlxName = argv[1];
        Space.LogFName = Space.GlxName + "/comp.log";

        mkdir(string("./"+Space.GlxName).c_str(), S_IRWXU | S_IRWXG);                                 // whole model folder
        mkdir(string("./"+Space.GlxName + "/" + Space.FramesFilesFolder).c_str(), S_IRWXU | S_IRWXG); // frames folder

        Space.sayLog("computation session " + Space.sessionID);

		string b_s;
		if(sayError(Space.LogFName, Space.SafeLoadParams(&b_s), "LoadParams")) return 1;
		Space.thermostatOn = (Space.dissipK < 0);
		Space.dissipK = std::abs(Space.dissipK);

        b_s = "./" + Space.GlxName + "_" + Space.SGname;
        if(sayError(Space.LogFName, Space.LoadParticles(b_s, Space.binOutF, CompLoadMode),b_s)) return 2;

        if(sayError(Space.LogFName, Space.main_compute(),"main_compute")) return 3;

        b_s = "./"+Space.GlxName+"/"+Space.PrmFName;
        if(sayError(Space.LogFName, Space.SaveParams(b_s), b_s)) return 4;
        cout << Space.PrmFName << " was copied to '" << b_s << "'\n";

        b_s = "./"+Space.GlxName+"/"+Space.HeadFName;
        if(sayError(Space.LogFName, Space.CreateHeadF(b_s), b_s)) return 5;
        cout << Space.HeadFName << " was created in '" << b_s << "'\n";
    } else{
    	cout << "Usage:\n"
    		 << "./script_name galaxy_name\n";
    	return 1;
    }

    Space.sayLog("all is ok | calc " + Space.sessionID);
    return 0;
}
