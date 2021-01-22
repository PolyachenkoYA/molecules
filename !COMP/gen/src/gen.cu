#include "../../Space.cuh"

int main(int argc, char* argv[])
{
	TSpace Space;
	if((argc == 2) || (argc == 3)){
		//int er;
		Space.GlxName = string(argv[1]);
		Space.NewGlxName = (argc == 3 ? string(argv[2]) : Space.GlxName);
		Space.LogFName = Space.NewGlxName + "_gen.log";
		Space.sayLog("generation session " + Space.sessionID);
		string b_s;
		if(sayError(Space.LogFName, Space.SafeLoadParams(&b_s, Space.GlxName), "LoadParams")) return 1;

		Space.thermostatOn = (Space.dissipK > 0) && (Space.TmpStabGap > 0);
		Space.dissipK = std::abs(Space.dissipK);

		if(sayError(Space.LogFName, Space.CreateParticles(),"SetupStars")) return 1;

		if(Space.thermostatOn){
			if(sayError(Space.LogFName, Space.stabilizeTmp(),"stabilizeTmp")) return 1;
		} else {
			Space.computeEk();
			Space.Tmp = Space.Tmp_curr;
		}

        b_s = "./" + Space.NewGlxName + "_" + Space.SGname;
		if(sayError(Space.LogFName, Space.SaveParticles(b_s,Space.binOutF,PostProcLoadMode),b_s)) return 1;
		cout << "New system was saved in '" << b_s << "'\n";

		b_s = "./" + (argc == 3 ? string(argv[2]) : Space.NewGlxName) + "_" + Space.PrmFName;
		if(sayError(Space.LogFName, Space.SaveParams(b_s), b_s)) return 1;
		cout << "Parameters were resaved in '" << b_s << "'\n";

		Space.sayLog("all is ok | gen " + Space.sessionID);
	} else {
    	cout << "Format is:\n"
    		 <<	argv[0] << "   model_name   [new_model_name]\n";
    	print_sys_info(cout);
    	return 1;
	}
    return 0;
}

