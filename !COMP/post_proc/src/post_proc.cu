/*
 ============================================================================
 Name        : post_proc.cu
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

    if(argc >= 3){
        // define vars
        int i, j;
        int spPerNum = 20;
        ofstream Fout;
        string buf;
        time_t rStime;

        Space.GlxName = argv[1];
        Space.LogFName = "./" + Space.GlxName + "/post_proc.log";

        int procMod = fromString<int>(argv[2])+PostProcModID;

        // check proc mode
        if((procMod != EPostProcMod) &&
           (procMod != FromToBinPostProcMod) &&
           (procMod != ConditionPostProcMod)){
        		return sayError(Space.LogFName, WrongPostProcMod, toString(procMod)); // toString(procMod) FOR SERVER COMPILATION;
        }
        if(procMod == FromToBinPostProcMod){
        	sayError(Space.LogFName, SayIt, "bin <-> txt isn't supported yet\n supported modes:\n E(1)\n conditionFnc(3)\n ");
        	//cout << "bin <-> txt isn't supported yet. See post_proc.log for more\n";
        	return WrongPostProcMod;
        }
        double dr = -1;
        if(procMod == ConditionPostProcMod){
        	if(argc < 4){
            	sayError(Space.LogFName, WrongCondProcInp, "0"); // "0"  FOR SERVER COMPILATION
            	return WrongCondProcInp;
        	}
        	Space.Nslice = fromString<int>(argv[3]);
        }

        // create log file & check for stop file
        Space.sayLog("post processing session " + Space.sessionID + "   procMod = " + toString(procMod-PostProcModID) + "\n");
    	if(!access(Space.StpFName.c_str(),0)){
    		Space.sayLog(Space.StpFName + " file exists\n");
    		return 0;
    	}

    	// read head file if necessary
       	buf = "./"+Space.GlxName+"/"+Space.HeadFName;
       	ifstream Finh(buf.c_str());
       	if(!Finh){
       		sayError(Space.LogFName, CantOpenFile, buf);
       		return CantOpenFile;
       	}
       	Finh >> Space.Nfrm;
       	Finh.close();

       	// read params
       	buf = "./" + Space.GlxName + "/" + Space.PrmFName;
       	if(sayError(Space.LogFName, Space.LoadParams(buf),buf))	return 1;
    	dr = Space.R/Space.Nslice;

       	// create result-file
        Space.sayLog("   loading/create head/single file{\n");
        string outFname;
        switch(procMod){
            case EPostProcMod:
            	outFname = "./" + Space.GlxName + "/" + Space.EFile;
                break;
            case ConditionPostProcMod:
            	outFname = "./" + Space.GlxName + "/" + Space.ConditionFNameBase + "_" + toString(Space.Nslice) + ".dat";
            	break;
        }
        Fout.open(outFname.c_str());
        if(!Fout){
            sayError(Space.LogFName, CantCreateFile, outFname);
            return CantCreateFile;
        }
        Space.sayLog("   }loaded/created (" + outFname + ")\n");
        Fout.precision(PrintPres);

        string currentFileName, job_name;

        switch(procMod){
        	case EPostProcMod:
        		Space.sayLog("   E computation\n");
            	job_name = "post processing E";
        		break;
        	case FromToBinPostProcMod:
        		Space.sayLog("   changing binarity from " + toString(Space.binOutF) + // REMOVE FOR SERVER COMPILATION
        					 " to " + toString(!Space.binOutF) + "\n");               // REMOVE FOR SERVER COMPILATION
        		//Space.sayLog("   changing binarity from " + string(Space.binOutF ? "1" : "0") + // REMOVE FOR SERVER COMPILATION
        		//			 " to " + string(!Space.binOutF ? "1" : "0") + "\n");               // REMOVE FOR SERVER COMPILATION
            	job_name = "post processing changing binarity to " + toString(!Space.binOutF);
        		break;
        	case ConditionPostProcMod:
        		Space.sayLog("   Condition computation\n");
            	for(j = 0; j < Space.Nslice; ++j){ // print X for plot
            		Fout << setw(spPerNum) << (j + 0.5)*dr;
            	}
            	Fout << "\n";

            	Space.gCondFnc = (double*)malloc(sizeof(double)*Space.Nslice);
            	if(Space.gCondFnc == NULL){
            		sayError(Space.LogFName,NULLmalloc,"gCondFnc");
            		return NULLmalloc;
            	}
            	job_name = "post processing Condition";
        		break;
        }
        Space.sayLog("   post processing{\n");

        time(&rStime);
        for(i = 0; i < Space.Nfrm; ++i){
        	currentFileName = "./" + Space.GlxName + "/" + Space.FramesFilesFolder + "/" + toString(i) + "." + Space.ParticleFileExt;
            if(sayError(Space.LogFName, Space.LoadParticles(currentFileName,Space.binOutF,PostProcLoadMode),"LoadStars") > 0) return 3;
            switch(procMod){
                case EPostProcMod:
                	Space.sayLog("   E_" + toString(i) + " computation ...\n");
                	if(sayError(Space.LogFName, Space.useCuda ? Space.cudaDoEstep() : Space.doEstep(), "doEstep") > 0) return 5;
                	Space.sayLog("   ...computed\n");
                    Fout << setw(spPerNum) << Space.E.x << " "
                    	 << setw(spPerNum) << Space.E.y << " "
                    	 << setw(spPerNum) << Space.E.x + Space.E.y << endl;
                    break;
                case FromToBinPostProcMod:
                    if(sayError(Space.LogFName, Space.SaveParticles(currentFileName,!Space.binOutF,PostProcLoadMode),"SaveStars") > 0) return 6;
                    break;
                case ConditionPostProcMod:
                	Space.sayLog("   Condition_" + toString(i) + " computation ...\n");
                	if(sayError(Space.LogFName, Space.doCondStep(dr), "doCondStep") > 0) return 11;
                	for(j = 0; j < Space.Nslice; ++j){ // print Y for plot
                		//Fout << setw(spPerNum) << Space.gCondFnc[j];
                		Fout << Space.gCondFnc[j] << " ";
                	}
                	Fout << "\n";
                	Space.sayLog("   ...computed\n");
                	break;
            }
            time_progress(rStime, time(0), double(i+1) / Space.Nfrm, job_name);
        }
        Space.sayLog("   }post processing done\n");

        switch(procMod){
        	case EPostProcMod:
        		cout << "Post processing E is done (" << outFname << ")\n";
        		break;
        	case FromToBinPostProcMod:
        		cout << "Post processing changing binary to " << !Space.binOutF << " is done\n";
        		Space.binOutF = !Space.binOutF;
        		if(sayError(Space.LogFName, Space.SaveParams("./"+Space.GlxName+"/"+Space.PrmFName),"SaveParams") > 0) return 8;
        		if(sayError(Space.LogFName, Space.CreateHeadF("./"+Space.GlxName+"/"+Space.HeadFName),"CreateHeadF") > 0) return 9;
        		cout << "'" << Space.HeadFName << "' and '" << Space.PrmFName << "' were resaved in '" << Space.GlxName << "'\n";
        		break;
        	case ConditionPostProcMod:
        		cout << "Post processing Condition is done (" << outFname << ")\n";
        		break;
        }
        Fout.close();
    } else {
        cout << "Format is:\n"
        	 << "./script_name    galaxy_name    processing_mod    [dr(3)]\n";
        return 100;
    }

    Space.sayLog("all is ok |post_proc " + Space.sessionID);
    return 0;
}
