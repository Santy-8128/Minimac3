

#include <iostream>
#include <ctime>
#include <stdio.h>
#include "Parameters.h"
#include "StringBasics.h"
#include "HaplotypeSet.h"
#include "Unique.h"
#include "MarkovModel.h"
#include "MarkovParameters.h"
#include "Imputation.h"
#include "ImputationStatistics.h"

using namespace std;
void Minimac3Version();
void helpFile();

int main(int argc, char ** argv)
{
	// Parameter Options

	String refHaps = "";
	String haps = "", snps = "",removeSam="";
	String outfile = "Minimac3.Output";
	String format = "GT,DS";
	String recFile = "", errFile = "",chr="",golden="";
	int cpus=1,start=0, end=0, window=0, max_indiv = 0, max_marker = 0, rounds=5, states=200;
    vector<bool> formatVector(3,false);
    #ifdef _OPENMP
    cpus=5;
    #endif

	bool phased = false,passOnly = false, doseOutput = false, vcfOutput = true, gzip = true, nobgzip = false, rsid=false;
	bool processReference=false,updateModel=false, onlyRefMarkers=false, help = false, params = false;

	ParameterList inputParameters;
	PhoneHome::allThinning = 50;

	BEGIN_LONG_PARAMETERS(longParameterList)
		LONG_PARAMETER_GROUP("Reference Haplotypes")
		LONG_STRINGPARAMETER("refHaps", &refHaps)
		LONG_PARAMETER("passOnly", &passOnly)
		//LONG_PARAMETER("rsid", &rsid)
		LONG_PARAMETER_GROUP("Target Haplotypes")
		LONG_STRINGPARAMETER("haps", &haps)
//		LONG_STRINGPARAMETER("snps", &snps)
		LONG_PARAMETER_GROUP("Output Parameters")
		LONG_PARAMETER("processReference", &processReference)
		LONG_STRINGPARAMETER("prefix", &outfile)
		LONG_PARAMETER("updateModel", &updateModel)
		LONG_PARAMETER("nobgzip", &nobgzip)
//		LONG_PARAMETER("vcfOutput", &vcfOutput)
		LONG_PARAMETER("doseOutput", &doseOutput)
		LONG_PARAMETER("hapOutput", &phased)
		LONG_STRINGPARAMETER("format", &format)
		LONG_PARAMETER_GROUP("Subset Parameters")
		LONG_STRINGPARAMETER("chr", &chr)
		LONG_INTPARAMETER("start", &start)
		LONG_INTPARAMETER("end", &end)
		LONG_INTPARAMETER("window", &window)
		//LONG_INTPARAMETER("block", &max_block)
		LONG_PARAMETER_GROUP("Starting Parameters")
		LONG_STRINGPARAMETER("rec", &recFile)
		LONG_STRINGPARAMETER("err", &errFile)
		LONG_PARAMETER_GROUP("Estimation Parameters")
		LONG_INTPARAMETER("rounds", &rounds)
		LONG_INTPARAMETER("states", &states)
		LONG_PARAMETER_GROUP("Other Parameters")
		LONG_PARAMETER("help", &help)
		LONG_INTPARAMETER("cpus", &cpus)
		LONG_PARAMETER("params", &params)
		LONG_PHONEHOME(VERSION)
		BEGIN_LEGACY_PARAMETERS()
		LONG_PARAMETER("onlyRefMarkers", &onlyRefMarkers)
		LONG_STRINGPARAMETER("golden", &golden)
		LONG_INTPARAMETER("sample", &max_indiv)
		LONG_INTPARAMETER("marker", &max_marker)
		LONG_STRINGPARAMETER("remove", &removeSam)
		END_LONG_PARAMETERS();


	inputParameters.Add(new LongParameters(" Command Line Options: ",longParameterList));
//    PhoneHome::allThinning = 100;
	Minimac3Version();
    String compStatus;
	inputParameters.Read(argc, &(argv[0]));
	if (help)
	{
		helpFile();
		return(-1);
	}

	inputParameters.Status();

    #ifdef _OPENMP
    omp_set_num_threads(cpus);
    #endif

    if(nobgzip)
        gzip=false;


    if (rsid==true)
    {
        cout << " NOTE : If \"--refHaps\" is an OPTM file, \"--rsid\" parameter provided will be ignored !!! \n";
        cout << "        Program will read marker information from OPTM file : " <<refHaps<<endl<<endl;
    }

	if (refHaps == "")
	{
		cout<< " Missing \"--refHaps\", a required parameter.\n";
		cout<< " Type \"--help\" for more help.\n\n";
		compStatus="Command.Line.Error";
		PhoneHome::completionStatus(compStatus.c_str());
		return(-1);
	}


    if(processReference)
    {

        cout<<" NOTE: Since \"--processReference\" is ON, all options under \"Target Haplotypes\" \n";
        cout<<"       and \"Starting Parameters\" will be ignored !!!\n";
        cout<<"       Program will only estimate parameters and create OPTM file.\n";
        cout<<"       No imputation will be performed, hence other parameters are unnecessary !!!"<<endl<<endl;

        cout<<" NOTE: If \"--processReference\" is ON, Parameter Estimation will be done by default ! \n";
        cout<<"       Use \"--rounds 0\" to AVOID Parameter Estimation !!!\n"<<endl<<endl;

        if(updateModel)
        {

            cout<<" Handle \"--updateModel\" does NOT work with handle \"--processReference\" !!! \n";
            cout<<" Type \"--help\" for more help.\n\n";
            compStatus="Command.Line.Error";
            PhoneHome::completionStatus(compStatus.c_str());
            return(-1);
		}

    }

    if(updateModel)
    {

        cout<<" NOTE: Handle \"--updateModel\" works only on M3VCF files ! \n";
        cout<<"       Program will crash if \"--refHaps\" is a VCF file !!!\n"<<endl;

        if(rounds<=0)
        {
            cout << " Invalid input for \"--rounds\" = "<<rounds<<"\n";;
            cout << " Value must be POSITIVE if \"--updateModel\" is ON !!! \n\n";
            compStatus="Command.Line.Error";
            PhoneHome::completionStatus(compStatus.c_str());
            return(-1);
        }
        if(states<=0)
        {
            cout << " Invalid input for \"--states\" = "<<states<<"\n";;
            cout << " Value must be POSITIVE if \"--updateModel\" is ON !!! \n\n";
            compStatus="Command.Line.Error";
            PhoneHome::completionStatus(compStatus.c_str());
            return(-1);
        }
    }

    if(rounds<0)
    {
        cout << " Invalid input for \"--rounds\" = "<<rounds<<"\n";;
        cout << " Value must be non-negative !!! \n\n";
        compStatus="Command.Line.Error";
		PhoneHome::completionStatus(compStatus.c_str());
		return(-1);
    }
    if(states<0)
    {
        cout << " Invalid input for \"--states\" = "<<states<<"\n";;
        cout << " Value must be non-negative !!! \n\n";
        compStatus="Command.Line.Error";
		PhoneHome::completionStatus(compStatus.c_str());
		return(-1);
    }


    if(window<0)
    {
        cout << " Invalid input for \"--window\" = "<<window<<"\n";;
        cout << " Value must be non-negative !!! \n\n";
        compStatus="Command.Line.Error";
		PhoneHome::completionStatus(compStatus.c_str());
		return(-1);
    }

    if(start<0)
    {
        cout << " Invalid input for \"--start\" = "<<start<<"\n";;
        cout << " Value must be non-negative !!! \n\n";
        compStatus="Command.Line.Error";
		PhoneHome::completionStatus(compStatus.c_str());
		return(-1);
    }
    if(end<0)
    {
        cout << " Invalid input for \"--end\" = "<<end<<"\n";;
        cout << " Value must be non-negative !!! \n\n";
        compStatus="Command.Line.Error";
		PhoneHome::completionStatus(compStatus.c_str());
		return(-1);
    }
    if(start>0 && end> 0 && start>=end)
    {
        cout << " Invalid Input !!!\n Value of \"--start\" must be less than value of \"--end\"."<<endl;
        cout << " User Input \"--start\" = "<<start<<" and \"--end\" = " <<end<<" \n\n";
        compStatus="Command.Line.Error";
		PhoneHome::completionStatus(compStatus.c_str());
		return(-1);
    }
    if(start>0)
    {
        if(chr=="")
        {
            cout << "\n Missing \"--chr\", a required parameter if using \"--start\" parameter.\n";
            cout << " Try --help for more information.\n\n";
            compStatus="Command.Line.Error";
            PhoneHome::completionStatus(compStatus.c_str());
            return(-1);
        }
        if(end==0)
        {
            cout << "\n Non-zero value of \"--end\" required parameter if using \"--start\" parameter.\n";
            cout << " Try --help for more information.\n\n";
            compStatus="Command.Line.Error";
            PhoneHome::completionStatus(compStatus.c_str());
            return(-1);
        }
        if(window==0)
        {
            window=500000;
            cout<<" NOTE: Default \"--window\" parameter to be used = 500000\n";
        }
    }
    if(end>0)
    {
        if(chr=="")
        {
            cout << "\n Missing \"--chr\", a required parameter if using \"--end\" parameter.\n";
            cout << " Try --help for more information.\n\n";
            compStatus="Command.Line.Error";
            PhoneHome::completionStatus(compStatus.c_str());
            return(-1);
        }
        if(start==0)
        {
            cout << "\n Non-zero value of \"--start\" required parameter if using \"--end\" parameter.\n";
            cout << " Try --help for more information.\n\n";
            compStatus="Command.Line.Error";
            PhoneHome::completionStatus(compStatus.c_str());
            return(-1);
        }
        if(window==0)
        {
            window=500000;
            cout<<" NOTE: Default \"--window\" parameter to be used = 500000\n";
        }
    }
    if(window>0)
    {

        if(chr=="")
        {
            cout << "\n Missing \"--chr\", a required parameter if using \"--end\" parameter.\n";
            cout << " Try --help for more information.\n\n";
            compStatus="Command.Line.Error";
            PhoneHome::completionStatus(compStatus.c_str());
            return(-1);
        }
        if(start==0 && end==0)
        {
            cout << "\n Missing \"--start\" or  \"--end\", a required parameter if using \"--window\" parameter.\n";
            cout << " Try --help for more information.\n\n";
            compStatus="Command.Line.Error";
            PhoneHome::completionStatus(compStatus.c_str());
            return(-1);
        }
    }
    else
    {
        if(start>0 || end>0)
         {
            cout<<" NOTE: No \"--window\" parameter provided  !!! \n";
            cout<<"       No buffer region will be used on either side of the chunk"<<endl<<endl;
        }
    }


    if(!processReference)
    {

         if (haps == "")
        {
            cout <<" Missing \"--haps\", a required parameter (for imputation).\n";
            cout <<" OR use \"--processReference\" to just process the reference panel.\n";
            cout<< " Type \"--help\" for more help.\n\n";
            compStatus="Command.Line.Error";
            PhoneHome::completionStatus(compStatus.c_str());
            return(-1);
        }

        if (snps == "")
        {
            cout<<" NOTE : \"--snps\" parameter not provided ! Please verify that \"--haps\" is NOT a MaCH file. \n";
            cout<<"        Program will crash otherwise !!!"<<endl<<endl;
        }
        else
        {
            cout << " NOTE : If \"--haps\" is a VCF file, \"--snps\" parameter provided will be ignored !!! \n";
            cout << "        Program will read marker information from VCF file : " <<haps<<endl<<endl;
        }
//        if(recFile=="" && errFile!="")
//        {
//            cout << " Both \"--rec\" and \"--err\" parameters must be provided (or NONE) !!! \n";
//            cout << " Only \"--err\" is provided. \n\n";
//            compStatus="Command.Line.Error";
//            PhoneHome::completionStatus(compStatus.c_str());
//            return(-1);
//        }
//        if(errFile=="" && recFile!="")
//        {
//            cout << " Both \"--rec\" and \"--err\" parameters must be provided (or NONE) !!! \n";
//            cout << " Only \"--rec\" is provided. \n\n";
//            compStatus="Command.Line.Error";
//            PhoneHome::completionStatus(compStatus.c_str());
//            return(-1);
//        }

        typedef boost::tokenizer< boost::char_separator<char> > wsTokenizer;
        wsTokenizer::iterator i;
        string formatPiece,formatTemp=format.c_str();
        boost::char_separator<char> comSep(",");
        wsTokenizer t(formatTemp,comSep);
        for(i = t.begin();i!=t.end();++i)
        {

            formatPiece=i->c_str();
            if(formatPiece.compare("GT")==0)
                formatVector[0]=true;
            else if(formatPiece.compare("DS")==0)
                formatVector[1]=true;
            else if(formatPiece.compare("GP")==0)
                formatVector[2]=true;
            else
            {
                cout << " Cannot identify handle for \"--format\" parameter : "<<formatPiece<<endl;
                cout << " Available handles GT, DS and GP (for genotype, dosage and posterior probability). \n\n";
                cout << " Type \"--help\" for more help.\n\n";
                compStatus="Command.Line.Error";
                PhoneHome::completionStatus(compStatus.c_str());
                return(-1);
            }
        }
    }

    HaplotypeSet target,reference;

	if(!processReference)
    {
        cout<<" ------------------------------------------------------------------------------"<<endl;
        cout<<"                            PRELIMNARY FILE CHECK                              "<<endl;
        cout<<" ------------------------------------------------------------------------------"<<endl;


	    if (!target.BasicCheckForTargetHaplotypes(haps))
        {
            cout << "\n Program Exiting ... \n\n";
            compStatus="Target.Panel.Load.Error";
            PhoneHome::completionStatus(compStatus.c_str());
            return(-1);
        }

        std::cout << " Performing basic file check on Reference haplotype file ..." << endl;

        if(reference.DetectReferenceFileType(refHaps).compare("m3vcf")==0)
        {
            cout<<"\n Reference File Format = M3VCF (Minimac3 VCF File) "<<endl;

            cout <<"\n NOTE: For M3VCF files, if estimates are available in file \n";
               cout<<"       \"--updateModel\" must be ON for further parametrization.\n";
               cout<<"       Else it will use the estimates available in the file.\n";
               cout<<"       If estimates are NOT available in file, it will estimate by default."<<endl;


            if(rounds==0)
            {

                cout <<"\n NOTE: User has specified \"--rounds\" = 0 !!!\n";
                cout<<"       Please verify that the M3VCF file has parameter estimates in it.\n";
                cout<<"       Otherwise program will use default estimates leading to possibly inaccurate estimates."<<endl;


            }
            else
            {
                if(!updateModel)
                {
                 cout <<"\n NOTE: For M3VCF files, if estimates are available in file \n";
                    cout<<"       \"--updateModel\" must be ON to use \"--rounds\">0.\n";
                    cout<<"       Program will ignore value of \"--rounds\" otherwise !!!"<<endl;
                }
            }


        }
        if(reference.DetectReferenceFileType(refHaps).compare("vcf")==0)
        {
            cout<<"\n Reference File Format = VCF (Variant Call Format)"<<endl;

            cout <<"\n NOTE: For VCF files, parameter estimation will be done by default (unless \"--rounds\" = 0)."<< endl;

            if(updateModel && recFile=="" && errFile=="")
            {
                cout << "\n Handle \"--updateModel\" does NOT work for VCF reference file.\n";
                cout << " This works only for M3VCF files or when \"--rec\" or \"--err\" is provided.\n";
                cout << " For VCF files, parameter estimation will be done by default (unless \"--rounds\" = 0).\n";
                cout << " Please turn OFF \"--updateModel\" or use M3VCF file.\n";
                cout << " Try --help for more information.\n\n";
                compStatus="Command.Line.Error";
                PhoneHome::completionStatus(compStatus.c_str());
                return(-1);
            }

            if(rounds==0)
            {

                cout <<"\n NOTE: User has specified \"--rounds\" = 0 !!!\n";
                cout<<"       No parameter estimation will be done on VCF file.\n";
                cout<<"       Program will use default estimates leading to possibly inaccurate estimates."<<endl;
            }
        }
        cout<<endl;

    }

	int start_time = time(0);
	int time_prev = start_time;

    cout<<" ------------------------------------------------------------------------------"<<endl;
    cout<<"                           REFERENCE HAPLOTYPE PANEL                           "<<endl;
    cout<<" ------------------------------------------------------------------------------"<<endl;


	reference.removeSam=removeSam;

	if (!reference.FastLoadHaplotypes(refHaps, max_indiv, max_marker,chr,start,end,window,rsid,processReference,passOnly))
	{
		cout << "\n Program Exiting ... \n\n";
		compStatus="Reference.Panel.Load.Error";
        PhoneHome::completionStatus(compStatus.c_str());
        return(-1);

	}

//
//abort();

    cout<<endl;
	if(processReference)
    {
	    Imputation thisDataFast(reference, reference, outfile, errFile, recFile, phased
                             , gzip, rounds, states, vcfOutput, doseOutput
                             , onlyRefMarkers,formatVector);
        thisDataFast.createEstimates(reference, reference, reference.optEndPoints, true);
	}

	int time_load = time(0) - time_prev;
	time_prev = time(0);



    cout << "\n Time taken to load reference haplotype set = " << time_load << " seconds."<<endl<<endl;


//    MINIMIZED COST              = " << sum<< endl;
//    cout<<" FOLD INCREASE IN SPEEED     = "<< (double)rHap.numHaplotypes*(double)rHap.numMarkers /(double)sum<<endl;
//    cout<<" MAX_BLOCK_LENGTH            = "<< maxmax <<endl;


	if(!processReference)
    {
        cout<<" ------------------------------------------------------------------------------"<<endl;
        cout<<"                          TARGET/GWAS HAPLOTYPE PANEL                         "<<endl;
        cout<<" ------------------------------------------------------------------------------"<<endl;


	    if (!target.LoadTargetHaplotypes(haps, snps, reference.markerName,reference))
        {

            cout << "\n Program Exiting ... \n\n";
            compStatus="Target.Panel.Load.Error";
            PhoneHome::completionStatus(compStatus.c_str());
            return(-1);
        }
        time_load = time(0) - time_prev;
        time_prev = time(0);
        cout << "\n Time taken to load target haplotype set = " << time_load << " seconds. "<<endl<<endl;
    }




    if(!processReference)
	{
	    Imputation thisDataFast(target, reference, outfile, errFile, recFile, phased
                             , gzip, rounds, states, vcfOutput, doseOutput
                             , onlyRefMarkers,formatVector,updateModel);
        thisDataFast.performImputation(target, reference, golden);
	}


    cout<<" ------------------------------------------------------------------------------"<<endl;
    cout<<"                                END OF PROGRAM                                 "<<endl;
    cout<<" ------------------------------------------------------------------------------"<<endl;

    time_load = time(0) - time_prev;
	int time_tot = time(0) - start_time;

//	cout << "\n\n Time taken for Imputation = " << time_load << " seconds. \n";
//

    cout << "\n Program Successfully Implemented... \n ";


	printf("\n Total Run completed in %d hours, %d mins, %d seconds.\n",
		time_tot / 3600, (time_tot % 3600) / 60, time_tot % 60);

    cout<<"\n Thank You for using Minimac3 !!! "<<endl<<endl;



    compStatus="Success";
    PhoneHome::completionStatus(compStatus.c_str());

	return 0;

}




void Minimac3Version()
{
	printf("\n\n -------------------------------------------------------------------------------- \n");
	printf("          Minimac3 - Fast Imputation Based on State Space Reduction HMM\n");
	printf(" --------------------------------------------------------------------------------\n");
    printf(" (c) 2014 - Sayantan Das, Christian Fuchsberger, Goncalo Abecasis, Mary Kate Wing \n");
//	printf(" Version	: Undocumented Release\n");
//	printf(" Built		: sayantan\n\n");
	cout<<"\n Version: " << VERSION<< ";\n Built: " << DATE << " by " << USER << std::endl;
}

void helpFile()
{


	printf("\n URL = http://genome.sph.umich.edu/wiki/Minimac3\n");


    printf("\n\n\t Minimac3 is a lower memory and more computationally efficient implementation of \"minimac\".\n");

    printf("\t It is an algorithm for genotypic imputation that works on phased genotypes (say from MaCH).\n");
    printf("\t Minimac3 is designed to handle very large reference panels in a more computationally efficient\n");
    printf("\t way with no loss of accuracy.\n\n");

    printf("\t Please visit web-page <http://genome.sph.umich.edu/wiki/Minimac3> for further details on \n");
    printf("\t documentation and usage\n\n\n");

    printf(" The most commonly used parameters are explained below:\n\n");

    printf("   List of Usage Options  :\n\n");

    printf("     --refHaps filename   : VCF file or M3VCF  file containing haplotype data for reference panel.\n");
  printf("\n        --haps filename   : File containing haplotype data for target (gwas) samples. Must be VCF \n");
    printf("                            file. Zipped versions allowed.\n");
  printf("\n        --prefix output   : Prefix for all output files generated. By default: [Minimac3.Output]\n");
    printf("              --nobgzip   : If ON, output files will NOT be gzipped.\n");
    printf("     --processReference   : This option will only convert an input VCF file to M3VCF format\n");
    printf("                            (maybe for a later run of imputation). If this option is ON, \n");
    printf("                            no imputation would be performed and thus all target parameters \n");
    printf("                            will be ignored (of course, except for parameters on Reference \n");
    printf("                            Haplotypes and Subsetting Options).\n");
    printf("          --updateModel   : This option will allow users to use the gwas panel to update the\n");
    printf("                            parameter estimates (if found in the M3VCF file) during imputation.\n");
    printf("                            By default, if parameter estimates are found inthe M3VCF file, they\n");
    printf("                            will be used without being updated. This default setting should give\n");
    printf("                            users accurate enough estimates and this handle need NOT be used unless \n");
    printf("                            the user has strong reasons to do so.\n");
 // printf("\n       --onlyRefMarkers   : If ON, markers which were PRESENT in the gwas panel but ABSENT \n");
 //   printf("                            in reference panel will be omitted from the output. By default,\n");
//    printf("                            it is OFF and these markers, although NOT important for the impu-\n");
 //   printf("                            -tation, will be present in the output files\n");
 // printf("\n            --vcfOutput   : If ON, imputed data will be output as VCF file as well. [Default: ON]. \n");
    printf("           --doseOutput   : If ON, imputed data will be output as dosage file as well [Default: OFF].\n");
    printf("            --hapOutput   : If ON, phased imputed data will be output as well [Default: OFF]. \n");
    printf("               --format   : Specifies which fields to output for the FORMAT field in output \n");
    printf("                            VCF file. Available handles: GT,DS,GP [Default: GT,DS].\n");
  printf("\n               --chr 22   : Chromosome number for which we will carry out imputation.\n");
    printf("         --start 100000   : Start position for imputation by chunking.\n");
    printf("           --end 200000   : End position for imputation by chunking. \n");
    printf("         --window 20000   : Length of buffer region on either side of --start and --end.\n");
  
  printf("\n             --rounds 5   : Rounds of optimization for model parameters, which describe population \n");
    printf("                            recombination rates and per SNP error rates. By default = 5.\n");
    printf("           --states 200   : Maximum number of reference (or target) haplotypes to be examined  \n");
    printf("                            during parameter optimization. By default = 200.\n");
  
  printf("\n                  --rec   : Recombination file (.recom) from previous run of imputation. \n"); 
   printf("                  --err   : Error file (.erate) from previous run of imputation. \n");
    
    printf("               --cpus 5   : Number of cpus for parallel computing. Works only with Minimac3-omp.\n\n");




//
//      printf(" ");
//    printf(" ");
//    printf(" ");
//    printf(" ");
//    printf(" ");
//    printf(" ");
//    printf(" ");
//    printf(" ");
//    printf(" ");
//    printf(" ");
//    printf(" ");
//    printf(" ");
//    printf(" ");
//    printf(" ");


	return;
}


