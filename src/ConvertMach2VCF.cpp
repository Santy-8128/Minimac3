

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
	String haps = "", snps = "";
	String outfile = "Minimac3.Output";
	String format = "GT,DS";
	String recFile = "", errFile = "",chr="",golden="";
	int cpus=1,start=0, end=0, window=0, max_indiv = 0, max_marker = 0, rounds=5, states=200;
    vector<bool> formatVector(3,false);
    #ifdef _OPENMP
    cpus=5;
    #endif

	bool phased = false, doseOutput = false, vcfOutput = false, gzip = true, nogzip = false, rsid=false;
	bool compressOnly=false, onlyRefMarkers=false, help = false, params = false;

	ParameterList inputParameters;
	PhoneHome::allThinning = 10;

	BEGIN_LONG_PARAMETERS(longParameterList)
		LONG_PARAMETER_GROUP("Reference Haplotypes")
		LONG_STRINGPARAMETER("refHaps", &refHaps)
		//LONG_PARAMETER("rsid", &rsid)
		LONG_PARAMETER_GROUP("Target Haplotypes")
		LONG_STRINGPARAMETER("haps", &haps)
		LONG_STRINGPARAMETER("snps", &snps)
		LONG_PARAMETER_GROUP("Output Parameters")
		LONG_STRINGPARAMETER("prefix", &outfile)
		LONG_PARAMETER("nogzip", &nogzip)
		LONG_PARAMETER("compressOnly", &compressOnly)
		LONG_PARAMETER("onlyRefMarkers", &onlyRefMarkers)
		LONG_PARAMETER("vcfOutput", &vcfOutput)
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
		LONG_STRINGPARAMETER("golden", &golden)
		LONG_INTPARAMETER("sample", &max_indiv)
		LONG_INTPARAMETER("marker", &max_marker)
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

    if(nogzip)
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


    if(compressOnly)
    {

        cout<<" NOTE: If \"--compressOnly\" is ON, all options under \"Target Haplotypes\", \"Starting Parameters\",";
        cout<<"       and \"Estimation Parameters\" will be ignored !!!\n";
        cout<<"       Program will only create optm file.\n";
        cout<<"       No imputation will be performed, hence other parameters are unnecessary !!!"<<endl<<endl;

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


    if(!compressOnly)
    {

         if (haps == "")
        {
            cout << "Missing \"--haps\", a required parameter.\n";
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
        if(recFile=="" && errFile!="")
        {
            cout << " Both \"--rec\" and \"--err\" parameters must be provided (or NONE) !!! \n";
            cout << " Only \"--err\" is provided. \n\n";
            compStatus="Command.Line.Error";
            PhoneHome::completionStatus(compStatus.c_str());
            return(-1);
        }
        if(errFile=="" && recFile!="")
        {
            cout << " Both \"--rec\" and \"--err\" parameters must be provided (or NONE) !!! \n";
            cout << " Only \"--rec\" is provided. \n\n";
            compStatus="Command.Line.Error";
            PhoneHome::completionStatus(compStatus.c_str());
            return(-1);
        }

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

    HaplotypeSet target;

	if(!compressOnly)
    {
        cout<<" ------------------------------------------------------------------------------"<<endl;
        cout<<"                       PRELIMNARY TARGET PANEL FILE CHECK                      "<<endl;
        cout<<" ------------------------------------------------------------------------------"<<endl;


	    if (!target.BasicCheckForTargetHaplotypes(haps))
        {
            cout << "\n Program Exiting ... \n\n";
            compStatus="Target.Panel.Load.Error";
            PhoneHome::completionStatus(compStatus.c_str());
            return(-1);
        }
    }

	int start_time = time(0);
	int time_prev = start_time;

    cout<<" ------------------------------------------------------------------------------"<<endl;
    cout<<"                           REFERENCE HAPLOTYPE PANEL                           "<<endl;
    cout<<" ------------------------------------------------------------------------------"<<endl;

	HaplotypeSet reference;
	if (!reference.FastLoadHaplotypes(refHaps, max_indiv, max_marker,chr,start,end,window,rsid))
	{
		cout << "\n Program Exiting ... \n\n";
		compStatus="Reference.Panel.Load.Error";
        PhoneHome::completionStatus(compStatus.c_str());
        return(-1);

	}

	int time_load = time(0) - time_prev;
	time_prev = time(0);

    if(reference.vcfType)
        reference.writeOptmFile(outfile,gzip);

    cout << "\n Time taken to load reference haplotype set = " << time_load << " seconds."<<endl<<endl;


//    MINIMIZED COST              = " << sum<< endl;
//    cout<<" FOLD INCREASE IN SPEEED     = "<< (double)rHap.numHaplotypes*(double)rHap.numMarkers /(double)sum<<endl;
//    cout<<" MAX_BLOCK_LENGTH            = "<< maxmax <<endl;


	if(!compressOnly)
    {
        cout<<" ------------------------------------------------------------------------------"<<endl;
        cout<<"                          TARGET/GWAS HAPLOTYPE PANEL                         "<<endl;
        cout<<" ------------------------------------------------------------------------------"<<endl;


	    if (!target.LoadTargetHaplotypes(haps, snps, reference.markerName,reference.VariantList))
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




    if(!compressOnly)
	{
	    Imputation thisDataFast(target, reference, outfile, errFile, recFile, phased
                             , gzip, rounds, states, vcfOutput, doseOutput
                             , onlyRefMarkers,formatVector);
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


	printf("\n  URL = https://statgen.sph.umich.edu/wiki/Minimac3\n");


    printf("\n\n\t  Minimac3 is a lower memory and more computationally efficient implementation of \"minimac\".\n");

    printf("\t It is an algorithm for genotypic imputation that works on phased genotypes (say from MaCH).\n");
    printf("\t The name is derived from \"Minimac eXpedited\". Minimac3 is designed to handle very large re-\n");
    printf("\t ference panels in a more computationally efficient way with no loss of accuracy. This algo-\n");
    printf("\t rithm analyzes only the unique sets of haplotypes in small genomic segments, thereby saving\n");
    printf("\t on time-complexity, computational memory but no loss in degree of accuracy.\n\n\n\n");

    printf(" The usage has the following main five segments (List of Usage Options Given Below):\n\n");


    printf("     Reference Haplotypes : Minimac3 can handle either VCF files or OPTM files as input for the\n");
    printf("                            reference panel. The program can itself identify the type of file,\n");
    printf("                            and no handle is necessary for that. OPTM files are customized files \n");
    printf("                            created by Minimac3 (possibly in some previous run) that stores large \n");
    printf("                            reference panels in a compact form so as to save memory and computation \n");
    printf("                            time involved in reading large files. OPTM files must be generated \n");
    printf("                            in some previous run of Minimac3 and can be saved and used in later \n");
    printf("                            runs for faster loading of data. See section on OPTM files and examples\n");
    printf("                            below on how to use them.\n");


  printf("\n                            \"--refHaps\" denotes the main input reference file which is either a VCF \n");
    printf("                            file or OPTM file. No handle is necessary for denoting type of file, \n");
    printf("                            program will detect it itself. \n\n");



    printf("     Target Haplotypes    : Minimac3 can handle either VCF files or MaCH files as input for\n");
    printf("                            the target/gwas data. The program can itself identify the type of \n");
    printf("                            file, and no handle is necessary for that. Note that input VCF files \n");
    printf("                            would be automatically assumed to be phased\n");
  printf("\n                            \"--haps\" denotes the main input target file which is either a VCF file\n");
    printf("                            (.vcf or .vcf.gz) or a MaCH phased file (.mach or .haps). The extensions\n");
    printf("                            are not mandatory. Zipped formats are supported.\n");
 printf("\n                            \"--snps\" denotes the marker name input file for target data panel. This \n");
    printf("                            panel parameter is mandatory if user is using MaCH phased files for the   \n");
    printf("                            target and will be ignored if VCF files are used. Markers which are in  \n");
    printf("                            the target panel and NOT in the reference panel would be excluded from \n");
    printf("                            the output files. User must merge these extra markers back to the original  \n");
    printf("                            data in order to analyze them. \n\n");

    printf("   List of Usage Options  :\n\n");

    printf("     --refHaps filename   : VCF file or OPTM file containing haplotype data for reference panel.\n");
    printf("                 --rsid   : This option specifies that the marker names should be read from ID \n");
    printf("                            column of VCF file provided in --refHaps.\n");
    printf("\n        --haps filename   : File containing haplotype data for target (gwas) samples. Can be VCF \n");
    printf("                            file or MaCH haplotype file. Zipped versions allowed.\n");
    printf("        --snps filename   : File containing marker information for the target (gwas) samples.\n");
    printf("                            It is mandatory for MaCH files but will be ignored for VCF files.\n");
    printf("\n        --prefix output   : Prefix for all output files generated. By default: [Minimac3.Output]\n");
    printf("                 --gzip   : If ON, output files will be gzipped.\n");
    printf("         --compressOnly   : This option will only convert an input VCF file to OPTM format\n");
    printf("                            (maybe for a later run of imputation). If this option is ON, \n");
    printf("                            no imputation would be performed and thus all other parameters \n");
    printf("                            will be ignored (of course, except for parameters on Reference \n");
    printf("                            Haplotypes and Subsetting Options).\n");
  printf("\n       --onlyRefMarkers   : If ON, markers which were PRESENT in the gwas panel but ABSENT \n");
    printf("                            in reference panel will be omitted from the output. By default,\n");
    printf("                            it is OFF and these markers, although NOT important for the impu-\n");
    printf("                            -tation, will be present in the output files\n");
  printf("\n            --vcfOutput   : If ON, imputed data will be output as VCF file as well. [Default: ON]. \n");
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


