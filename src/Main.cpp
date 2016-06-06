

#include <stdio.h>
#include <iostream>
#include <ctime>
#include "Parameters.h"
#include "StringBasics.h"
#include "HaplotypeSet.h"
#include "Unique.h"
#include "MarkovModel.h"
#include "MarkovParameters.h"
#include "Imputation.h"
#include "ImputationStatistics.h"
#include <unistd.h>
#include <fstream>
#include <ostream>
int transFactor = 3;
int cisFactor = 2;


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

	bool log = false, lowMemory=false, duplicates=false, unphasedOutput=false, phased = false,passOnly = false, doseOutput = false, vcfOutput = true, gzip = true, nobgzip = false, rsid=false;
	bool processReference=false,updateModel=false, typedOnly=false, help = false, params = false;
    String MyChromosome="";


	ParameterList inputParameters;
	PhoneHome::allThinning = 50;

	BEGIN_LONG_PARAMETERS(longParameterList)
		LONG_PARAMETER_GROUP("Reference Haplotypes")
		LONG_STRINGPARAMETER("refHaps", &refHaps)
		LONG_PARAMETER("passOnly", &passOnly)
		LONG_PARAMETER("rsid", &rsid)
		LONG_PARAMETER_GROUP("Target Haplotypes")
		LONG_STRINGPARAMETER("haps", &haps)
//		LONG_STRINGPARAMETER("snps", &snps)
		LONG_PARAMETER_GROUP("Output Parameters")
		LONG_STRINGPARAMETER("prefix", &outfile)
		LONG_PARAMETER("processReference", &processReference)
		LONG_PARAMETER("updateModel", &updateModel)
		LONG_PARAMETER("nobgzip", &nobgzip)
		LONG_PARAMETER("vcfOutput", &vcfOutput)
		LONG_PARAMETER("doseOutput", &doseOutput)
		LONG_PARAMETER("hapOutput", &phased)
		LONG_STRINGPARAMETER("format", &format)
		LONG_PARAMETER("allTypedSites", &typedOnly)
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
		LONG_PARAMETER("log", &log)
		LONG_PARAMETER("lowMemory", &lowMemory)
		LONG_PARAMETER("help", &help)
		LONG_INTPARAMETER("cpus", &cpus)
		LONG_PARAMETER("params", &params)
		LONG_PHONEHOME(VERSION)
		BEGIN_LEGACY_PARAMETERS()
		LONG_STRINGPARAMETER("MyChromosome", &MyChromosome)
//		LONG_INTPARAMETER("transFactor", &transFactor)
//		LONG_INTPARAMETER("cisFactor", &cisFactor)
		LONG_INTPARAMETER("sample", &max_indiv)
		LONG_INTPARAMETER("marker", &max_marker)
		LONG_PARAMETER("duplicates", &duplicates)
		LONG_PARAMETER("unphasedOutput", &unphasedOutput)
		END_LONG_PARAMETERS();


	inputParameters.Add(new LongParameters(" Command Line Options: ",longParameterList));

    String compStatus;
	inputParameters.Read(argc, &(argv[0]));

    FILE *LogFile=NULL;
    if(log)
        LogFile=freopen(outfile+".logfile","w",stdout);
    dup2(fileno(stdout), fileno(stderr));


    Minimac3Version();
	if (help)
	{
		helpFile();
		return(-1);
	}

	inputParameters.Status();

    #ifdef _OPENMP
        omp_set_num_threads(cpus);
    #else
        cpus=1;
    #endif


    if(nobgzip)
        gzip=false;

    cout<<endl<<endl;
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
        cout<<"       Program will only estimate parameters and create M3VCF file.\n";
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
        cout<<"       Program will NOT run if \"--refHaps\" is a VCF file !!!\n"<<endl;

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

        string formatPiece,formatTemp=format.c_str();
        char *end_str1;

        for(char * pch = strtok_r ((char*)formatTemp.c_str(),",", &end_str1);
            pch!=NULL;
            pch = strtok_r (NULL, ",", &end_str1))
        {

            formatPiece=(string)pch;
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



    if(lowMemory)
     {
        cout<<" Low Memory Version of Minimac3 initiated  !!! \n"<<endl;
    }



    HaplotypeSet target,reference;
    target.MyChromosome=(string)MyChromosome;
    reference.MyChromosome=(string)MyChromosome;
    reference.CPU=cpus;
    reference.Duplicates=duplicates;





	if(!processReference)
    {
        cout<<" ------------------------------------------------------------------------------"<<endl;
        cout<<"                       PRELIMINARY GWAS/TARGET FILE CHECK                      "<<endl;
        cout<<" ------------------------------------------------------------------------------"<<endl;


	    if (!target.BasicCheckForTargetHaplotypes(haps))
        {
            cout << "\n Program Exiting ... \n\n";
            compStatus="Target.Panel.Load.Error";
            PhoneHome::completionStatus(compStatus.c_str());
            return(-1);
        }

        cout<<" ------------------------------------------------------------------------------"<<endl;
        cout<<"                       PRELIMINARY REFERENCE FILE CHECK                        "<<endl;
        cout<<" ------------------------------------------------------------------------------"<<endl<<endl;


        cout << " Performing basic file check on Reference haplotype file ..." << endl;

        cout << "\n Checking File ..." << endl;

        if(reference.DetectReferenceFileType(refHaps).compare("NA")==0)
        {
            cout << "\n Program could NOT open file : " << refHaps << endl;
            cout << "\n Program Exiting ... \n\n";
            compStatus="Reference.Panel.Load.Error";
            PhoneHome::completionStatus(compStatus.c_str());
            return(-1);
        }
        cout << " File Exists ..." << endl;
        cout << "\n Checking File Format ..." << endl;

        if(reference.DetectReferenceFileType(refHaps).compare("Invalid")==0)
        {
            cout << "\n Reference File provided by \"--refHaps\" must be a VCF/M3VCF file !!! \n";
            cout << " Please check the following file : "<<refHaps<<endl;
            cout << "\n Program Exiting ... \n\n";
            compStatus="Reference.Panel.Load.Error";
            PhoneHome::completionStatus(compStatus.c_str());
            return(-1);
        }
        else if(reference.DetectReferenceFileType(refHaps).compare("m3vcf")==0)
        {
            cout<<"\n Reference File Format = M3VCF (Minimac3 VCF File) "<<endl;

            cout <<"\n NOTE: For M3VCF files, if parameter estimates are available in the file, \n";
               cout<<"       they will be used by default (RECOMMENDED !). If the user has reasons\n";
               cout<<"       to believe that updating the parameters would increase accuracy, they\n";
             cout<<"       should use handle \"--updateModel\" (not required in typical GWAS studies).\n";
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
                    cout<<"       value of \"--rounds\" will be ignored unless user has\n";
                    cout<<"       \"--updateModel\" ON (since, otherwise estimates are\n";
                    cout<<"       not going to be updated and value of \"--rounds\" would\n";
                    cout<<"       not make sense) !!!"<<endl;

                }
            }


        }
        else if(reference.DetectReferenceFileType(refHaps).compare("vcf")==0)
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



    reference.updateCoeffs(transFactor,cisFactor);


	if (!reference.FasterLoadHaplotypes(refHaps, max_indiv, max_marker,chr,start,end,window,rsid,processReference,passOnly,outfile,gzip))
	{
		cout << "\n Program Exiting ... \n\n";
		compStatus="Reference.Panel.Load.Error";
        PhoneHome::completionStatus(compStatus.c_str());
        return(-1);

	}


    cout<<endl;
	if(processReference)
    {
	    Imputation thisDataFast(reference, reference, outfile, errFile, recFile, phased
                             , gzip, rsid, rounds, states, vcfOutput, doseOutput
                             , typedOnly,formatVector,lowMemory);
        thisDataFast.createEstimates(reference, reference, reference.optEndPoints, true);
	}

	int time_load = time(0) - time_prev;
	time_prev = time(0);



    cout << "\n Time taken to load reference haplotype set = " << time_load << " seconds."<<endl<<endl;


	if(!processReference)
    {
        cout<<" ------------------------------------------------------------------------------"<<endl;
        cout<<"                          TARGET/GWAS HAPLOTYPE PANEL                         "<<endl;
        cout<<" ------------------------------------------------------------------------------"<<endl;


	    if (!target.LoadTargetHaplotypes(haps, snps, reference.markerName,reference,typedOnly,passOnly))
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
                             , gzip, rsid, rounds, states, vcfOutput, doseOutput
                             , typedOnly,formatVector,updateModel,unphasedOutput,lowMemory);

        thisDataFast.performImputation(target, reference, golden);
	}


    cout<<" ------------------------------------------------------------------------------"<<endl;
    cout<<"                                END OF PROGRAM                                 "<<endl;
    cout<<" ------------------------------------------------------------------------------"<<endl;

    time_load = time(0) - time_prev;
	int time_tot = time(0) - start_time;

    cout << "\n Program Successfully Implemented... \n ";


	printf("\n Total Run completed in %d hours, %d mins, %d seconds.\n",
		time_tot / 3600, (time_tot % 3600) / 60, time_tot % 60);

    cout<<"\n Thank You for using Minimac3 !!! "<<endl<<endl;

    if(log)
        fclose (LogFile);


    compStatus="Success";
    PhoneHome::completionStatus(compStatus.c_str());

	return 0;

}




void Minimac3Version()
{
	printf("\n\n -------------------------------------------------------------------------------- \n");
	printf("          Minimac3 - Fast Imputation Based on State Space Reduction HMM\n");
	printf(" --------------------------------------------------------------------------------\n");
    printf("           (c) 2014 - Sayantan Das, Christian Fuchsberger, David Hinds\n");
    printf("                             Mary Kate Wing, Goncalo Abecasis \n");
//	printf(" Version	: Undocumented Release\n");
//	printf(" Built		: sayantan\n\n");
	cout<<"\n Version: " << VERSION<< ";\n Built: " << DATE << " by " << USER << std::endl;
    printf("\n URL = http://genome.sph.umich.edu/wiki/Minimac3\n");



}

void helpFile()
{

    printf("\n\n\t  Minimac3 is a lower memory and more computationally efficient implementation of \"minimac2\".\n");


    printf("\t It is an algorithm for genotypic imputation that works on phased genotypes (say from MaCH).\n");
    printf("\t Minimac3 is designed to handle very large reference panels in a more computationally efficient \n");
    printf("\t way with no loss of accuracy. This algorithm analyzes only the unique sets of haplotypes in \n");
    printf("\t small genomic segments, thereby saving on time-complexity, computational memory but no loss\n");
    printf("\t in degree of accuracy.\n");

printf("\n\n ----------------------------------------------------------------------------------------- \n");
	printf("                            Minimac3 - List of Usage Options \n");
	printf(" -----------------------------------------------------------------------------------------\n\n");

    printf(" --------- Reference Haplotypes --------- \n");
  printf("\n              --refHaps   : VCF file or M3VCF file containing haplotype data for reference panel.\n");
    printf("             --passOnly   : This option only imports variants with FILTER = PASS.\n");
    printf("                 --rsid   : This option only imports RS ID of variants from ID column (if available).\n");



  printf("\n --------- GWAS Haplotypes --------- \n");
  printf("\n                 --haps   : File containing haplotype data for target (gwas) samples. Must be VCF \n");
    printf("                            file. Zipped versions allowed.\n");

  printf("\n --------- Output Parameters --------- \n");
  printf("\n               --prefix   : Prefix for all output files generated. By default: [Minimac3.Output]\n");
    printf("     --processReference   : This option will only convert an input VCF file to M3VCF format\n");
    printf("                            (maybe for a later run of imputation). If this option is ON, \n");
    printf("                            no imputation would be performed.\n");
    printf("              --nobgzip   : If ON, output files will NOT be gzipped.\n");
    printf("          --updateModel   : If ON, the GWAS panel will also be used to update the parameter \n");
    printf("                            estimates (if and when estimates are found in M3VCF files)\n");
    printf("            --vcfOutput   : If ON, imputed data will NOT be output as VCF dosage file [Default: ON].\n");
    printf("           --doseOutput   : If ON, imputed data will be output as MaCH dosage file    [Default: OFF].\n");
    printf("            --hapOutput   : If ON, phased imputed data will be output as well         [Default: OFF]. \n");
    printf("               --format   : Specifies which fields to output for the FORMAT field in output \n");
    printf("                            VCF file. Available handles: GT,DS,GP [Default: GT,DS].\n");
    printf("        --allTypedSites   : If ON, sites available ONLY in GWAS panel will also be output [Default: OFF]. \n");



  printf("\n --------- Subset Parameters --------- \n");
  printf("\n                  --chr   : Chromosome number for which to carry out imputation.\n");
    printf("                --start   : Start position for imputation by chunking.\n");
    printf("                  --end   : End position for imputation by chunking. \n");
    printf("               --window   : Length of buffer region on either side of --start and --end.\n");


  printf("\n --------- Estimation Parameters --------- \n");
  printf("\n               --rounds   : Rounds of optimization for model parameters, which describe population \n");
    printf("                            recombination rates and per SNP error rates. By default = 5.\n");
    printf("               --states   : Maximum number of reference (or target) haplotypes to be examined  \n");
    printf("                            during parameter optimization. By default = 200.\n");

  printf("\n --------- Starting Parameters --------- \n");
  printf("\n                  --rec   : Recombination estimates file from a previous run of Minimac3.\n");
    printf("                  --err   : Error estimates file from a previous run of Minimac3.\n");

  printf("\n --------- Other Parameters --------- \n");
  printf("\n                  --log   : If ON, log will be written to $prefix.logfile.\n");
    printf("                 --help   : If ON, detailed help on options and usage.\n");
    printf("            --lowMemory   : If ON, a low memory version of Minimac3 will be run.\n");
    printf("                 --cpus   : Number of cpus for parallel computing. Works only with Minimac3-omp.\n\n");


  printf("\n Please visit <http://genome.sph.umich.edu/wiki/Minimac3> for detailed documentation ...\n\n");
    cout<<endl;

	return;
}


