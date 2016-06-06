#include "Imputation.h"


void Imputation::ImputeTraverse(HaplotypeSet &rHap,HaplotypeSet &tHap,int hapID,
                              MarkovModel &MM,int group, vector<float> &recomProb,
                              vector<float> &PrevRightFoldedProb,
                              vector<float> &CurrentRightProb,
                              vector<float> &CurrentNoRecoRightProb,
                              HaplotypeSet &rHapMain)
{

    PrevRightFoldedProb.resize(rHap.ReducedStructureInfo[group-1].RepSize);

    MM.foldProbabilities(PrevRightFoldedProb,group-1,
                            rHap.ReducedStructureInfo[group-1],1,rHap.numHaplotypes);

    MM.Impute(tHap,hapID,group-1,
                            PrevRightFoldedProb,CurrentRightProb,CurrentNoRecoRightProb,
                        rHap.ReducedStructureInfo[group-1],rHapMain.AlleleFreq);

    splitFoldedProb(recomProb,CurrentRightProb,CurrentNoRecoRightProb);
    MM.unfoldProbabilities(group-1,recomProb,CurrentNoRecoRightProb,
                           PrevRightFoldedProb,1,
                           rHap.ReducedStructureInfo,rHap.numHaplotypes);

}


void Imputation::EMTraverse(HaplotypeSet &rHap,HaplotypeSet &tHap,int hapID,
                              MarkovModel &MM,int group, vector<float> &recomProb,
                              vector<float> &PrevRightFoldedProb,
                              vector<float> &CurrentRightProb,
                              vector<float> &CurrentNoRecoRightProb,
                              HaplotypeSet &rHapMain)
{

    PrevRightFoldedProb.resize(rHap.ReducedStructureInfo[group-1].RepSize);

    MM.foldProbabilities(PrevRightFoldedProb,group-1,
                            rHap.ReducedStructureInfo[group-1],1,rHap.numHaplotypes);

    MM.CountExpected(tHap,hapID,group-1,
                        PrevRightFoldedProb,CurrentRightProb,CurrentNoRecoRightProb,
                        rHap.ReducedStructureInfo[group-1],rHapMain.AlleleFreq);



    splitFoldedProb(recomProb,CurrentRightProb,CurrentNoRecoRightProb);
    MM.unfoldProbabilities(group-1,recomProb,CurrentNoRecoRightProb,
                           PrevRightFoldedProb,1,
                           rHap.ReducedStructureInfo,rHap.numHaplotypes);

}




void Imputation::LeftTraverse(HaplotypeSet &rHap,HaplotypeSet &tHap,int hapID,
                              MarkovModel &MM,int group, vector<float> &recomProb,
                              HaplotypeSet &rHapMain)
{
    vector<int> &optStructure=rHapMain.optEndPoints;

    int Start=optStructure[group-1];
    int End=optStructure[group];
    vector<vector<float> > &ThisBlockLeftProb=LowMemory ? MM.ThisBlockLeftProb : MM.leftProb[group-1];

    MM.foldProbabilities(ThisBlockLeftProb[0],group-1,rHap.ReducedStructureInfo[group-1],
                            0,rHap.numHaplotypes);
    MM.CurrentLeftNoRecoProb=ThisBlockLeftProb[0];

    MM.WalkLeft(tHap,hapID,group-1,
                   rHap.ReducedStructureInfo[group-1],rHapMain.AlleleFreq);

    splitFoldedProb(recomProb,ThisBlockLeftProb[End-Start],MM.CurrentLeftNoRecoProb);
    MM.unfoldProbabilities(group-1,recomProb,MM.CurrentLeftNoRecoProb,
                                       ThisBlockLeftProb[0],0,
                                       rHap.ReducedStructureInfo,rHap.numHaplotypes);
}


double Imputation::CalculateLikelihood(HaplotypeSet &rHap,MarkovModel &MM)
{
    double tempLogLikelihoodValue=0.0;
    vector<vector<float> > &LastBlock=LowMemory? MM.ThisBlockLeftProb:MM.leftProb[rHap.NoBlocks-1];
    ReducedHaplotypeInfo &TempBlock=rHap.ReducedStructureInfo[rHap.NoBlocks-1];
    vector<float> &LastColumn=LastBlock[TempBlock.BlockSize-1];

    for(int tempVar=0;tempVar<TempBlock.RepSize;tempVar++)
    {
        tempLogLikelihoodValue+=LastColumn[tempVar];
    }

    tempLogLikelihoodValue=log(tempLogLikelihoodValue);
    tempLogLikelihoodValue-=(MM.NoPrecisionJumps*log(1e15));
    return tempLogLikelihoodValue;
}





MarkovParameters* Imputation::createEstimates(HaplotypeSet &rHap,HaplotypeSet &tHap,vector<int> &optStructure,bool NoTargetEstimation)
{
    MarkovParameters *MP=new MarkovParameters(rHap.numMarkers);


    cout<<" ------------------------------------------------------------------------------"<<endl;
    cout<<"                             PARAMETER ESTIMATION                              "<<endl;
    cout<<" ------------------------------------------------------------------------------"<<endl;

    int time_prev = time(0);

    rHap.CreateSummary();
    cout<<endl;


    if(rHap.Recom.size()>0)
    {
        printf(" Reading pre-calculated estimates located in M3VCF file ...\n");
        cout<<endl;
        MP->Recom=rHap.Recom;
        MP->Error=rHap.Error;

        if(!updateModel)
        {
            EstimationRounds=0;

        }
        else
        {
            printf(" Updating pre-calculated estimates ...\n");
            cout<<endl;
        }
    }

    if (EstimationRounds > 0)
    {
        printf(" Setting up Markov Model for Parameter Estimation...\n");
        cout<<endl;
    }
    else
    {
        if(rHap.vcfType || rHap.Recom.size()==0)
            printf(" Parameter Estimation being skipped ...\n");
        else
            printf(" Using pre-calculated estimates in Markov Model ...\n");

        cout<<endl;
    }

    if(recFile!="")
    {
        MP->ReadCrossoverRates(recFile);
    }


    if(errFile!="")
    {
        MP->ReadErrorRates(errFile);
    }

    if(recFile!="" && errFile!="")
    {
        if(!updateModel)
        {
            EstimationRounds=0;
        }
        else
        {
            printf("\n Updating read estimates ...\n");
            cout<<endl;
        }

    }


    if (EstimationRounds > 0)
        {
            printf(" Initializing Model Parameters (using %s and up to %d haplotypes) ...",
               em ? "E-M" : "MCMC", EstimationStates);
        cout<<endl;
        }

    double LogLikelihoodValue=0.0;
    for (int round = 0; round < EstimationRounds; round++)
    {
        cout<<"\n Round "<<round+1<<" of Parameter Refinement ...";
        cout<<endl;
        int iterations = EstimationStates < rHap.numHaplotypes ? EstimationStates :  rHap.numHaplotypes;
        LogLikelihoodValue=0.0;
        #pragma omp parallel for
        for (int i = 0; i < iterations; i++)
        {
            // Reference leave one out (loo) panel
            HaplotypeSet rHap_loo;
            LooOptimalStructure(i,rHap,rHap_loo);

           // Target with one panel at position i
            HaplotypeSet tHap_loo;
            tHap_loo.Create(i,rHap);

            // Create Markov Model
            vector<float> PrevRightFoldedProb,CurrentRightProb,CurrentNoRecoRightProb,recomProb;
            PrevRightFoldedProb.reserve(rHap.maxRepSize);
            CurrentRightProb.reserve(rHap.maxRepSize);
            CurrentNoRecoRightProb.reserve(rHap.maxRepSize);
            recomProb.reserve(rHap.maxRepSize);


            vector<bool> missing_loo(markerCount,false);
                     MarkovModel MM(tHap_loo,rHap_loo,missing_loo,rHap.major,LowMemory);
            MM.CopyParameters(MP);
            MM.initializeMatrices(rHap_loo,tHap_loo);
            MM.ReinitializeMatrices();


            if(!missing_loo[0])
            {
                if(!tHap_loo.getMissingScaffoldedHaplotype(0,0))
                {
                    ConditionJunctionProb(rHap_loo,0,MM.junctionLeftProb[0],
                                          MM.Error[0],
                                          tHap_loo.getScaffoldedHaplotype(0,0)? rHap.AlleleFreq[0] : 1-rHap.AlleleFreq[0],
                                          tHap_loo.getScaffoldedHaplotype(0,0),
                                          MM.backgroundError,
                                          rHap_loo.ReducedStructureInfo[0]);
                            }
            }

            for(int group=1;group<=rHap.NoBlocks;group++)
            {
                LeftTraverse(rHap_loo,tHap_loo,0,MM,group,recomProb,rHap);
            }

            double tempLogLikelihoodValue=CalculateLikelihood(rHap_loo,MM);

            if(em)
            {
                for(int group=rHap.NoBlocks;group>0;group--)
                {
                    EMTraverse(rHap_loo,tHap_loo,0,MM,group,
                              recomProb,PrevRightFoldedProb,
                              CurrentRightProb,CurrentNoRecoRightProb,rHap);
                }
            MM.empiricalCount++;
            }

            #pragma omp critical
            {
                (*MP)+=MM;
                LogLikelihoodValue+=tempLogLikelihoodValue;
            }

        }

        if(!NoTargetEstimation && round>=EstimationRounds/2)
        {
            iterations = EstimationStates < (int)tHap.numHaplotypes ? EstimationStates : (int)tHap.numHaplotypes;

            #pragma omp parallel for
            for (int i = 0; i < iterations; i++)
            {
                MarkovModel MM(tHap,rHap,tHap.missing,rHap.major,LowMemory);
                vector<float> PrevRightFoldedProb,CurrentRightProb,CurrentNoRecoRightProb,recomProb;

                MM.CopyParameters(MP);
                //initializeMatrices(tHap,rHap,optStructure,rHap.ReducedStructureInfo);
                MM.initializeMatrices(rHap,tHap);
                MM.ReinitializeMatrices();

                if(!tHap.missing[0])
                {
                    if(!tHap.getMissingScaffoldedHaplotype(i,0))
                    {

                        ConditionJunctionProb(rHap,0,MM.junctionLeftProb[0],
                                              MM.Error[0],
                                              tHap.getScaffoldedHaplotype(i,0)? rHap.AlleleFreq[0] : 1-rHap.AlleleFreq[0],
                                              tHap.getScaffoldedHaplotype(i,0),
                                              MM.backgroundError,
                                              rHap.ReducedStructureInfo[0]);

                    }
                }


                for(int group=1;group<=rHap.NoBlocks;group++)
                {
                    LeftTraverse(rHap,tHap,i,MM,group,recomProb,rHap);
                }

                double tempLogLikelihoodValue=CalculateLikelihood(rHap,MM);
                if(em)
                {
                    for(int group=rHap.NoBlocks;group>0;group--)
                    {
                        EMTraverse(rHap,tHap,i,MM,group,
                              recomProb,PrevRightFoldedProb,
                              CurrentRightProb,CurrentNoRecoRightProb,rHap);
                    }
                    MM.empiricalCount++;
                }

                #pragma omp critical
                 {
                    (*MP)+=MM;
                    LogLikelihoodValue+=tempLogLikelihoodValue;
                }
            }
        }


        MP->UpdateModel();

        double crossovers = 0;
        for (int i = 0; i < rHap.numMarkers - 1; i++)
            crossovers += MP->Recom[i];

        double errors = 0;
        for (int i = 0; i <  rHap.numMarkers ; i++)
        {
            double heterozygosity = 1.0 - pow(rHap.AlleleFreq[i],2)
                                         - pow(1-rHap.AlleleFreq[i],2);
            errors += MP->Error[i] * heterozygosity;
        }
        errors /= (double) rHap.numMarkers  + 1e-30;

        printf("      %.0f mosaic crossovers expected per haplotype\n", crossovers);
        printf("      %.3g errors in mosaic expected per marker\n", errors);
        printf("    Log-Likelihood of this Iteration : %.5g", LogLikelihoodValue);
        cout<<endl;

        rHap.Recom=MP->Recom;
        rHap.Error=MP->Error;

    }



    printf("\n Saving estimated parameters for future use/reference to ...");
    cout<<endl;
    MP->WriteParameters(rHap.markerName, outFile, false);
    int time_load = time(0) - time_prev;

    cout << "\n Time taken for parameter estimation = " << time_load << " seconds. "<<endl;


    if(rHap.vcfType || EstimationRounds>0)
    {
        std::cout << "\n Writing final reduced haplotype information to [.m3vcf] file  : " <<outFile+".m3vcf" + (gzip ? ".gz" : "")<< endl<<endl;
        rHap.writem3vcfFile(outFile,gzip);

        if(rHap.vcfType)
       {
            remove(outFile+".draft"+".m3vcf" + (gzip ? ".gz" : ""));
            std::cout << "\n Temporary Draft [.m3vcf] file "<<outFile+".draft"+".m3vcf" + (gzip ? ".gz" : "")<<" deleted ... "<< endl<<endl;
       }

    }


    return MP;


}


void Imputation::performImputation(HaplotypeSet &tHap,HaplotypeSet &rHap, String Golden)
{

    vector<int> optStructure=rHap.optEndPoints;
    int time_prev = time(0),time_load,vcfSampleIndex=0;;
    includeGwas=true;


    MarkovParameters* MP=createEstimates(rHap,tHap,rHap.optEndPoints,1-includeGwas);

    cout<<" ------------------------------------------------------------------------------"<<endl;
    cout<<"                              MAIN IMPUTATION                                  "<<endl;
    cout<<" ------------------------------------------------------------------------------"<<endl;


    ImputationStatistics stats(rHap.numMarkers );
    IFILE dosages=NULL, hapdose=NULL, haps=NULL,vcfdosepartial=NULL;
    HaplotypeSet DosageForVcfPartial;
    DosageForVcfPartial.unphasedOutput=unphasedOutput;
    DosageForVcfPartial.TypedOnly=tHap.TypedOnly;
    DosageForVcfPartial.GWASOnlycounter=tHap.GWASOnlycounter;

    if(tHap.TypedOnly)
    {
        printf("\n Calculating Allele Frequency for Typed-Only variants ... ");
        cout<<endl;
        tHap.CalculateGWASOnlyFreq();

    }

    cout << "\n Starting Imputation ...";
    printf("\n\n Setting up Markov Model for Imputation ...");
    cout<<endl<<endl;


    if (phased && !unphasedOutput)
    {

        hapdose = ifopen(outFile + ".hapDose" + (gzip ? ".gz" : ""), "wb", gzip ?InputFile::BGZF:InputFile::UNCOMPRESSED);
        haps = ifopen(outFile + ".hapLabel" + (gzip ? ".gz" : ""), "wb", gzip ?InputFile::BGZF:InputFile::UNCOMPRESSED);

    }

    int maxVcfSample=200,NumVcfWritten=0,NumVcfCreated=0,NovcfParts=1;

    if((maxVcfSample)>=tHap.numSamples)
        maxVcfSample=tHap.numSamples;

    if(vcfOutput)
    {


        vcfdosepartial = ifopen(outFile + ".dose.vcf" + (gzip ? ".gz" : ""), "wb", gzip ?InputFile::BGZF:InputFile::UNCOMPRESSED);
        ifprintf(vcfdosepartial,"##fileformat=VCFv4.1\n");
        time_t t = time(0);
        struct tm * now = localtime( & t );
        ifprintf(vcfdosepartial,"##filedate=%d.%d.%d\n",(now->tm_year + 1900),(now->tm_mon + 1) ,now->tm_mday);
        ifprintf(vcfdosepartial,"##source=Minimac3\n");
        ifprintf(vcfdosepartial,"##contig=<ID=%s>\n",rHap.finChromosome.c_str());

        if(GT)
                ifprintf(vcfdosepartial,"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
        if(tHap.AllMaleTarget)
        {
            if(DS)
                ifprintf(vcfdosepartial,"##FORMAT=<ID=DS,Number=1,Type=Float,Description=\"Estimated Alternate Allele Dosage (For Male Chr: X) : [P(Alt Allele)]\">\n");
            if(GP)
                ifprintf(vcfdosepartial,"##FORMAT=<ID=GP,Number=2,Type=Float,Description=\"Estimated Posterior Probabilities for Genotypes 0 and 1 (For Male Chr: X) \">\n");
        }
        else
        {
            if(DS)
                ifprintf(vcfdosepartial,"##FORMAT=<ID=DS,Number=1,Type=Float,Description=\"Estimated Alternate Allele Dosage : [P(0/1)+2*P(1/1)]\">\n");
            if(GP)
                ifprintf(vcfdosepartial,"##FORMAT=<ID=GP,Number=3,Type=Float,Description=\"Estimated Posterior Probabilities for Genotypes 0/0, 0/1 and 1/1 \">\n");
        }


        ifprintf(vcfdosepartial,"##INFO=<ID=MAF,Number=1,Type=Float,Description=\"Estimated Alternate Allele Frequency\">\n");
        ifprintf(vcfdosepartial,"##INFO=<ID=R2,Number=1,Type=Float,Description=\"Estimated Imputation Accuracy\">\n");
        ifprintf(vcfdosepartial,"##INFO=<ID=ER2,Number=1,Type=Float,Description=\"Empirical (Leave-One-Out) R-square (available only for genotyped variants)\">\n");
        ifprintf(vcfdosepartial,"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
        ifclose(vcfdosepartial);

        if(!tHap.AllMaleTarget)
            DosageForVcfPartial.InitializePartialDosageForVcfOutput((2*maxVcfSample),rHap.numMarkers,format);
        else
            DosageForVcfPartial.InitializePartialDosageForVcfOutputMaleSamples(maxVcfSample<MaxSample?maxVcfSample:MaxSample,rHap.numMarkers,format);
    }

    if(doseOutput)
        dosages = ifopen(outFile + ".dose" + (gzip ? ".gz" : ""), "wb",(gzip ? InputFile::BGZF:InputFile::UNCOMPRESSED) );


    #pragma omp parallel for
    for(int hapId=0;hapId<MaxSample;hapId++)
    {


        if (hapId %2==1)
        {
            if(rHap.finChromosome!="X")
                continue;
            else if(!tHap.AllMaleTarget)
                continue;
        }

        vector<float> foldedProb,noRecomProb, rightProb,FirstHapDoseForVCF;
        vector<bool> padded(rHap.numMarkers),paddedMissing(rHap.numMarkers), FirstHapAlleleForVCF;
        vector<float> PrevRightFoldedProb,CurrentRightProb,CurrentNoRecoRightProb,recomProb;

        MarkovModel MM(tHap,rHap,tHap.missing,rHap.major,LowMemory);
        MM.CopyParameters(MP);
        MM.initializeMatrices(rHap,tHap);

        int hapIdIndiv=hapId;
        do{

            printf("  Processing Haplotype %d of %d ...", hapIdIndiv + 1, MaxSample);
            cout<<endl;

            MM.ThisHapId=hapIdIndiv;
            MM.ReinitializeMatrices();

            if(!tHap.missing[0])
            {
                if(!tHap.getMissingScaffoldedHaplotype(hapIdIndiv,0))
                {

                    ConditionJunctionProb(rHap,0,MM.junctionLeftProb[0],
                                          MM.Error[0],
                                          tHap.getScaffoldedHaplotype(hapIdIndiv,0)? rHap.AlleleFreq[0] : 1-rHap.AlleleFreq[0],
                                          tHap.getScaffoldedHaplotype(hapIdIndiv,0),
                                          MM.backgroundError,
                                          rHap.ReducedStructureInfo[0]);
                }
            }

            for(int group=1;group<=rHap.NoBlocks;group++)
            {
                LeftTraverse(rHap,tHap,hapIdIndiv,MM,group,recomProb,rHap);
            }

            for(int group=rHap.NoBlocks;group>0;group--)
            {
                ImputeTraverse(rHap,tHap,hapIdIndiv,MM,group,
                              recomProb,PrevRightFoldedProb,CurrentRightProb,CurrentNoRecoRightProb,rHap);
            }

            for(int jjj=0;jjj<rHap.numMarkers;jjj++)
            {
                padded[jjj]=tHap.getScaffoldedHaplotype(hapIdIndiv,jjj);
                paddedMissing[jjj]=tHap.getMissingScaffoldedHaplotype(hapIdIndiv,jjj);
            }

            if(vcfOutput)
            {
                if(hapIdIndiv%2==0)
                {
                   FirstHapDoseForVCF= MM.imputedHap;
                   FirstHapAlleleForVCF= MM.imputedAlleleNumber;
                }
            }
            #pragma omp critical
            {
                stats.Update(MM.imputedHap, MM.leaveOneOut,padded,paddedMissing,rHap.major);
            }

            #pragma omp critical
            if (phased && !unphasedOutput)
            {

                PrintHaplotypeData(rHap, tHap, hapdose, haps,
                                    MM.imputedHap, MM.imputedAlleleNumber,
                                    hapIdIndiv, tHap.AllMaleTarget?hapId:hapId/2);
            }


            if(tHap.AllMaleTarget)
                break;
            hapIdIndiv++;
        }while(hapIdIndiv<MaxSample && hapIdIndiv%2==1);

        #pragma omp critical
        if(doseOutput)
        {
            PrintDosageData(rHap, tHap, dosages, MM.imputedDose, tHap.AllMaleTarget?hapId:hapId/2);
        }
         #pragma omp critical
        if(vcfOutput)
        {

            printf("    Saving Individual %s for VCF File...\n",  tHap.individualName[tHap.AllMaleTarget?hapId:hapId/2].c_str());
            if(!tHap.AllMaleTarget)
                DosageForVcfPartial.SaveDosageForVcfOutputSampleWise(NumVcfCreated-NumVcfWritten,
                                                                 tHap.individualName[tHap.AllMaleTarget?hapId:hapId/2],
                                                                 FirstHapDoseForVCF,MM.imputedHap,
                                                                 FirstHapAlleleForVCF,MM.imputedAlleleNumber);
            else
                DosageForVcfPartial.SaveDosageForVcfOutputSampleWiseChrX(NumVcfCreated-NumVcfWritten,
                                                                 tHap.individualName[tHap.AllMaleTarget?hapId:hapId/2],
                                                                  MM.imputedHap,
                                                                 MM.imputedAlleleNumber);

            if(DosageForVcfPartial.TypedOnly)
            {

                DosageForVcfPartial.SaveIndexForGWASOnlyForVcfOutput(NumVcfCreated-NumVcfWritten,
                                                                     tHap.AllMaleTarget?hapId:hapId/2);
            }



            NumVcfCreated++;
            vcfSampleIndex++;

            if(NumVcfCreated%maxVcfSample==0 || NumVcfCreated==(tHap.AllMaleTarget?MaxSample:MaxSample/2))
            {

                string PartialVcfFileName(outFile),tempFileIndex1(outFile);
                stringstream strs;
                strs<<(NovcfParts);


                PartialVcfFileName+=(".dose.vcf.part." +
                                      (string)(strs.str())
                                     +(gzip ? ".gz" : ""));
                if(!tHap.AllMaleTarget)
                    printf("\n    --->>> Saving samples %d-%d in VCF file : %s ...\n\n",
                       (NumVcfWritten)+1,(MaxSample/2<(NumVcfWritten+maxVcfSample)?MaxSample/2:(NumVcfWritten+maxVcfSample)),
                       PartialVcfFileName.c_str());
                else
                    printf("\n    --->>> Saving samples %d-%d in VCF file : %s ...\n\n",
                       (NumVcfWritten)+1,(MaxSample<(NumVcfWritten+maxVcfSample)?MaxSample:(NumVcfWritten+maxVcfSample)),
                       PartialVcfFileName.c_str());


                FlushPartialVcf(rHap,tHap,DosageForVcfPartial,PartialVcfFileName,NovcfParts);
                if(NumVcfCreated<(tHap.AllMaleTarget?MaxSample:MaxSample/2))
                {
                    NovcfParts++;
                    NumVcfWritten+=maxVcfSample;

                    if(!tHap.AllMaleTarget)
                        DosageForVcfPartial.InitializePartialDosageForVcfOutput(maxVcfSample<(MaxSample/2-NumVcfWritten)?2*maxVcfSample:2*(MaxSample/2-NumVcfWritten),rHap.numMarkers,format);
                    else
                        DosageForVcfPartial.InitializePartialDosageForVcfOutputMaleSamples(maxVcfSample<(MaxSample-NumVcfWritten)?maxVcfSample:(MaxSample-NumVcfWritten),rHap.numMarkers,format);

                }
            }
        }
    }

    cout<<endl<<" Imputation Finished ... "<<endl;


    if (phased && !unphasedOutput)
    {
        ifclose(hapdose);
        ifclose(haps);

        cout<<endl<<" Haplotype Dosage information written to : "<<
            outFile + ".hapDose" + (gzip ? ".gz" : "")<<endl;
        cout<<endl<<" Haplotype Allele information written to : "<<
        outFile + ".hapLabel" + (gzip ? ".gz" : "")<<endl;
    }



    if(doseOutput)
    {
        ifclose(dosages);
        cout<<endl<<" Dosage information written to           : "<<
        outFile + ".dose" + (gzip ? ".gz" : "")<<endl;
    }

    PrintInfoFile(rHap,tHap,stats);

    time_load = time(0) - time_prev;
    cout << "\n Time taken for imputation = " << time_load << " seconds."<<endl<<endl;


    if(vcfOutput)
        MergeFinalVcfAllVariants(rHap,tHap,stats,NovcfParts);

}








void Imputation::MergeFinalVcfAllVariants(HaplotypeSet &rHap,HaplotypeSet &tHap,ImputationStatistics &stats,int MaxIndex)
{
    cout<<" ------------------------------------------------------------------------------"<<endl;
    cout<<"                                FINAL VCF MERGE                                "<<endl;
    cout<<" ------------------------------------------------------------------------------"<<endl;

    printf("\n Merging partial VCF files to final output VCF File :  %s ",(outFile + ".dose.vcf" + (gzip ? ".gz" : "")).c_str() );
    cout<<endl<<endl;

    IFILE vcfdosepartial = ifopen(outFile + ".dose.vcf" + (gzip ? ".gz" : ""),  "a", gzip ?InputFile::BGZF:InputFile::UNCOMPRESSED);

    vector<IFILE> vcfdosepartialList(MaxIndex);

    for(int i=1;i<=MaxIndex;i++)
    {
        string tempFileIndex(outFile);
        stringstream strs;
        strs<<(i);
        tempFileIndex+=(".dose.vcf.part." +
                         (string)(strs.str())+(gzip ? ".gz" : ""));
        vcfdosepartialList[i-1] = ifopen(tempFileIndex.c_str(), "r");
    }
    string line;
    for(int i=1;i<=MaxIndex;i++)
    {
        line.clear();
        vcfdosepartialList[i-1]->readLine(line);
        ifprintf(vcfdosepartial,"%s",line.c_str());
    }

    int i=0;
    for (int index =0; index < rHap.RefTypedTotalCount; index++)
    {

        if(index%10000==0)
        {
            printf("    Merging marker %d of %d [%.1f%%] to VCF File ...", index + 1, rHap.RefTypedTotalCount,100*(double)(index + 1)/(int)rHap.RefTypedTotalCount);
            cout<<endl;
        }


        if(rHap.RefTypedIndex[index]==-1)
        {

            if(i>=rHap.PrintStartIndex && i <= rHap.PrintEndIndex)
            {

                ifprintf(vcfdosepartial,"\n%s\t%d\t%s\t%s\t%s\t.\tPASS\tMAF=%.5f;R2=%.5f",
                rHap.VariantList[i].chr.c_str(),rHap.VariantList[i].bp,
                RsId?rHap.VariantList[i].rsid.c_str():rHap.VariantList[i].name.c_str(),rHap.VariantList[i].refAlleleString.c_str(),
                rHap.VariantList[i].altAlleleString.c_str(),stats.AlleleFrequency(i) > 0.5 ? 1.0 - stats.AlleleFrequency(i) : stats.AlleleFrequency(i),stats.Rsq(i));


                if(!tHap.missing[i])
                    ifprintf(vcfdosepartial,";ER2=%.5f",stats.EmpiricalRsq(i));

                ifprintf(vcfdosepartial,"\t%s",GT?(DS?(GP?"GT:DS:GP":"GT:DS"):(GP?"GT:GP":"GT")):(DS?(GP?"DS:GP":"DS"):(GP?"GP":"")));

                for(int j=1;j<=MaxIndex;j++)
                {
                    string tempFileIndex(outFile);
                    stringstream strs;
                    strs<<(j);
                    tempFileIndex+=(".dose.vcf.part."
                                    + (string)(strs.str())
                                    +(gzip ? ".gz" : ""));
                    line.clear();
                    vcfdosepartialList[j-1]->readLine(line);
                    ifprintf(vcfdosepartial,"%s",line.c_str());
                }

            }



            i++;
        }
        else
        {


            variant ThisTypedVariant =tHap.TypedOnlyVariantList[rHap.RefTypedIndex[index]];
            ifprintf(vcfdosepartial,"\n%s\t%d\t%s\t%s\t%s\t.\tPASS\t",
                     ThisTypedVariant.chr.c_str(),
                     ThisTypedVariant.bp,
                     RsId? ThisTypedVariant.rsid.c_str():ThisTypedVariant.name.c_str(),
                     ThisTypedVariant.refAlleleString.c_str(),
                     ThisTypedVariant.altAlleleString.c_str());


            ifprintf(vcfdosepartial,"GENOTYPED_ONLY;AN=%d;MAF=%.5f",
                     tHap.TotalSample[rHap.RefTypedIndex[index]],
                     tHap.AlleleFreq[rHap.RefTypedIndex[index]]);

//cout<<rHap.RefTypedIndex[index]<<" " <<tHap.TotalSample[rHap.RefTypedIndex[index]]<<" " << tHap.AlleleFreq[rHap.RefTypedIndex[index]]/(double)tHap.TotalSample[rHap.RefTypedIndex[index]]<< endl;

            ifprintf(vcfdosepartial,"\t%s",GT?(DS?(GP?"GT:DS:GP":"GT:DS"):(GP?"GT:GP":"GT")):(DS?(GP?"DS:GP":"DS"):(GP?"GP":"")));

            for(int j=1;j<=MaxIndex;j++)
            {
                string tempFileIndex(outFile);
                stringstream strs;
                strs<<(j);
                tempFileIndex+=(".dose.vcf.part."
                                + (string)(strs.str())
                                +(gzip ? ".gz" : ""));
                line.clear();
                vcfdosepartialList[j-1]->readLine(line);
                ifprintf(vcfdosepartial,"%s",line.c_str());

            }

//            ifprintf(vcfdosepartial,"\n");



        }

    }


    for(int i=1;i<=MaxIndex;i++)
    {
        ifclose(vcfdosepartialList[i-1]);
        string tempFileIndex(outFile);
        stringstream strs;
        strs<<(i);
        tempFileIndex+=(".dose.vcf.part." +
                        (string)(strs.str())+
                        (gzip ? ".gz" : ""));
        remove(tempFileIndex.c_str());
    }

    ifclose(vcfdosepartial);

    printf("\n Merging Finished ..." );
    cout<<endl <<endl;
}





void Imputation::FlushPartialVcf(HaplotypeSet &rHap,HaplotypeSet &tHap,HaplotypeSet &PartialDosage, string &filename,int &Index)
{

    string tempFileIndex(outFile),tempFileIndex1(outFile);
    IFILE vcfdosepartial = ifopen(filename.c_str(), "wb", InputFile::BGZF);

    for(int hapId=0;hapId<(int)PartialDosage.individualName.size();hapId++)
    {
        ifprintf(vcfdosepartial,"\t%s",PartialDosage.individualName[hapId].c_str());
    }
    ifprintf(vcfdosepartial,"\n");

    int i=0;
    for (int index =0; index < rHap.RefTypedTotalCount; index++)
    {

        if(rHap.RefTypedIndex[index]==-1)
        {

            if(i>=rHap.PrintStartIndex && i <= rHap.PrintEndIndex)
            {
                bool majorIsReference=false;
                if(!rHap.major[i])
                    majorIsReference=true;

                if(!tHap.AllMaleTarget)
                    PartialDosage.PrintDosageForVcfOutputForID(vcfdosepartial,i, majorIsReference,rHap.VariantList[i].refAllele);
                else
                    PartialDosage.PrintDosageForVcfOutputForIDMaleSamples(vcfdosepartial,i, majorIsReference,rHap.VariantList[i].refAllele);

                ifprintf(vcfdosepartial,"\n");

            }
            i++;

        }
        else
        {


            if(!tHap.AllMaleTarget)
                PartialDosage.PrintDosageGWASOnlyForVcfOutputForID
                (tHap,vcfdosepartial,rHap.RefTypedIndex[index]);
            else
                PartialDosage.PrintDosageGWASOnlyForVcfOutputForIDMaleSamples
                (tHap,vcfdosepartial,rHap.RefTypedIndex[index]);
            ifprintf(vcfdosepartial,"\n");
        }

    }

    ifclose(vcfdosepartial);



}



void Imputation::PrintDosageData(HaplotypeSet &rHap,HaplotypeSet &tHap,
                                 IFILE dosages, vector<float> &ThisDosage,
                                 int ThisSampleId)
{

    printf("    Outputting Individual %s for Dosage file...",  tHap.individualName[ThisSampleId].c_str());
    cout<<endl;
    ifprintf(dosages, "%s\tDOSE",tHap.individualName[ThisSampleId].c_str());
    int i=0;
    for (int index =0; index < rHap.RefTypedTotalCount; index++)
    {

        if(rHap.RefTypedIndex[index]==-1)
        {

            if(i>=rHap.PrintStartIndex && i <= rHap.PrintEndIndex)
            {

                 ifprintf(dosages, "\t%.3f", ThisDosage[i]);
            }
            i++;

        }
        else
        {
            int MarkerIndex=rHap.RefTypedIndex[index];
            bool a1,a2;
            double outAllele1=0.0,outAllele2=0.0;


            if(tHap.AllMaleTarget)
            {
                a1=tHap.GWASOnlyMissingSampleUnscaffolded[ThisSampleId][MarkerIndex];

                if(a1)
                {
                    outAllele1=tHap.AlleleFreq[MarkerIndex];
                    a1=round(outAllele1)==1?true:false;
                }
                else
                {
                    a1=tHap.GWASOnlyhaplotypesUnscaffolded[ThisSampleId][MarkerIndex];
                    if(a1)
                        outAllele1=1.0;
                }

                ifprintf(dosages, "\t%.3f", outAllele1);
            }
            else
            {
                a1=tHap.GWASOnlyMissingSampleUnscaffolded[2*ThisSampleId][MarkerIndex];
                a2=tHap.GWASOnlyMissingSampleUnscaffolded[2*ThisSampleId+1][MarkerIndex];

                if(a1 || a2)
                {
                    outAllele1=tHap.AlleleFreq[MarkerIndex];
                    outAllele2=outAllele1;
                    a1=round(outAllele1)==1?true:false;
                    a2=a1;
                }
                else
                {
                    a1=tHap.GWASOnlyhaplotypesUnscaffolded[2*ThisSampleId][MarkerIndex];
                    a2=tHap.GWASOnlyhaplotypesUnscaffolded[2*ThisSampleId+1][MarkerIndex];
                    if(a1)
                        outAllele1=1.0;
                    if(a2)
                        outAllele2=1.0;

                }

                ifprintf(dosages, "\t%.3f", outAllele1+outAllele2);
            }

        }

    }

    ifprintf(dosages,"\n");

}



void Imputation::PrintHaplotypeData(HaplotypeSet &rHap,HaplotypeSet &tHap,
                                 IFILE hapdose, IFILE haps,
                                 vector<float> &ThisimputedHap,vector<bool> ThisimputedAlleles,
                                 int ThisHapId, int ThisSampleId)
{

    char labels[]= {0, 'A', 'C', 'G', 'T', 'D', 'I', 'R'};


    printf("    Outputting HAPLO%d of Individual %s for Haplotype File...",
           tHap.AllMaleTarget?1:(ThisHapId%2+1) ,tHap.individualName[ThisSampleId].c_str());
    cout<<endl;
    ifprintf(hapdose, "%s\tHAPLO%d",  tHap.individualName[ThisSampleId].c_str(), tHap.AllMaleTarget?1:(ThisHapId%2+1) );
    ifprintf(haps, "%s\tHAPLO%d\t", tHap.individualName[ThisSampleId].c_str(), tHap.AllMaleTarget?1:(ThisHapId%2+1) );
    int i=0;
    for (int index =0; index < rHap.RefTypedTotalCount; index++)
    {

        if(rHap.RefTypedIndex[index]==-1)
        {

            if(i>=rHap.PrintStartIndex && i <= rHap.PrintEndIndex)
            {
                ifprintf(hapdose, "\t%.5f", ThisimputedHap[i]);
                ifprintf(haps, "%c", labels[(int) (ThisimputedAlleles[i]?
                        rHap.VariantList[i].altAllele
                            :rHap.VariantList[i].refAllele)]);

            }
            i++;
        }
        else
        {
            int MarkerIndex=rHap.RefTypedIndex[index];
            bool a1;
            a1=tHap.GWASOnlyMissingSampleUnscaffolded[ThisHapId][MarkerIndex];
            double outAllele1=0.0;

            if(a1)
            {
                outAllele1=tHap.AlleleFreq[MarkerIndex];
                a1=round(outAllele1)==1?true:false;
            }
            else
            {
                a1=tHap.GWASOnlyhaplotypesUnscaffolded[ThisHapId][MarkerIndex];
                 if(a1)
                        outAllele1=1.0;
            }

            ifprintf(haps, "%c", labels[(int) (a1?
                     tHap.TypedOnlyVariantList[MarkerIndex].altAllele
                        :tHap.TypedOnlyVariantList[MarkerIndex].refAllele)]);

            ifprintf(hapdose, "\t%.5f",outAllele1);


        }
    }

    ifprintf(hapdose, "\n");
    ifprintf(haps, "\n");

}


void Imputation::PrintInfoFile(HaplotypeSet &rHap,HaplotypeSet &tHap,  ImputationStatistics &stats)

{
    cout<<endl<<" Writing summary (.info) files ... "<<endl;
    IFILE info = ifopen(outFile + ".info", "wb");
    ifprintf(info, "SNP\tREF(0)\tALT(1)\tALT_Frq\tMAF\tAvgCall\tRsq\tGenotyped\tLooRsq\tEmpR\tEmpRsq\tDose0\tDose1\n");


    int i=0;
    for (int index =0; index < rHap.RefTypedTotalCount; index++)
    {

        if(rHap.RefTypedIndex[index]==-1)
        {

            if(i>=rHap.PrintStartIndex && i <= rHap.PrintEndIndex)
            {
                ifprintf(info, "%s\t%s\t%s\t%.5f\t%.5f\t%.5f\t%.5f\t",
                RsId? rHap.VariantList[i].rsid.c_str(): rHap.VariantList[i].name.c_str(),
                rHap.VariantList[i].refAlleleString.c_str(),
                rHap.VariantList[i].altAlleleString.c_str(),
                stats.AlleleFrequency(i),
                stats.AlleleFrequency(i) > 0.5 ? 1.0 - stats.AlleleFrequency(i) : stats.AlleleFrequency(i),
                stats.AverageCallScore(i),
                stats.Rsq(i));

                if (!tHap.missing[i])
                {
                    ifprintf(info, "Genotyped\t%.3f\t%.3f\t%.5f\t%.5f\t%.5f\n",
                      stats.LooRsq(i), stats.EmpiricalR(i), stats.EmpiricalRsq(i),
                      stats.LooMajorDose(i), stats.LooMinorDose(i));
                }
                else
                 ifprintf(info, "Imputed\t-\t-\t-\t-\t-\n");
            }
            i++;
        }
        else
        {
            variant ThisTypedVariant =tHap.TypedOnlyVariantList[rHap.RefTypedIndex[index]];

            ifprintf(info, "%s\t%s\t%s\t%.5f\t%.5f\t-\t-\tTyped_Only\t-\t-\t-\t-\t-\n",
            RsId? ThisTypedVariant.rsid.c_str(): ThisTypedVariant.name.c_str(),
            ThisTypedVariant.refAlleleString.c_str(),
            ThisTypedVariant.altAlleleString.c_str(),
            tHap.AlleleFreq[rHap.RefTypedIndex[index]],
            tHap.AlleleFreq[rHap.RefTypedIndex[index]] > 0.5 ?
                        1.0 - tHap.AlleleFreq[rHap.RefTypedIndex[index]] : tHap.AlleleFreq[rHap.RefTypedIndex[index]]);

        }
    }
    ifclose(info);


    cout<<endl<<" Summary information written to          : "<<outFile<<".info"<<endl;
   }


void Imputation::splitFoldedProb(vector<float> &SplitProb, vector<float> &totalProb, vector<float> &noRecomProb)
{
    SplitProb.resize(totalProb.size());

    for(int i=0;i<(int)totalProb.size();i++)
    {
        SplitProb[i]=totalProb[i]-noRecomProb[i];
    }
}




void Imputation::Condition(HaplotypeSet &rHap, int markerPos,vector<float> &Prob,
                            vector<float> &noRecomProb,
                            double e, double freq, bool observed, double backgroundError, int NoRedStates, ReducedHaplotypeInfo &Info)
{

    double prandom = e*freq+backgroundError;
    double pmatch = (1.0 - e)+e*freq+backgroundError;

    for (int i = 0; i<NoRedStates; i++)
    {

        bool allele=Info.returnHapAtPosition(i,markerPos);
        if(allele==observed)
        {
            Prob[i]*=pmatch;
            noRecomProb[i]*=pmatch;
        }
        else
        {
            Prob[i]*=prandom;
            noRecomProb[i]*=prandom;
        }
    }
}



void Imputation::ConditionJunctionProb(HaplotypeSet &rHap, int markerPos,vector<float> &Prob,
                            double e, double freq, bool observed, double backgroundError,
                            ReducedHaplotypeInfo &Info)
{

    double prandom = e*freq+backgroundError;
    double pmatch = (1.0 - e)+e*freq+backgroundError;
    int NoStates=rHap.numHaplotypes;

    for (int i = 0; i<NoStates; i++)
    {
        Prob[i]*=Info.returnHapAtPosition(Info.uniqueIndexMap[i],markerPos)==observed ?
                 pmatch:prandom;
    }
}




void Imputation::LooOptimalStructure(int loo,HaplotypeSet &rHap, HaplotypeSet &HapLoo)
{

    HapLoo.numHaplotypes=rHap.numHaplotypes-1;
    HapLoo.numMarkers=rHap.numMarkers;
    HapLoo.maxBlockSize=rHap.maxBlockSize;
    HapLoo.NoBlocks=rHap.NoBlocks;
    HapLoo.maxRepSize=rHap.maxRepSize;
    HapLoo.ReducedStructureInfo.clear();
    HapLoo.ReducedStructureInfo.resize(rHap.NoBlocks);

    for(int i=0;i<rHap.NoBlocks;i++)
    {

        ReducedHaplotypeInfo &ToBlock=HapLoo.ReducedStructureInfo[i];
        ReducedHaplotypeInfo &FromBlock=rHap.ReducedStructureInfo[i];


        int looIndex=FromBlock.uniqueIndexMap[loo];

        ToBlock.uniqueCardinality=FromBlock.uniqueCardinality;
        ToBlock.startIndex=FromBlock.startIndex;
        ToBlock.endIndex=FromBlock.endIndex;
        ToBlock.BlockSize=FromBlock.BlockSize;

        ToBlock.uniqueCardinality[looIndex]--;
        if(ToBlock.uniqueCardinality[looIndex]==0)
        {

            ToBlock.uniqueHaps.resize(FromBlock.RepSize-1);
            ToBlock.uniqueCardinality.clear();
            ToBlock.uniqueCardinality.resize(FromBlock.RepSize-1);
            for (int in = 0, out = 0; in < FromBlock.RepSize; in++)
            {
                if (in != looIndex)
                {
                    ToBlock.uniqueHaps[out] = FromBlock.uniqueHaps[in];
                    ToBlock.uniqueCardinality[out++] = FromBlock.uniqueCardinality[in];
                }
            }

            ToBlock.uniqueIndexMap.resize(rHap.numHaplotypes-1);
            for (int in = 0, out = 0; in <(rHap.numHaplotypes); in++)
            {

                if(FromBlock.uniqueIndexMap[in]<looIndex)
                    ToBlock.uniqueIndexMap[out++] = FromBlock.uniqueIndexMap[in];
                else if(FromBlock.uniqueIndexMap[in]>looIndex)
                    ToBlock.uniqueIndexMap[out++] = FromBlock.uniqueIndexMap[in]-1;
            }
        }
        else
        {
            ToBlock.uniqueHaps=FromBlock.uniqueHaps;
            ToBlock.uniqueIndexMap.resize(rHap.numHaplotypes-1);
            for (int in = 0, out = 0; in <(rHap.numHaplotypes); in++)
            {
                if(in!=loo)
                    ToBlock.uniqueIndexMap[out++] = FromBlock.uniqueIndexMap[in];

            }

        }

        ToBlock.RepSize=ToBlock.uniqueCardinality.size();
        ToBlock.InvuniqueCardinality.resize(ToBlock.RepSize);

        for (int in = 0; in < ToBlock.RepSize; in++)
        {
            ToBlock.InvuniqueCardinality[in]=1.0/(float)ToBlock.uniqueCardinality[in];
        }

    }

}
