#include "Imputation.h"

MarkovParameters* Imputation::createEstimates(HaplotypeSet &rHap,HaplotypeSet &tHap,vector<int> &optStructure,bool NoTargetEstimation)
{
    MarkovParameters *MP=new MarkovParameters(rHap.numMarkers);
    rHap.CalculateFreq();

    cout<<" ------------------------------------------------------------------------------"<<endl;
    cout<<"                             PARAMETER ESTIMATION                              "<<endl;
    cout<<" ------------------------------------------------------------------------------"<<endl;

    int time_prev = time(0);

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
            rHap_loo.numHaplotypes=rHap.numHaplotypes-1;
            rHap_loo.numMarkers=rHap.numMarkers;
            vector<float> foldedProb,recomProb,noRecomProb,
                            rightProbTemp,probAlleleNoStandardize(8,0.0);

           // Target with one panel at position i
            HaplotypeSet tHap_loo;

            vector<bool> tempHap(rHap.numMarkers);
            rHap.reconstructHaplotype(tempHap,i);
            tHap_loo.Create(tempHap);

            // New Strcture with State Space Reduction
            vector<ReducedHaplotypeInfo> StructureInfo_loo;
            LooOptimalStructure(StructureInfo_loo,i,rHap);

            // Create Markov Model

            vector<bool> missing_loo(markerCount,false);

            MarkovModel MM(tHap_loo,rHap_loo,missing_loo,rHap.major);

            MM.CopyParameters(MP);

            MM.initializeMatrices(tHap_loo,rHap_loo,optStructure,StructureInfo_loo);

            for(int group=1;group<(int)optStructure.size();group++)
            {

                MM.foldProbabilities(foldedProb,group-1,StructureInfo_loo[group-1],0,refCount-1);
                MM.leftNoRecoProb[group-1][0]=foldedProb;


                if(group==1 && !missing_loo[0])
                    {
                        if(!tHap_loo.getMissingScaffoldedHaplotype(0,0))
                            {

                                Condition(rHap,0,foldedProb,MM.leftNoRecoProb[group-1][0],MM.Error[0],
                                        tHap_loo.getScaffoldedHaplotype(0,0)? rHap.AlleleFreq[0] : 1-rHap.AlleleFreq[0],
                                        tHap_loo.getScaffoldedHaplotype(0,0),MM.backgroundError, foldedProb.size(),StructureInfo_loo[0]);
                            }
                    }


                int hapID=0;
                MM.WalkLeft(tHap_loo,hapID,MM.leftProb[group-1],MM.leftNoRecoProb[group-1],
                            foldedProb,optStructure[group-1],optStructure[group],
                            StructureInfo_loo[group-1],rHap.AlleleFreq);
                splitFoldedProb(recomProb,MM.leftProb[group-1][optStructure[group]-optStructure[group-1]],MM.leftNoRecoProb[group-1][optStructure[group]-optStructure[group-1]]);
                MM.unfoldProbabilities(group-1,recomProb,MM.leftNoRecoProb[group-1][optStructure[group]-optStructure[group-1]],foldedProb,0,StructureInfo_loo,refCount-1);

            }


            if(em)
            {
                for(int group=optStructure.size()-1;group>0;group--)
                {
                    MM.foldProbabilities(foldedProb,group-1,StructureInfo_loo[group-1],1,refCount-1);
                    noRecomProb=foldedProb;
                    MM.CountExpected(tHap_loo,0,foldedProb,MM.leftProb[group-1],MM.leftNoRecoProb[group-1],rightProbTemp,noRecomProb,
                                     MM.junctionLeftProb[group-1],MM.junctionRightProb[group], optStructure[group-1],
                                     optStructure[group],StructureInfo_loo[group-1],rHap.AlleleFreq);
                    splitFoldedProb(recomProb,rightProbTemp,noRecomProb);
                    MM.unfoldProbabilities(group-1,recomProb,noRecomProb,foldedProb,1,StructureInfo_loo,refCount-1);
                }

            MM.empiricalCount++;

            }
            #pragma omp critical
            {
                (*MP)+=MM;
                double tempLogLikelihoodValue=0.0;
                for(int tempVar=0;tempVar<(int)MM.leftProb[MM.leftProb.size()-1][MM.leftProb[MM.leftProb.size()-1].size()-1].size();tempVar++)
                {
                    tempLogLikelihoodValue+=MM.leftProb[MM.leftProb.size()-1][MM.leftProb[MM.leftProb.size()-1].size()-1][tempVar];
                }



                tempLogLikelihoodValue=log(tempLogLikelihoodValue);
                tempLogLikelihoodValue-=(MM.NoPrecisionJumps*log(1e15));
                LogLikelihoodValue+=tempLogLikelihoodValue;
            }

        }



        if(!NoTargetEstimation && round>=EstimationRounds/2)
        {
            iterations = EstimationStates < (int)tHap.numHaplotypes ? EstimationStates : (int)tHap.numHaplotypes;

            #pragma omp parallel for
            for (int i = 0; i < iterations; i++)
            {
                MarkovModel MM(tHap,rHap,tHap.missing,rHap.major);
                vector<float> foldedProb,recomProb,noRecomProb, rightProbTemp,probAlleleNoStandardize(8,0.0);


                MM.CopyParameters(MP);
                //initializeMatrices(tHap,rHap,optStructure,rHap.ReducedStructureInfo);
                MM.initializeMatrices(tHap,rHap,optStructure,rHap.ReducedStructureInfo);
                for(int group=1;group<(int)optStructure.size();group++)
                {
                    MM.foldProbabilities(foldedProb,group-1,rHap.ReducedStructureInfo[group-1],0,refCount);
                    MM.leftNoRecoProb[group-1][0]=foldedProb;

                    if(group==1 && !tHap.missing[0])
                        if(!tHap.getMissingScaffoldedHaplotype(i,0))
                            Condition(rHap,0,foldedProb,MM.leftNoRecoProb[group-1][0],MM.Error[0],
                                tHap.getScaffoldedHaplotype(i,0)? rHap.AlleleFreq[0] : 1-rHap.AlleleFreq[0],
                                tHap.getScaffoldedHaplotype(i,0),MM.backgroundError, foldedProb.size(),rHap.ReducedStructureInfo[0]);



                    MM.WalkLeft(tHap,i,MM.leftProb[group-1],MM.leftNoRecoProb[group-1],foldedProb,optStructure[group-1],optStructure[group],rHap.ReducedStructureInfo[group-1],rHap.AlleleFreq);

                    splitFoldedProb(recomProb,MM.leftProb[group-1][optStructure[group]-optStructure[group-1]],MM.leftNoRecoProb[group-1][optStructure[group]-optStructure[group-1]]);

                    MM.unfoldProbabilities(group-1,recomProb,MM.leftNoRecoProb[group-1][optStructure[group]-optStructure[group-1]],foldedProb,0,rHap.ReducedStructureInfo,refCount);


                }

                if(em)
                {
                    for(int group=optStructure.size()-1;group>0;group--)
                    {
                        MM.foldProbabilities(foldedProb,group-1,rHap.ReducedStructureInfo[group-1],1,refCount);
                        noRecomProb=foldedProb;



                        MM.CountExpected(tHap,i,foldedProb,MM.leftProb[group-1],MM.leftNoRecoProb[group-1],rightProbTemp,noRecomProb,
                                         MM.junctionLeftProb[group-1],MM.junctionRightProb[group],optStructure[group-1],optStructure[group],rHap.ReducedStructureInfo[group-1],rHap.AlleleFreq);



                        splitFoldedProb(recomProb,rightProbTemp,noRecomProb);
                        MM.unfoldProbabilities(group-1,recomProb,noRecomProb,foldedProb,1,rHap.ReducedStructureInfo,refCount);
                    }



                    MM.empiricalCount++;
                }


                #pragma omp critical
                 {
                    (*MP)+=MM;
                    double tempLogLikelihoodValue=0.0;
                    for(int tempVar=0;tempVar<(int)MM.leftProb[MM.leftProb.size()-1][MM.leftProb[MM.leftProb.size()-1].size()-1].size();tempVar++)
                    {
                        tempLogLikelihoodValue+=MM.leftProb[MM.leftProb.size()-1][MM.leftProb[MM.leftProb.size()-1].size()-1][tempVar];
                    }
                    tempLogLikelihoodValue=log(tempLogLikelihoodValue);
                    tempLogLikelihoodValue-=(MM.NoPrecisionJumps*log(1e15));
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

    if(rHap.vcfType)
        rHap.writem3vcfFile(outFile,gzip);

    cout << "\n Time taken for parameter estimation = " << time_load << " seconds. "<<endl<<endl;

    return MP;


}


void Imputation::MergeFinalVcf(HaplotypeSet &rHap,HaplotypeSet &tHap,ImputationStatistics &stats,int MaxIndex)
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

    for (int i = rHap.PrintStartIndex; i <= rHap.PrintEndIndex; i++)
   {


        if(i%10000==0)
            {
                printf("    Merging marker %d of %d [%.1f%%] to VCF File ...", i + 1, rHap.numMarkers,100*(double)(i + 1)/(int)rHap.numMarkers);
                cout<<endl;
            }


        ifprintf(vcfdosepartial,"\n%s\t%d\t%s\t%s\t%s\t.\tPASS\tMAF=%.5f;R2=%.5f",
        rHap.VariantList[i].chr.c_str(),rHap.VariantList[i].bp,
        rHap.VariantList[i].name.c_str(),rHap.VariantList[i].refAlleleString.c_str(),
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

    for (int i = rHap.PrintStartIndex; i <= rHap.PrintEndIndex; i++)
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


    ifclose(vcfdosepartial);

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
//    HaplotypeSet DosageForVcf;
    HaplotypeSet DosageForVcfPartial;
    DosageForVcfPartial.unphasedOutput=unphasedOutput;


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
            DosageForVcfPartial.InitializePartialDosageForVcfOutput((tHap.AllMaleTarget?maxVcfSample:2*maxVcfSample),rHap.numMarkers,format);
        else
            DosageForVcfPartial.InitializePartialDosageForVcfOutputMaleSamples((tHap.AllMaleTarget?maxVcfSample:2*maxVcfSample),rHap.numMarkers,format);
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

        vector<float> foldedProb,recomProb,noRecomProb, rightProb,probAlleleNoStandardize(8,0.0),tempDoseHap1;
        vector<bool> tempHap(rHap.numMarkers),tempMissHap(rHap.numMarkers);
        vector<bool> tempDoseAlleleHap1;

        MarkovModel MM(tHap,rHap,tHap.missing,rHap.major);

        MM.CopyParameters(MP);

        int hapIdIndiv=hapId;

        do{

            MM.initializeMatrices(tHap,rHap,optStructure,rHap.ReducedStructureInfo);
            printf("  Processing Haplotype %d of %d ...", hapIdIndiv + 1, MaxSample);
            cout<<endl;


            MM.ThisHapId=hapIdIndiv;


            for(int group=1;group<(int)optStructure.size();group++)
            {

                MM.foldProbabilities(foldedProb,group-1,rHap.ReducedStructureInfo[group-1],0,refCount);
                MM.leftNoRecoProb[group-1][0]=foldedProb;


                if(group==1 && !tHap.missing[0])
                        if(!tHap.getMissingScaffoldedHaplotype(hapIdIndiv,0))
                            {

                                Condition(rHap,0,foldedProb,MM.leftNoRecoProb[group-1][0],MM.Error[0],
                                tHap.getScaffoldedHaplotype(hapIdIndiv,0)? rHap.AlleleFreq[0] : 1-rHap.AlleleFreq[0],
                                tHap.getScaffoldedHaplotype(hapIdIndiv,0),MM.backgroundError,
                                      foldedProb.size(),rHap.ReducedStructureInfo[0]);
                            }



                MM.WalkLeft(tHap,hapIdIndiv,MM.leftProb[group-1],MM.leftNoRecoProb[group-1],
                            foldedProb,optStructure[group-1],optStructure[group],
                            rHap.ReducedStructureInfo[group-1],rHap.AlleleFreq);

                splitFoldedProb(recomProb,MM.leftProb[group-1][optStructure[group]-optStructure[group-1]],MM.leftNoRecoProb[group-1][optStructure[group]-optStructure[group-1]]);

                MM.unfoldProbabilities(group-1,recomProb,MM.leftNoRecoProb[group-1][optStructure[group]-optStructure[group-1]],foldedProb,0,rHap.ReducedStructureInfo,refCount);

            }



            for(int group=optStructure.size()-1;group>0;group--)
            {

                MM.foldProbabilities(foldedProb,group-1,rHap.ReducedStructureInfo[group-1],1,refCount);
                rightProb=foldedProb;
                noRecomProb=foldedProb;

                MM.Impute(tHap,foldedProb,hapIdIndiv,MM.leftProb[group-1],MM.leftNoRecoProb[group-1],rightProb,noRecomProb,MM.junctionLeftProb[group-1],
                          MM.junctionRightProb[group],optStructure[group-1], optStructure[group],rHap.ReducedStructureInfo[group-1],1,rHap.AlleleFreq);

                splitFoldedProb(recomProb,rightProb,noRecomProb);
                MM.unfoldProbabilities(group-1,recomProb,noRecomProb,foldedProb,1,rHap.ReducedStructureInfo,refCount);
            }

            for(int jjj=0;jjj<rHap.numMarkers;jjj++)
                {
                    tempHap[jjj]=tHap.getScaffoldedHaplotype(hapIdIndiv,jjj);
                    tempMissHap[jjj]=tHap.getMissingScaffoldedHaplotype(hapIdIndiv,jjj);

                }

            if(vcfOutput)
            {
                if(hapIdIndiv%2==0)
                {
                   tempDoseHap1= MM.imputedHap;
                   tempDoseAlleleHap1= MM.imputedAlleleNumber;
                }
            }
            #pragma omp critical
            {
                stats.Update(MM.imputedHap, MM.leaveOneOut,tempHap,tempMissHap,rHap.major);
            }

            #pragma omp critical
            if (phased && !unphasedOutput)
            {

                printf("    Outputting HAPLO%d of Individual %s for Haplotype File...", tHap.AllMaleTarget?1:(hapIdIndiv%2+1) ,tHap.individualName[tHap.AllMaleTarget?hapId:hapId/2].c_str());

                cout<<endl;
                ifprintf(hapdose, "%s\tHAPLO%d",  tHap.individualName[tHap.AllMaleTarget?hapId:hapId/2].c_str(), tHap.AllMaleTarget?1:(hapIdIndiv%2+1) );
                ifprintf(haps, "%s\tHAPLO%d\t", tHap.individualName[tHap.AllMaleTarget?hapId:hapId/2].c_str(), tHap.AllMaleTarget?1:(hapIdIndiv%2+1) );
                for (int j = rHap.PrintStartIndex; j <= rHap.PrintEndIndex; j++)
                {
                    ifprintf(hapdose, "\t%.5f", MM.imputedHap[j]);
                    ifprintf(haps, "%c", MM.imputedAlleles[j]?
                                    rHap.VariantList[j].refAllele
                                        :rHap.VariantList[j].altAllele);
                }

            ifprintf(hapdose, "\n");
            ifprintf(haps, "\n");
            }


            if(tHap.AllMaleTarget)
                break;
            hapIdIndiv++;
        }while(hapIdIndiv<MaxSample && hapIdIndiv%2==1);

        #pragma omp critical
        if(doseOutput)
        {
            printf("    Outputting Individual %s for Dosage file...",  tHap.individualName[tHap.AllMaleTarget?hapId:hapId/2].c_str());
            cout<<endl;
            ifprintf(dosages, "%s\tDOSE",tHap.individualName[tHap.AllMaleTarget?hapId:hapId/2].c_str());
            for (int j = rHap.PrintStartIndex; j <= rHap.PrintEndIndex; j++)
                ifprintf(dosages, "\t%.3f", MM.imputedDose[j]);
            ifprintf(dosages, "\n");
        }
         #pragma omp critical
        if(vcfOutput)
        {

            printf("    Saving Individual %s for VCF File...\n",  tHap.individualName[tHap.AllMaleTarget?hapId:hapId/2].c_str());
            if(!tHap.AllMaleTarget)
                DosageForVcfPartial.SaveDosageForVcfOutputSampleWise(NumVcfCreated-NumVcfWritten,
                                                                 tHap.individualName[tHap.AllMaleTarget?hapId:hapId/2],
                                                                 tempDoseHap1,MM.imputedHap,
                                                                 tempDoseAlleleHap1,MM.imputedAlleleNumber);
            else
                DosageForVcfPartial.SaveDosageForVcfOutputSampleWiseChrX(NumVcfCreated-NumVcfWritten,
                                                                 tHap.individualName[tHap.AllMaleTarget?hapId:hapId/2],
                                                                 MM.imputedHap,
                                                                 MM.imputedAlleleNumber);
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
                       (NumVcfWritten)+1,(MaxSample/2<(NumVcfWritten+maxVcfSample)?MaxSample:(NumVcfWritten+maxVcfSample)),
                       PartialVcfFileName.c_str());



                FlushPartialVcf(rHap,tHap,DosageForVcfPartial,PartialVcfFileName,NovcfParts);
                if(NumVcfCreated<(tHap.AllMaleTarget?MaxSample:MaxSample/2))
                {
                    NovcfParts++;
                    NumVcfWritten+=maxVcfSample;
                    DosageForVcfPartial.InitializePartialDosageForVcfOutput(maxVcfSample<(((tHap.AllMaleTarget?MaxSample:MaxSample/2))-NumVcfWritten)?2*maxVcfSample:2*(((tHap.AllMaleTarget?MaxSample:MaxSample/2))-NumVcfWritten),rHap.numMarkers,format);

    if(!tHap.AllMaleTarget)
        DosageForVcfPartial.InitializePartialDosageForVcfOutput(maxVcfSample<(((tHap.AllMaleTarget?MaxSample:MaxSample/2))-NumVcfWritten)?2*maxVcfSample:2*(((tHap.AllMaleTarget?MaxSample:MaxSample/2))-NumVcfWritten),rHap.numMarkers,format);
    else
        DosageForVcfPartial.InitializePartialDosageForVcfOutputMaleSamples(maxVcfSample<(((tHap.AllMaleTarget?MaxSample:MaxSample/2))-NumVcfWritten)?2*maxVcfSample:2*(((tHap.AllMaleTarget?MaxSample:MaxSample/2))-NumVcfWritten),rHap.numMarkers,format);


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

    cout<<endl<<" Writing summary (.info) files ... "<<endl;

    IFILE info = ifopen(outFile + ".info", "wb");

    ifprintf(info, "SNP\tREF\tALT\tMajor\tMinor\tMAF\tAvgCall\tRsq\tGenotyped\tLooRsq\tEmpR\tEmpRsq\tDose1\tDose2\n");
    for (int i = rHap.PrintStartIndex; i <= rHap.PrintEndIndex; i++)
    {

        ifprintf(info, "%s\t%s\t%s\t%s\t%s\t%.5f\t%.5f\t%.5f\t",
            rHap.VariantList[i].name.c_str(),
            rHap.VariantList[i].refAlleleString.c_str(),
            rHap.VariantList[i].altAlleleString.c_str(),
            rHap.VariantList[i].MajAlleleString.c_str(),
            rHap.VariantList[i].MinAlleleString.c_str(),
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
         ifprintf(info, "-\t-\t-\t-\t-\t-\n");
    }
    ifclose(info);

    cout<<endl<<" Summary information written to          : "<<outFile<<".info"<<endl;
    time_load = time(0) - time_prev;
    cout << "\n Time taken for imputation = " << time_load << " seconds."<<endl<<endl;


    if(vcfOutput)
        MergeFinalVcf(rHap,tHap,stats,NovcfParts);

}




void Imputation::splitFoldedProb(vector<float> &SplitProb, vector<float> &totalProb, vector<float> &noRecomProb)
{
    SplitProb.resize(totalProb.size());

    for(int i=0;i<(int)totalProb.size();i++)
    {
        SplitProb[i]=totalProb[i]-noRecomProb[i];
    }
}




void Imputation::normalize(vector<float> &x)
{
    float sum=0.0;
    for(int i=0;i<(int)x.size();i++)
        sum+=x[i];

    for(int i=0;i<(int)x.size();i++)
        x[i]/=sum;

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



void Imputation::LooOptimalStructure( vector<ReducedHaplotypeInfo> &StructureInfo_loo,int loo,HaplotypeSet &rHap)
{
    StructureInfo_loo.clear();
    StructureInfo_loo.resize(rHap.ReducedStructureInfo.size());

    for(int i=0;i<(int)rHap.ReducedStructureInfo.size();i++)
    {


        int looIndex=rHap.ReducedStructureInfo[i].uniqueIndexMap[loo];

        StructureInfo_loo[i].uniqueCardinality=rHap.ReducedStructureInfo[i].uniqueCardinality;
        StructureInfo_loo[i].InvuniqueCardinality=rHap.ReducedStructureInfo[i].InvuniqueCardinality;
        StructureInfo_loo[i].startIndex=rHap.ReducedStructureInfo[i].startIndex;
        StructureInfo_loo[i].endIndex=rHap.ReducedStructureInfo[i].endIndex;



        StructureInfo_loo[i].uniqueCardinality[looIndex]--;
        if(StructureInfo_loo[i].uniqueCardinality[looIndex]==0)
        {

            StructureInfo_loo[i].uniqueHaps.resize(rHap.ReducedStructureInfo[i].uniqueHaps.size()-1);
            StructureInfo_loo[i].uniqueCardinality.clear();
            StructureInfo_loo[i].uniqueCardinality.resize(rHap.ReducedStructureInfo[i].uniqueCardinality.size()-1);
            for (int in = 0, out = 0; in <(int) rHap.ReducedStructureInfo[i].uniqueCardinality.size(); in++)
            {
                if (in != looIndex)
                {
                    StructureInfo_loo[i].uniqueHaps[out] = rHap.ReducedStructureInfo[i].uniqueHaps[in];
                    StructureInfo_loo[i].uniqueCardinality[out++] = rHap.ReducedStructureInfo[i].uniqueCardinality[in];
                }

            }


            StructureInfo_loo[i].uniqueIndexMap.resize(rHap.ReducedStructureInfo[i].uniqueIndexMap.size()-1);
            for (int in = 0, out = 0; in <(int) rHap.ReducedStructureInfo[i].uniqueIndexMap.size(); in++)
            {

                if(rHap.ReducedStructureInfo[i].uniqueIndexMap[in]<looIndex)
                    StructureInfo_loo[i].uniqueIndexMap[out++] = rHap.ReducedStructureInfo[i].uniqueIndexMap[in];
                else if(rHap.ReducedStructureInfo[i].uniqueIndexMap[in]>looIndex)
                    StructureInfo_loo[i].uniqueIndexMap[out++] = rHap.ReducedStructureInfo[i].uniqueIndexMap[in]-1;

            }
        }
        else
        {
            StructureInfo_loo[i].uniqueHaps=rHap.ReducedStructureInfo[i].uniqueHaps;
            StructureInfo_loo[i].uniqueIndexMap.resize(rHap.ReducedStructureInfo[i].uniqueIndexMap.size()-1);
            for (int in = 0, out = 0; in <(int) rHap.ReducedStructureInfo[i].uniqueIndexMap.size(); in++)
            {
                if(in!=loo)
                    StructureInfo_loo[i].uniqueIndexMap[out++] = rHap.ReducedStructureInfo[i].uniqueIndexMap[in];
            }

        }


        for (int in = 0; in < (int)StructureInfo_loo[i].uniqueCardinality.size(); in++)
        {
            StructureInfo_loo[i].InvuniqueCardinality[in]=1.0/(float)StructureInfo_loo[i].uniqueCardinality[in];
        }

    }

}
