#ifndef IMPUTATION_H_INCLUDED
#define IMPUTATION_H_INCLUDED
#include "MarkovParameters.h"

#include "HaplotypeSet.h"
#include "MemoryInfo.h"
#include "MemoryAllocators.h"
#include "MarkovModel.h"
#include "MarkovParameters.h"
#include <ctime>
#include "Unique.h"
#include <cstdio>
#include <cmath>


#ifdef _OPENMP
#include <omp.h>
#endif
using namespace std;




class Imputation
    {
        public:

            int targetCount,refCount,markerCount,noReducedStatesCurrent;
            String outFile;
            int MaxSample;
            bool full,makeOptmxOnly;
            bool em,vcfOutput,doseOutput,onlyRefMarkers;
            bool phased;
            bool gzip;
            bool GT,DS,GP;
            vector<bool> format;
            bool includeGwas,updateModel;
            String errFile,recFile;
            int EstimationRounds,EstimationStates;

        Imputation(HaplotypeSet &rHap,String &out,bool gz)
        {
            refCount=rHap.numHaplotypes;
            markerCount=rHap.numMarkers;
            outFile=out;
            gzip=gz;
            updateModel=true;
//            MaxSample=20000000;
//            MaxSample=(MaxSample>(int)targetCount?(int)targetCount:MaxSample);

        }

        Imputation(HaplotypeSet &tHap,HaplotypeSet &rHap,String &out,String &err,
                   String &rec,bool Ph,bool gz,int rounds, int states,bool vcfoutput,
                   bool doseoutput, bool onlyrefmarkers,vector<bool> &Format )
        {

            targetCount=tHap.numHaplotypes;
            refCount=rHap.numHaplotypes;
            markerCount=rHap.numMarkers;
            outFile=out;
            MaxSample=2000000;
            MaxSample=(MaxSample>(int)targetCount?(int)targetCount:MaxSample);
            errFile=err;
            recFile=rec;
            EstimationRounds=rounds;
            EstimationStates=states;
            em=true;
            phased=Ph;
            gzip=gz;
            vcfOutput=vcfoutput;
            doseOutput=doseoutput;
            onlyRefMarkers=onlyrefmarkers;
            format=Format;
            GT=Format[0];
            DS=Format[1];
            GP=Format[2];

            updateModel=false;


        };
       Imputation(HaplotypeSet &tHap,HaplotypeSet &rHap,String &out,String &err,
                   String &rec,bool Ph,bool gz,int rounds, int states,bool vcfoutput,
                   bool doseoutput, bool onlyrefmarkers,vector<bool> &Format, bool updateMODEL)
        {

            updateModel=updateMODEL;
            targetCount=tHap.numHaplotypes;
            refCount=rHap.numHaplotypes;
            markerCount=rHap.numMarkers;
            outFile=out;
            MaxSample=2000000;
            MaxSample=(MaxSample>(int)targetCount?(int)targetCount:MaxSample);
            errFile=err;
            recFile=rec;
            EstimationRounds=rounds;
            EstimationStates=states;
            em=true;
            phased=Ph;
            gzip=gz;
            vcfOutput=vcfoutput;
            doseOutput=doseoutput;
            onlyRefMarkers=onlyrefmarkers;
            format=Format;
            GT=Format[0];
            DS=Format[1];
            GP=Format[2];


        };
        void                            createImputationStatistics  (HaplotypeSet &rHap,HaplotypeSet &tHap);
        void                            CreateOptmFileOnly          (HaplotypeSet &rHap);

        MarkovParameters*               createEstimates             (HaplotypeSet &rHap,HaplotypeSet &tHap,vector<int> &optStructure,bool NoTargetEstimation);
        void                            splitFoldedProb             (vector<float> &SplitProb,vector<float> &totalProb, vector<float> &noRecomProb);
        void                            normalize                   (vector<float> &x);
        void                            performImputation           (HaplotypeSet &tHap,HaplotypeSet &rHap, String Golden);
        void                            initializeMatrices          (HaplotypeSet &tHap,HaplotypeSet &rHap, vector<int> &optStructure,vector<ReducedHaplotypeInfo> &StructureInfo);
        void                            initializeMatrices1         (HaplotypeSet &tHap,HaplotypeSet &rHap, vector<int> &optStructure);
        void                            calculateFreq               (HaplotypeSet &rHap);
        void                            LooOptimalStructure         ( vector<ReducedHaplotypeInfo> &StructureInfo_loo,int loo,HaplotypeSet &rHap);
        void                            Condition                   (HaplotypeSet &rHap, int markerPos,vector<float> &Prob,
                                                                    vector<float> &noRecomProb,
                                                                    double e, double freq, char observed, double backgroundError, int NoRedStates, ReducedHaplotypeInfo &Info);
        void                            FlushPartialVcf             (HaplotypeSet &rHap,HaplotypeSet &tHap,HaplotypeSet &PartialDosage, string &filename,int &Index);
        void                            MergeFinalVcf               (HaplotypeSet &rHap,HaplotypeSet &tHap,ImputationStatistics &stats,int MaxIndex);
};




#endif // IMPUTATION_H_INCLUDED
