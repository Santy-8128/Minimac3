#ifndef __MARKOVMODEL_H__
#define __MARKOVMODEL_H__


#include "MarkovParameters.h"
#include "HaplotypeSet.h"
#include "StringBasics.h"
#include "MathVector.h"
#include "MarkovModel.h"
#include "Random.h"

#include "Unique.h"

#include "ImputationStatistics.h"



class MarkovModel : public MarkovParameters
{
    public:

        int ThisHapId;
        bool LowMemory;
        double backgroundError;
        int refCount,tarCount,noReducedStatesCurrent;

        vector<vector<vector<float> > > leftProb;
        vector<vector<float> > ThisBlockLeftNoRecoProb ;
        vector<vector<float> > ThisBlockLeftProb ;
        vector<float> CurrentLeftNoRecoProb;


        vector<vector<float> > junctionLeftProb;
        vector<vector<float> > junctionLeftFoldProb;
        vector<float> PrevjunctionRightProb;


        vector<float> probHap;
        vector<float> Constants;
        vector<float> tempRightProb;

        int NoPrecisionJumps;
        vector<bool> major;
        vector<bool> missing;
        vector<bool> PrecisionJump;
        vector<float>    imputedDose, imputedHap, leaveOneOut;
        vector<bool> imputedAlleleNumber;






        bool            Transpose                       (vector<float> &from,vector<float> &to,  vector<float> &noRecomProb, double reco,vector<int> &uniqueCardinality);
        void            Condition                       (int markerPos,vector<float> &Prob, vector<float> &noRecomProb,
                                                        bool observed,double e,double freq,ReducedHaplotypeInfo &Info);
        void            WalkLeft                        (HaplotypeSet &tHap, int &hapID,
                                                        int group, ReducedHaplotypeInfo &Info,
                                                        vector<double> &alleleFreq);
        void            ReCreateBothLeftProb            (HaplotypeSet &tHap, int &hapID,
                                                        int group, ReducedHaplotypeInfo &Info,
                                                        vector<double> &alleleFreq);
        void            ReCreateLeftNoRecoProb          (HaplotypeSet &tHap, int &hapID,
                                                        int group, ReducedHaplotypeInfo &Info,
                                                        vector<double> &alleleFreq);

        void            Impute                          (HaplotypeSet &tHap,int hapID,int group,
                                                        vector<float> &PrevRightFoldedProb,
                                                        vector<float> &CurrentRightProb, vector<float> &CurrentNoRecoRightProb,
                                                        ReducedHaplotypeInfo &Info,vector<double> &alleleFreq);


        void            Impute                          (int position, bool observed, bool observedMiss,
                                                        vector<float> &leftProb,vector<float> &rightProb,
                                                        vector<float> &leftNoRecoProb,vector<float> &rightNoRecoProb,
                                                        vector<float> &leftEndProb,vector<float> &rightEndProb,
                                                        vector<float> &Constants,ReducedHaplotypeInfo &Info,
                                                        vector<double> &alleleFreq);
        void            initializeMatrices              (HaplotypeSet &rHap,HaplotypeSet &tHap);
        void            ReinitializeMatrices            ();


        double          CountErrors                     (vector<float> &probHap,
                                                        int position, bool observed, double e,double freq, ReducedHaplotypeInfo &Info);
        double          CountRecombinants               (vector<float> &from, vector<float> &to,
                                                        vector<float> &probHap,
                                                        double r,bool PrecisionMultiply);
        void            foldProbabilities               (vector<float> &foldProb,int bridgeIndex,ReducedHaplotypeInfo &Info,int direction,int noReference);
        void            unfoldProbabilities             (int bridgeIndex,vector<float> &recomProb, vector<float> &noRecomProb,
                                                        vector<float> &PrevFoldedProb,int direction,vector<ReducedHaplotypeInfo> &StructureInfo,int noReference);

        void            CountExpected                   (HaplotypeSet &tHap,int hapID,int group,
                                                        vector<float> &PrevRightFoldedProb,
                                                        vector<float> &CurrentRightProb, vector<float> &CurrentNoRecoRightProb,
                                                        ReducedHaplotypeInfo &Info,vector<double> &alleleFreq);

        void            CreatePosteriorProb             (vector<float> &leftProb,vector<float> &rightProb,
                                                        vector<float> &leftNoRecoProb,vector<float> &rightNoRecoProb,
                                                        vector<float> &leftEndProb,vector<float> &rightEndProb,
                                                        vector<float> &Constants,vector<float> &probHap,ReducedHaplotypeInfo &Info);

        void            CheckSize                       (HaplotypeSet & rHap,HaplotypeSet &tHap);


                        MarkovModel                     (HaplotypeSet &tHap,HaplotypeSet &rHap, vector<bool> &miss,vector<bool> &maj,bool lowM)
                                                        {
                                                            major=maj;
                                                            missing=miss;
                                                            PrecisionJump.resize(miss.size(),false);
                                                            refCount=rHap.numHaplotypes;
                                                            tarCount=tHap.numHaplotypes;
                                                            noMarker=rHap.numMarkers;
                                                            imputedDose.resize(noMarker);
                                                            imputedHap.resize(noMarker);
                                                            leaveOneOut.resize(noMarker);
                                                            imputedAlleleNumber.resize(noMarker);
                                                            NoPrecisionJumps=0;
                                                            LowMemory=lowM;
                                                            backgroundError = 1e-5;

                                                        }


    };

#endif
