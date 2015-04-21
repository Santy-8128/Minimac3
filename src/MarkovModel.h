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

        double backgroundError;
        int refCount,tarCount,noReducedStatesCurrent;
        vector<vector<vector<float> > > leftProb;
        vector<vector<vector<float> > > leftNoRecoProb;
        vector<vector<float> > junctionLeftProb;
        vector<vector<float> > junctionRightProb;
        int NoPrecisionJumps;
        vector<bool> major;
        vector<bool> missing;
        vector<bool> PrecisionJump;
        vector<float>    imputedDose, imputedHap, leaveOneOut;
        vector<bool> imputedAlleles,imputedAlleleNumber;


        bool            Transpose       (vector<float> &from,vector<float> &to,  vector<float> &noRecomProb, double reco,vector<int> &uniqueCardinality);
        void            Condition       (int markerPos,vector<float> &Prob, vector<float> &noRecomProb,
                                         bool observed,double e,double freq,ReducedHaplotypeInfo &Info);
        void            WalkLeft        (HaplotypeSet &tHap, int &hapID,
                                        vector<vector<float> > &leftProb,vector<vector<float> > &noRecomProb,vector<float> &foldedProb,
                                        int start,int end,ReducedHaplotypeInfo &Info,vector<double> &alleleFreq);

        void            Impute          (HaplotypeSet &tHap, vector<float> &rightFoldProb,
                                         int &hapID,
                                         vector<vector<float> > &leftProb, vector<vector<float> > &leftNoRecomProb,
                                         vector<float> &rightProb, vector<float> &rightNoRecomProb,
                                         vector<float> &juncLeftProb,vector<float> &juncRightProb,
                                         int start,int end,ReducedHaplotypeInfo &Info,int type,vector<double> &alleleFreq);

        void            Impute          (int position, bool observed, bool observedMiss,
                                        vector<float> &leftProb,vector<float> &rightProb,
                                        vector<float> &leftNoRecoProb,vector<float> &rightNoRecoProb,
                                        vector<float> &leftEndProb,vector<float> &rightEndProb,
                                        vector<float> &Constants,ReducedHaplotypeInfo &Info,
                                        vector<double> &alleleFreq);

        void            ImputeOld          (int position, bool observed, bool observedMiss,
                                        vector<float> &leftProb,vector<float> &rightProb,
                                        vector<float> &leftNoRecoProb,vector<float> &rightNoRecoProb,
                                        vector<float> &leftEndProb,vector<float> &rightEndProb,
                                        vector<float> &Constants,ReducedHaplotypeInfo &Info,
                                        vector<double> &alleleFreq);
        void            initializeMatrices
                                        (HaplotypeSet &tHap,HaplotypeSet &rHap,
                                         vector<int> &optStructure,
                                         vector<ReducedHaplotypeInfo> &StructureInfo);

        double          CountErrors     (vector<float> &probHap,
                                        int position, bool observed, double e,double freq, ReducedHaplotypeInfo &Info);


        double          CountRecombinants(vector<float> &from, vector<float> &to,
                                        vector<float> &probHap,
                                        double r,bool PrecisionMultiply);

        void            foldProbabilities           (vector<float> &foldProb,int bridgeIndex,ReducedHaplotypeInfo &Info,int direction,int noReference);
        void            unfoldProbabilities         (int bridgeIndex,vector<float> &recomProb, vector<float> &noRecomProb,
                                                                    vector<float> &foldedProb,int direction,vector<ReducedHaplotypeInfo> &StructureInfo,int noReference);

        void            unfoldProbabilitiesOld         (int bridgeIndex,vector<float> &recomProb, vector<float> &noRecomProb,
                                                                    vector<float> &foldedProb,int direction,vector<ReducedHaplotypeInfo> &StructureInfo,int noReference);

        void            CountExpected   (HaplotypeSet &tHap,int hapID,vector<float> &rightFoldProb,
                                         vector<vector<float> > &leftProb, vector<vector<float> > &leftNoRecomProb,
                                         vector<float> &rightProb, vector<float> &rightNoRecomProb,
                                         vector<float> &juncLeftProb,vector<float> &juncRightProb,
                                         int start,int end,ReducedHaplotypeInfo &Info,vector<double> &alleleFreq);

        void            CreatePosteriorProb(vector<float> &leftProb,vector<float> &rightProb,
                                         vector<float> &leftNoRecoProb,vector<float> &rightNoRecoProb,
                                         vector<float> &leftEndProb,vector<float> &rightEndProb,
                                         vector<float> &Constants,vector<float> &probHap,ReducedHaplotypeInfo &Info);

                        MarkovModel     (HaplotypeSet &tHap,HaplotypeSet &rHap, vector<bool> &miss, vector<double> &err, vector<double> &reco,vector<bool> &maj)
                        {
                            major=maj;
                            missing=miss;
                            PrecisionJump.resize(miss.size(),false);
                            Recom=reco;
                            Error=err;

                            noMarker=rHap.numMarkers;
                            imputedDose.resize(noMarker,0.0);
                            imputedHap.resize(noMarker);
                            leaveOneOut.resize(noMarker);
                            imputedAlleles.resize(noMarker);
                            imputedAlleleNumber.resize(noMarker);

                            refCount=rHap.numHaplotypes;
                            tarCount=tHap.numHaplotypes;
                            backgroundError = 1e-5;
                            NoPrecisionJumps=0;
                        }
                        MarkovModel     (HaplotypeSet &tHap,HaplotypeSet &rHap, vector<bool> &miss,vector<bool> &maj)
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
                            imputedAlleles.resize(noMarker);
                            imputedAlleleNumber.resize(noMarker);
                            NoPrecisionJumps=0;


                            backgroundError = 1e-5;


                        }





    };

#endif
