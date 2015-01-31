#ifndef __MARKOVMODEL_H__
#define __MARKOVMODEL_H__


#include "MarkovParameters.h"
#include "HaplotypeSet.h"
#include "StringBasics.h"
#include "MathVector.h"
#include "MarkovModel.h"

#include "ImputationStatistics.h"
#include<boost/random/uniform_real.hpp>
#include<boost/random/variate_generator.hpp>
#include<boost/random/mersenne_twister.hpp>



class MarkovModel : public MarkovParameters
{
    public:


        double backgroundError;
        int refCount,tarCount,noReducedStatesCurrent;
        vector<vector<vector<double> > > leftProb;
        vector<vector<vector<double> > > leftNoRecoProb;
        vector<vector<double> > junctionLeftProb;
        vector<vector<double> > junctionRightProb;
        int NoPrecisionJumps;
        vector<char> major;
        vector<bool> missing;
        vector<bool> PrecisionJump;
        vector<double>    imputedDose, imputedHap, leaveOneOut;
        vector<char> imputedAlleles,imputedAlleleNumber;


        bool            Transpose       (vector<double> &from,vector<double> &to,  vector<double> &noRecomProb, double reco,vector<int> &uniqueCardinality);
        void            Condition       (int markerPos,vector<double> &Prob, vector<double> &noRecomProb,
                                         char observed,double e,double freq,ReducedHaplotypeInfo &Info);
        void            WalkLeft        (HaplotypeSet &tHap, int &hapID,
                                        vector<vector<double> > &leftProb,vector<vector<double> > &noRecomProb,vector<double> &foldedProb,
                                        int start,int end,ReducedHaplotypeInfo &Info,vector<vector<double> > &alleleFreq);

        void            Impute          (HaplotypeSet &tHap, vector<double> &rightFoldProb,
                                         int &hapID,
                                         vector<vector<double> > &leftProb, vector<vector<double> > &leftNoRecomProb,
                                         vector<double> &rightProb, vector<double> &rightNoRecomProb,
                                         vector<double> &juncLeftProb,vector<double> &juncRightProb,
                                         int start,int end,ReducedHaplotypeInfo &Info,int type,vector<vector<double> > &alleleFreq);

        void            Impute          (int position, char observed,
                                        vector<double> &leftProb,vector<double> &rightProb,
                                        vector<double> &leftNoRecoProb,vector<double> &rightNoRecoProb,
                                        vector<double> &leftEndProb,vector<double> &rightEndProb,
                                        vector<double> &Constants,ReducedHaplotypeInfo &Info,
                                        vector<vector<double> > &alleleFreq);
        void            initializeMatrices
                                        (HaplotypeSet &tHap,HaplotypeSet &rHap,
                                         vector<int> &optStructure,
                                         vector<ReducedHaplotypeInfo> &StructureInfo);


         void            ProfileModel    (HaplotypeSet &tHap, HaplotypeSet &rHap,
                                        vector<vector<double> > &leftProb, vector<vector<double> > &leftNoRecomProb,
                                        vector<double> &juncLeftProbLeft,vector<double> &juncLeftProbRight,
                                        int start,int end,ReducedHaplotypeInfo &Info,int &currentState,boost::variate_generator<boost::mt19937&, boost::uniform_real<> > &uni);

        double          CountErrors     (vector<double> &probHap,
                                        int position, char observed, double e,double freq, ReducedHaplotypeInfo &Info);


        double          CountRecombinants(vector<double> &from, vector<double> &to,
                                        vector<double> &probHap,
                                        double r,bool PrecisionMultiply);

        vector<double>                  foldProbabilities           (int bridgeIndex,ReducedHaplotypeInfo &Info,int direction,int noReference);
        void                            unfoldProbabilities         (int bridgeIndex,vector<double> &recomProb, vector<double> &noRecomProb,
                                                                    vector<double> &foldedProb,int direction,vector<ReducedHaplotypeInfo> &StructureInfo,int noReference);

        void            CountExpected   (HaplotypeSet &tHap,int hapID,vector<double> &rightFoldProb,
                                         vector<vector<double> > &leftProb, vector<vector<double> > &leftNoRecomProb,
                                         vector<double> &rightProb, vector<double> &rightNoRecomProb,
                                         vector<double> &juncLeftProb,vector<double> &juncRightProb,
                                         int start,int end,ReducedHaplotypeInfo &Info,vector<vector<double> > &alleleFreq);

        void            CreatePosteriorProb(vector<double> &leftProb,vector<double> &rightProb,
                                         vector<double> &leftNoRecoProb,vector<double> &rightNoRecoProb,
                                         vector<double> &leftEndProb,vector<double> &rightEndProb,
                                         vector<double> &Constants,vector<double> &probHap,ReducedHaplotypeInfo &Info);

                        MarkovModel     (HaplotypeSet &tHap,HaplotypeSet &rHap, vector<bool> &miss, vector<double> &err, vector<double> &reco,vector<char> &maj)
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
                        MarkovModel     (HaplotypeSet &tHap,HaplotypeSet &rHap, vector<bool> &miss,vector<char> &maj)
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
