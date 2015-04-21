#ifndef __IMPUTATIONSTATISTICS_H__
#define __IMPUTATIONSTATISTICS_H__

#include "MathVector.h"
#include "IntArray.h"

class ImputationStatistics
   {
   public:
      ImputationStatistics(int markers);
      ~ImputationStatistics();

      void Update(vector<float> &doses, vector<float> &loo,vector<bool> &observed,vector<bool> &Miss, vector<bool> &major);
//      void GoldenUpdate(vector<float> &doses, vector<float> &loo,vector<char> &observed, vector<char> &major);
//      void GoldenUpdate2(vector<float> &doses, vector<char> &observed1,vector<char> &observed2, vector<char> &major);

      double Rsq(int marker);
      double AlleleFrequency(int marker);
      double AverageCallScore(int marker);
      double LooRsq(int marker);
      double LooMajorDose(int marker);
      double LooMinorDose(int marker);
      double EmpiricalR(int marker);
      double EmpiricalRsq(int marker);
      double GoldenR(int marker);
      double GoldenRsq(int marker);
      double GoldenRgeno(int marker);
      double GoldenRsqgeno(int marker);

  // private:
      vector<double>   TargetPanelFreq, sum, sumSq, sumCall, looSum, looSumSq, looProduct, looObserved, looObservedSq;
      vector<int> count, looCount;
   };

#endif
