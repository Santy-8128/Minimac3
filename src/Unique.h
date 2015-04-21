#ifndef UNIQUE_H_INCLUDED
#define UNIQUE_H_INCLUDED
#include<cmath>
#include<fstream>
#include "StringBasics.h"
#include<vector>
#include "assert.h"
using namespace std;

class variant
{
public:

    string name;
    int bp;
    string chr;
    char refAllele,altAllele;
    string refAlleleString,altAlleleString;
    string MajAlleleString,MinAlleleString;

    variant(){};
    variant(string &id,string &CHR,int &BP)
    {
        name=id;
        bp=BP;
        chr=CHR;
    };
    void assignValues(string &id,string &CHR,int &BP)
    {
        name=id;
        bp=BP;
        chr=CHR;
    };
     void assignRefAlt(string &refe,string &alt)
    {
        refAlleleString=refe;
        altAlleleString=alt;
    };


};


class ReducedHaplotypeInfo
{
    public:
        int startIndex,endIndex;
        //vector<int> uniqueRep;  // in reduced state space (which is finally enumerated by this vector) ... indices of original haplotypes which are representatives
        vector<int> uniqueCardinality; // has number of representatives for each unique representative
        vector<float> InvuniqueCardinality; // has number of representatives for each unique representative
        vector<int> uniqueIndexMap; // maps to corresponding item in the uniqueRep... required to map Constants and unfold fold probs


        vector<vector<bool> > uniqueHaps;


        bool returnHapAtPosition(int i,int position)
        {
            //assert((position-startIndex)>=0);
            return uniqueHaps[i][position-startIndex];
        }



    ReducedHaplotypeInfo()
    {

        startIndex=0;
        endIndex=0;
        //uniqueRep.clear();
        uniqueCardinality.clear();
        uniqueIndexMap.clear();
        uniqueHaps.clear();
    }

};


class findUnique
{

    public:

        int N; // No. of Samples
        int max_block;  // Maximum Block Length to consider
        int M; // No. of Markers
        vector<vector<int> > uniqueCount; // matrix to store unique number of haplotypes between positions
        vector<int> minCost;
        vector<int> minAllocation;
        vector<int> optimalAllocation;
        int transFactor, cisFactor;

         findUnique()
         {


         }

         findUnique(vector<vector<char> > &hap)
        {
            N=hap.size();
            M=hap[0].size();
        }

        findUnique(vector<vector<char> > &hap,int block)
        {
            N=hap.size();
            M=hap[0].size();
            minAllocation.resize(M);
            minCost.resize(M,0);
            // Different options for choosing maximum block length
            // max_block=75;
            if(block==0)
                max_block=findMaxBlockLength(N);
            else
                max_block=block;

            uniqueCount.resize(M);
            //if fixed max block length is used, pre-allocate the memory for faster.

                for(int i=0;i<M;i++)
                uniqueCount[i].resize(max_block,-1);
        }

        int findMaxBlockLength(int total);
		void uniqueGivenPrevious(vector<vector<char> > &hap, int start, int end, vector<int> &prev_permute_index, vector<int> &prev_group_index, vector<int>& new_permute_index, vector<int>& new_group_index, vector<int>& final_permute_index);
		//MTW void uniqueGivenPrevious(vector<vector<char> > &hap, int start, int end, vector<int> &prev_permute_index, vector<int> &prev_group_index);
        ReducedHaplotypeInfo uniqueGivenPrevious(vector<vector<char> > &hap,int start,int end,vector<int> &prev_permute_index, vector<int> &prev_group_index,int posFirst);
        void createOutput(String filename);
        void unique(vector<vector<char> > &hap);
        void swap(int &a,int &b);
        void printAllocation(int position);
        void UpdateDeltaMatrix(vector<String> & haplotypes, vector<int> & index,
          vector<int> & firstDifference, int length, int blockSize,
          vector<int> & oldIndex,  vector<int> & previousPredecessor,  vector<int> & previousDifference);
         void AnalyzeBlocks(
         vector<int> & index, vector<int> & firstDifference, int length, int blockSize,
         vector<int> & cost,
         vector<int> & bestSlice, vector<int> & bestComplexity, vector<vector<int> > &bestIndex);

         double FlushBlocks(vector<int> &optEndPoints,vector<ReducedHaplotypeInfo> &HapInfo, vector<variant> &varList, int LastflushPos,vector<String> & haplotypes, vector<int> & cost,
                   vector<int> & bestComplexity, vector<int> & bestSlice, vector<vector<int> > &bestIndex);

  double FlushBlocks(vector<ReducedHaplotypeInfo> &HapInfo, vector<variant> &varList, int LastflushPos,vector<String> & haplotypes, vector<int> & cost,
                   vector<int> & bestComplexity, vector<int> & bestSlice, vector<vector<int> > &bestIndex);


        void updateCoeffs(int trans,int cis)
        {
            transFactor = trans;
            cisFactor = cis;

        }


};






#endif // UNIQUE_H_INCLUDED
