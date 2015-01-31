
#include "Unique.h"
#include "MarkovModel.h"
#include "Random.h"
#include "MarkovParameters.h"
#include<boost/random/uniform_real.hpp>
#include<boost/random/variate_generator.hpp>
#include<boost/random/mersenne_twister.hpp>


void MarkovModel::initializeMatrices(HaplotypeSet &tHap,HaplotypeSet &rHap, vector<int> &optStructure,vector<ReducedHaplotypeInfo> &StructureInfo)
{

//    Leftprob.clear();
//    leftNoRecoProb.clear();

    NoPrecisionJumps=0;
    leftProb.resize(StructureInfo.size());
    leftNoRecoProb.resize(StructureInfo.size());


    vector<double> InitialHapProb;
    for(int i=0;i<(int)leftProb.size();i++)
    {
        leftProb[i].resize(optStructure[i+1]-optStructure[i]+1);
        leftNoRecoProb[i].resize(optStructure[i+1]-optStructure[i]+1);
        InitialHapProb.resize(StructureInfo[i].uniqueCardinality.size(),1.0);
        for(int j=0;j<(int)leftProb[i].size();j++)
            {

                leftProb[i][j]=InitialHapProb;
                leftNoRecoProb[i][j]=InitialHapProb;
            }
    }


    junctionLeftProb.clear();
    junctionRightProb.clear();


    junctionLeftProb.resize(optStructure.size());
    junctionRightProb.resize(optStructure.size());


    for(int i=0;i<(int)junctionLeftProb.size();i++)
    {
        junctionLeftProb[i].resize((int)rHap.numHaplotypes,1.0);
        junctionRightProb[i].resize((int)rHap.numHaplotypes,1.0);
        //cout<<"WELL " <<junctionLeftprob[i].size()<<endl;
    }
}


double MarkovModel::CountErrors(vector<double> &probHap,
                                int position, char observed, double e,double freq,
                                ReducedHaplotypeInfo &Info)
{
    if (observed == 0)
        return e;

    double match = 0;
    double mismatch = 0;
    double background = 0;


    for (int i = 0; i < noReducedStatesCurrent; i++)
{


        if(Info.returnHapAtPosition(i,position)==observed)
            match += probHap[i];
        else
            mismatch += probHap[i];
        }


    background = (match + mismatch) * backgroundError;
    mismatch = (match + mismatch) * e *freq;
    match *= 1.0 - e;

    //cout<<" GROUP "<<mismatch<<endl;
    return mismatch / (mismatch + match + background);
}


double MarkovModel::CountRecombinants(vector<double> &from, vector<double> &to,
                                      vector<double> &probHap,double r,bool PrecisionMultiply)
{
    if (r == 0)
      return 0.0;

    double fromSum = 0.0,toSum=0.0,totalSum=0.0;

    for (int i = 0; i < noReducedStatesCurrent; i++)
    {
        fromSum += from[i];
        toSum += to[i];
        totalSum+=probHap[i];
    }

    double rsum = fromSum*toSum*r/(double)refCount;
    //cout<<" THUS = "<<rsum<<endl;

    if(PrecisionMultiply)
        return (1e15*rsum / totalSum);
    else
        return (rsum / totalSum);
}



void MarkovModel::CountExpected(HaplotypeSet &tHap,int hapID,vector<double> &rightFoldProb,
                         vector<vector<double> > &Leftprob, vector<vector<double> > &leftNoRecomProb,
                         vector<double> &rightProb, vector<double> &rightNoRecomProb,
                         vector<double> &juncLeftprob,vector<double> &juncRightProb,
                         int start,int end,ReducedHaplotypeInfo &Info,vector<vector<double> > &alleleFreq)
{

    rightProb=rightFoldProb;
    noReducedStatesCurrent=(int)rightFoldProb.size();
    vector<double> probHap(8,0.0);
    vector<double> Constants(Info.uniqueCardinality.size(),0.0);
    vector<double> tempRightProb=rightFoldProb;

    for(int i=0;i<refCount;i++)
        Constants[Info.uniqueIndexMap[i]]+=(juncLeftprob[i]*juncRightProb[i]);

    CreatePosteriorProb(Leftprob[end-start],rightProb,leftNoRecomProb[end-start],
                        rightNoRecomProb,Leftprob[0],rightFoldProb,Constants,
                        probHap,Info);





    if(!missing[end])
        empError[end]+=CountErrors(probHap,end,tHap.getScaffoldedHaplotype(hapID,end),
                               Error[end], alleleFreq[tHap.getScaffoldedHaplotype(hapID,end)][end],Info);
    else
        empError[end]+=Error[end];

    for (int markerPos=end-1; markerPos>start; markerPos--)
    {
        if (!missing[markerPos+1])
        {
            Condition(markerPos+1,rightProb,rightNoRecomProb,
                      tHap.getScaffoldedHaplotype(hapID,markerPos+1),
                      Error[markerPos+1],
                      alleleFreq[tHap.getScaffoldedHaplotype(hapID,markerPos+1)][markerPos+1],Info);
        }

        tempRightProb=rightProb;

        empRecom[markerPos]+=CountRecombinants(Leftprob[markerPos-start],rightProb,probHap,
                                               Recom[markerPos]
                                               ,PrecisionJump[markerPos+1]);

        //bool temp=;
        Transpose(tempRightProb,rightProb,rightNoRecomProb,Recom[markerPos],Info.uniqueCardinality);


        CreatePosteriorProb(Leftprob[markerPos-start],rightProb,
                            leftNoRecomProb[markerPos-start],rightNoRecomProb,Leftprob[0],rightFoldProb,Constants,probHap,Info);

        if(!missing[markerPos])
            empError[markerPos]+=CountErrors(probHap,markerPos,tHap.getScaffoldedHaplotype(hapID,markerPos),
                                         Error[markerPos], alleleFreq[tHap.getScaffoldedHaplotype(hapID,markerPos)][markerPos],Info);
        else
            empError[markerPos]+=Error[markerPos];
    }


    if (!missing[start+1])
    {

        Condition(start+1,rightProb,rightNoRecomProb,tHap.getScaffoldedHaplotype(hapID,start+1),Error[start+1],alleleFreq[tHap.getScaffoldedHaplotype(hapID,start+1)][start+1],Info);
    }
    //cout<<" THUS = "<<end-2<<empRecom[end-2]<<endl;

    tempRightProb=rightProb;
    empRecom[start]+=CountRecombinants(Leftprob[0],rightProb,probHap,Recom[start],PrecisionJump[start+1]);
    //bool temp=;
    Transpose(tempRightProb,rightProb,rightNoRecomProb,Recom[start],Info.uniqueCardinality);

    if(start==0)
        {
            CreatePosteriorProb(Leftprob[0],rightProb,leftNoRecomProb[0],rightNoRecomProb,Leftprob[0],rightFoldProb,Constants,probHap,Info);
            if(!missing[start])
                empError[start]+=CountErrors(probHap,start,tHap.getScaffoldedHaplotype(hapID,start), Error[start],
                                         alleleFreq[tHap.getScaffoldedHaplotype(hapID,start)][start],Info);
            else
                empError[start]+=Error[start];
        }
}



void MarkovModel::Impute(HaplotypeSet &tHap,vector<double> &rightFoldProb,
                         int &hapID,
                         vector<vector<double> > &Leftprob, vector<vector<double> > &leftNoRecomProb,
                         vector<double> &rightProb, vector<double> &rightNoRecomProb,
                         vector<double> &juncLeftprob,vector<double> &juncRightProb,
                         int start,int end,ReducedHaplotypeInfo &Info,int type,vector<vector<double> > &alleleFreq)
{

    rightProb=rightFoldProb;
    vector<double> tempRightProb;
    vector<double> Constants(Info.uniqueCardinality.size(),0.0);
    noReducedStatesCurrent=(int)rightFoldProb.size();

    for(int i=0;i<refCount;i++)
            Constants[Info.uniqueIndexMap[i]]+=(juncLeftprob[i]*juncRightProb[i]);

    Impute(end,tHap.getScaffoldedHaplotype(hapID,end),Leftprob[end-start],rightProb,leftNoRecomProb[end-start],rightNoRecomProb,Leftprob[0],
               rightFoldProb,Constants,Info,alleleFreq);

    for (int markerPos=end-1; markerPos>start; markerPos--)
    {
        if (!missing[markerPos+1])
        {
            Condition(markerPos+1,rightProb,
                      rightNoRecomProb,tHap.getScaffoldedHaplotype(hapID,markerPos+1),
                      Error[markerPos+1],alleleFreq[tHap.getScaffoldedHaplotype(hapID,markerPos+1)][markerPos+1],Info);
        }

        tempRightProb=rightProb;
        //bool temp=;
        Transpose(tempRightProb,rightProb,rightNoRecomProb,Recom[markerPos],
                            Info.uniqueCardinality);
        Impute(markerPos,tHap.getScaffoldedHaplotype(hapID,markerPos),Leftprob[markerPos-start],rightProb,leftNoRecomProb[markerPos-start],
                   rightNoRecomProb,Leftprob[0],rightFoldProb,Constants,Info,alleleFreq);
    }

    if (!missing[start+1])
    {
        Condition(start+1,rightProb,rightNoRecomProb,
                    tHap.getScaffoldedHaplotype(hapID,start+1),Error[start+1],
                    alleleFreq[tHap.getScaffoldedHaplotype(hapID,start+1)][start+1],Info);
    }

    tempRightProb=rightProb;
    //bool temp=;
    Transpose(tempRightProb,rightProb,rightNoRecomProb,Recom[start],Info.uniqueCardinality);


    if(start==0)
        Impute(start,tHap.getScaffoldedHaplotype(hapID,start),Leftprob[0],rightProb,leftNoRecomProb[0],
                       rightNoRecomProb,Leftprob[0],rightFoldProb,Constants,Info,alleleFreq);

}


void MarkovModel::CreatePosteriorProb(vector<double> &Leftprob,vector<double> &rightProb,
                         vector<double> &leftNoRecoProb,vector<double> &rightNoRecoProb,
                         vector<double> &leftEndProb,vector<double> &rightEndProb,
                         vector<double> &Constants,vector<double> &probHap,ReducedHaplotypeInfo &Info)
{
    probHap.clear();
    probHap.resize(Info.uniqueCardinality.size());

    double value=0.0;
    for(int i=0;i<noReducedStatesCurrent;i++)
    {
        double LeftR=Leftprob[i]-leftNoRecoProb[i];
        double RightR=rightProb[i]-rightNoRecoProb[i];
        if(Info.uniqueCardinality[i]>0)
            value    =(Constants[i]*leftNoRecoProb[i]*rightNoRecoProb[i]/(leftEndProb[i]*rightEndProb[i]))
                        +(((leftNoRecoProb[i]*RightR)+(LeftR*rightNoRecoProb[i])+(LeftR*RightR))/((double)Info.uniqueCardinality[i]));
        probHap[i]=value;
    }
}




void MarkovModel::Impute(int position, char observed,
                         vector<double> &Leftprob,vector<double> &rightProb,
                         vector<double> &leftNoRecoProb,vector<double> &rightNoRecoProb,
                         vector<double> &leftEndProb,vector<double> &rightEndProb,
                         vector<double> &Constants,ReducedHaplotypeInfo &Info,
                         vector<vector<double> > &alleleFreq)
{


    double P[8] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double value    =0.0;
    for(int i=0;i<noReducedStatesCurrent;i++)
    {
        double LeftR=Leftprob[i]-leftNoRecoProb[i];
        double RightR=rightProb[i]-rightNoRecoProb[i];
        if(Info.uniqueCardinality[i]>0)
            value    =(Constants[i]*leftNoRecoProb[i]*rightNoRecoProb[i]/(leftEndProb[i]*rightEndProb[i]))
                        +(((leftNoRecoProb[i]*RightR)+(LeftR*rightNoRecoProb[i])+(LeftR*RightR))/((double)Info.uniqueCardinality[i]));
		P[(int)Info.returnHapAtPosition(i,position)]+=value;
    }

    //cout<<"WELL = "<<Leftprob[0]<<endl;//-leftNoRecoProb[0]<<endl;

    double ptotal = P[1] + P[2] + P[3] + P[4] + P[5] + P[6] + P[7];
    double pmajor = P[(int)major[position]];
    int mle = 1;
    for (int i = 2; i <= 7; i++)
        if (P[i] >= P[mle])
            mle = i;


    char labels[]= {0, 'A', 'C', 'G', 'T', 'D', 'I', 'R'};

    imputedDose[position] += imputedHap[position]= (pmajor / ptotal);

    imputedAlleles[position] = labels[mle];
    imputedAlleleNumber[position] = mle;

    double fmatch = 1.0 / (1. - Error[position] + Error[position] * alleleFreq[major[position]][position] + backgroundError);
    double fmismatch = 1.0 / (Error[position] * alleleFreq[major[position]][position] + backgroundError);

    for (int i = 1; i <= 7; i++)
        if (observed== i)
            P[i] *= fmatch;
        else
            P[i] *= fmismatch;

   ptotal = P[1] + P[2] + P[3] + P[4] + P[5] + P[6] + P[7];


   pmajor = P[(int)major[position]];
   leaveOneOut[position] = pmajor / ptotal;
}



vector<double> MarkovModel::foldProbabilities(int bridgeIndex,ReducedHaplotypeInfo &Info,int direction,int noReference) //0 - left; 1 - right
{

    vector<double> foldProb(Info.uniqueCardinality.size(),0.0);

    if(direction==0)
    {
        for(int i=0;i<noReference;i++)
        {
            //cout<<" MAP = "<<MM.junctionLeftProb[bridgeIndex][i]<<endl;
            foldProb[Info.uniqueIndexMap[i]]+=junctionLeftProb[bridgeIndex][i];
        }


    }
    else if(direction==1)
    {
        for(int i=0;i<noReference;i++)
        {
            foldProb[Info.uniqueIndexMap[i]]+=junctionRightProb[bridgeIndex+1][i];
        }
    }

    return foldProb;
}





void MarkovModel::unfoldProbabilities(int bridgeIndex,vector<double> &recomProb, vector<double> &noRecomProb,vector<double> &foldedProb,
                                     int direction,vector<ReducedHaplotypeInfo> &StructureInfo,int noReference)
{


    for(int i=0;i<noReference;i++)
    {

        int tempIndex=StructureInfo[bridgeIndex].uniqueIndexMap[i];

        if(direction==0)
            {

                junctionLeftProb[bridgeIndex+1][i]=((recomProb[tempIndex])/(double)StructureInfo[bridgeIndex].uniqueCardinality[tempIndex])+(noRecomProb[tempIndex]*(junctionLeftProb[bridgeIndex][i]/foldedProb[tempIndex]));
            // cout<<" MAP = "<<recomProb[tempIndex]<<"\t"<<StructureInfo[bridgeIndex].uniqueCardinality[tempIndex]<<"\t"<<noRecomProb[tempIndex]<<"\t"<<MM.junctionLeftProb[bridgeIndex][i]<<"\t"<<foldedProb[tempIndex]<<endl;

            }
        if(direction==1)
            junctionRightProb[bridgeIndex][i]=((recomProb[tempIndex])/(double)StructureInfo[bridgeIndex].uniqueCardinality[tempIndex])+(noRecomProb[tempIndex]*(junctionRightProb[bridgeIndex+1][i]/foldedProb[tempIndex]));

    }


}





void MarkovModel::WalkLeft(HaplotypeSet &tHap, int &hapID,
                           vector<vector<double> > &Leftprob,vector<vector<double> > &noRecomProb,
                           vector<double> &foldedProb,int start,int end,ReducedHaplotypeInfo &Info,
                           vector<vector<double> > &alleleFreq)
{

    Leftprob[0] = foldedProb;
    noReducedStatesCurrent=(int)foldedProb.size();

    for (int markerPos=start+1; markerPos<=end; markerPos++)
    {
        noRecomProb[markerPos-start]=noRecomProb[markerPos-start-1];
        PrecisionJump[markerPos]=Transpose(Leftprob[markerPos-start-1],
                  Leftprob[markerPos-start],noRecomProb[markerPos-start],
                  Recom[markerPos-1],Info.uniqueCardinality);

        if (!missing[markerPos])
            {
               Condition(markerPos,Leftprob[markerPos-start],
                         noRecomProb[markerPos-start],
                         tHap.getScaffoldedHaplotype(hapID,markerPos),
                         Error[markerPos],
                         alleleFreq[tHap.getScaffoldedHaplotype(hapID,markerPos)][markerPos],Info);            }
//cout<<" MAP = "<< Error[markerPos]<<"\t"<<alleleFreq[tHap.getScaffoldedHaplotype(hapID][markerPos]][markerPos)<<"\n";//<<noRecomProb[tempIndex]<<"\t"<<junctionLeftprob[bridgeIndex][i]<<"\t"<<foldedProb[tempIndex]<<endl;

    }

}



void MarkovModel::Condition(int markerPos,vector<double> &Prob,
                            vector<double> &noRecomProb, char observed,double e,double freq,ReducedHaplotypeInfo &Info)
{

    if (observed ==  0)
    {
        return;
    }

    double prandom = e*freq+backgroundError;
    double pmatch = (1.0 - e)+e*freq+backgroundError;

    double P[8] = { prandom, prandom, prandom, prandom,prandom, prandom, prandom, prandom };

	P[(int)observed] = pmatch;

    for (int i = 0; i<noReducedStatesCurrent; i++)
    {

        char allele=Info.returnHapAtPosition(i,markerPos);
        Prob[i]*=P[(int)allele];
        noRecomProb[i]*=P[(int)allele];
    }
}








bool MarkovModel::Transpose(vector<double> &from,
                            vector<double> &to, vector<double> &noRecomProb,
                            double reco,vector<int> &uniqueCardinality)
{
    bool tempPrecisionJumpFlag=false;
    if (reco == 0)
    {
        to=from;
        return false;
    }

    double sum = 0.0;
    for (int i = 0; i <noReducedStatesCurrent; i++)
    {
        sum += from[i];
        noRecomProb[i]*=(1.-reco);
    }

    sum*=(reco/(double)refCount);
    double complement = 1. - reco;

    // avoid underflows
    if (sum < 1e-10)
    {
        tempPrecisionJumpFlag=true;
        sum*= 1e15;
        complement *= 1e15;
        for(int i=0;i<noReducedStatesCurrent;i++)
            noRecomProb[i]*=1e15;
        NoPrecisionJumps++;
    }

    for (int i = 0; i <noReducedStatesCurrent; i++)
    {
        to[i]=from[i]*complement+(uniqueCardinality[i]*sum);
    }
    return tempPrecisionJumpFlag;
 }





//
//
//
//void MarkovModel::ProfileModel  (HaplotypeSet &tHap, HaplotypeSet &rHap,
//                                 vector<vector<double> > &Leftprob, vector<vector<double> > &leftNoRecomProb,
//                                 vector<double> &juncLeftprobLeft,vector<double> &juncLeftprobRight,
//                                  int start,int end,ReducedHaplotypeInfo &Info,int &currentState,
//                                  boost::variate_generator<boost::mt19937&, boost::uniform_real<> > &uni)
//{
//
//    int newCurrentState=currentState;
//    int new_end=end;
//
//
//    if(end==(int)rHap.haplotypes[0].size()-1) // if last marker, needs to be separately sampled.
//    {
//        // Cumulative probability
//        double sum = 0.0;
//
//        // Sample state at the last position
//        for (int i = 0; i<(int)juncLeftprobRight.size(); i++)
//            {
//                sum += juncLeftprobRight[i];
//            }
//
//        double r = uni() * sum;
//        newCurrentState = 0;
//
//        for (sum = 0.0; newCurrentState<(int)juncLeftprobRight.size()-1 && sum<r; newCurrentState++)
//           {
//               sum = sum + juncLeftprobRight[newCurrentState];
//            }
//        if(!missing[end])
//            if(tHap.getScaffoldedHaplotype(0][end]!=rHap.haplotypes[newCurrentState][end))
//                countError[end]++;
//
//
//
////        cout<<"\n END = "<<newCurrentState<<" ";
//        }
//
//
//        new_end--;
//    //if(!missing[end])
//
//
//    for (int currentMarker = new_end; currentMarker>=start; currentMarker--)
//    {
//
//        vector<double> recomProb(leftNoRecomProb.size());
//        for(int i=0;i<(int)recomProb.size();i++)
//        {
//            recomProb[i]=Leftprob[new_end-currentMarker][i]-leftNoRecomProb[new_end-currentMarker][i];
//        }
//
//        int tempIndex=Info.uniqueIndexMap[newCurrentState];
//
//        double currentUnfoldProb=((recomProb[tempIndex])/(double)Info.uniqueCardinality[tempIndex])+(leftNoRecomProb[new_end-currentMarker][tempIndex]*(juncLeftprobLeft[newCurrentState]/Leftprob[new_end-currentMarker][tempIndex]));
//
//
//
//        double sum = 0.0;
//
//        for (int i = 0; i < Leftprob[0].size(); i++)
//            sum += Leftprob[currentMarker-start][i];
//
//        double rec = sum*Recom[currentMarker]/(double)rHap.haplotypes.size();
//        double norec = currentUnfoldProb * (1.0 - Recom[currentMarker]);
//
//
//        sum=rec+norec;
//
//        double r = uni() * sum;
//
//
//        if(r>norec)
//        {
//            countRecom[currentMarker]++;
//
//
//            //cout<<" BREAK ";
//
//            r=uni()*(int)rHap.haplotypes.size();
//            newCurrentState = 0;
//            for ( sum = 0.0 ; newCurrentState < (int)rHap.haplotypes.size() - 1; newCurrentState++)
//            if ( (sum += 1) > r)
//            break;
//        }
//
//         if(!missing[currentMarker])
//            if(tHap.getScaffoldedHaplotype(0][currentMarker]!=rHap.haplotypes[newCurrentState][currentMarker))
//                countError[currentMarker]++;
//
//        //cout<<"  "<<newCurrentState<<" ";
//
//
//    }
//
//
//
//    currentState=newCurrentState;
//
//}
//






