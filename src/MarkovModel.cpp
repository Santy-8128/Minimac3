
#include "MarkovModel.h"

void MarkovModel::initializeMatrices(HaplotypeSet &tHap,HaplotypeSet &rHap, vector<int> &optStructure,
                                     vector<ReducedHaplotypeInfo> &StructureInfo)
{
    NoPrecisionJumps=0;
    leftProb.resize(StructureInfo.size());
    leftNoRecoProb.resize(StructureInfo.size());


    vector<float> InitialHapProb;
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
    }
}


double MarkovModel::CountErrors(vector<float> &probHap,
                                int position, bool observed, double e,double freq,
                                ReducedHaplotypeInfo &Info)
{

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

    return mismatch / (mismatch + match + background);
}


double MarkovModel::CountRecombinants(vector<float> &from, vector<float> &to,
                                      vector<float> &probHap,double r,bool PrecisionMultiply)
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

    double rsum = fromSum*r*toSum/(double)refCount;

    if(PrecisionMultiply)
        return (1e15*rsum / totalSum);
    else
        return (rsum / totalSum);
}



void MarkovModel::CountExpected(HaplotypeSet &tHap,int hapID,vector<float> &rightFoldProb,
                         vector<vector<float> > &Leftprob, vector<vector<float> > &leftNoRecomProb,
                         vector<float> &rightProb, vector<float> &rightNoRecomProb,
                         vector<float> &juncLeftprob,vector<float> &juncRightProb,
                         int start,int end,ReducedHaplotypeInfo &Info,vector<double> &alleleFreq)
{

    rightProb=rightFoldProb;
    noReducedStatesCurrent=(int)rightFoldProb.size();
    vector<float> probHap(8,0.0);
    vector<float> Constants(Info.uniqueCardinality.size(),0.0);
    vector<float> tempRightProb=rightFoldProb;

    for(int i=0;i<refCount;i++)
        Constants[Info.uniqueIndexMap[i]]+=(juncLeftprob[i]*juncRightProb[i]);

    CreatePosteriorProb(Leftprob[end-start],rightProb,leftNoRecomProb[end-start],
                        rightNoRecomProb,Leftprob[0],rightFoldProb,Constants,
                        probHap,Info);


    if(!missing[end] && !tHap.getMissingScaffoldedHaplotype(hapID,end))
         empError[end]+=CountErrors(probHap,end,tHap.getScaffoldedHaplotype(hapID,end),
                               Error[end],
                               tHap.getScaffoldedHaplotype(hapID,end)? alleleFreq[end] : 1-alleleFreq[end],Info);
    else
        empError[end]+=Error[end];


    for (int markerPos=end-1; markerPos>start; markerPos--)
    {
        if (!missing[markerPos+1] && !tHap.getMissingScaffoldedHaplotype(hapID,markerPos+1))
        {

              Condition(markerPos+1,rightProb,rightNoRecomProb,
                      tHap.getScaffoldedHaplotype(hapID,markerPos+1),
                      Error[markerPos+1],
                       tHap.getScaffoldedHaplotype(hapID,markerPos+1)?
                          alleleFreq[markerPos+1] : 1-alleleFreq[markerPos+1],Info);
        }


        tempRightProb=rightProb;

        empRecom[markerPos]+=CountRecombinants(Leftprob[markerPos-start],rightProb,probHap,
                                               Recom[markerPos]
                                               ,PrecisionJump[markerPos+1]);

        Transpose(tempRightProb,rightProb,rightNoRecomProb,Recom[markerPos],Info.uniqueCardinality);


        CreatePosteriorProb(Leftprob[markerPos-start],rightProb,
                            leftNoRecomProb[markerPos-start],rightNoRecomProb,Leftprob[0],rightFoldProb,Constants,probHap,Info);


        if(!missing[markerPos] && !tHap.getMissingScaffoldedHaplotype(hapID,markerPos))
            empError[markerPos]+=CountErrors(probHap,markerPos,tHap.getScaffoldedHaplotype(hapID,markerPos),
                               Error[markerPos],
                               tHap.getScaffoldedHaplotype(hapID,markerPos)? alleleFreq[markerPos] : 1-alleleFreq[markerPos]
                              ,Info);
        else
            empError[markerPos]+=Error[markerPos];

    }

    if (!missing[start+1] && !tHap.getMissingScaffoldedHaplotype(hapID,start+1))
        {

              Condition(start+1,rightProb,rightNoRecomProb,
                      tHap.getScaffoldedHaplotype(hapID,start+1),
                      Error[start+1],
                       tHap.getScaffoldedHaplotype(hapID,start+1)?
                          alleleFreq[start+1] : 1-alleleFreq[start+1],Info);
        }

    tempRightProb=rightProb;
    empRecom[start]+=CountRecombinants(Leftprob[0],rightProb,probHap,Recom[start],PrecisionJump[start+1]);
    Transpose(tempRightProb,rightProb,rightNoRecomProb,Recom[start],Info.uniqueCardinality);

    if(start==0)
        {
            CreatePosteriorProb(Leftprob[0],rightProb,leftNoRecomProb[0],rightNoRecomProb,Leftprob[0],rightFoldProb,Constants,probHap,Info);

            if(!missing[start] && !tHap.getMissingScaffoldedHaplotype(hapID,start))
                 empError[start]+=CountErrors(probHap,start,tHap.getScaffoldedHaplotype(hapID,start),
                                       Error[start],
                                       tHap.getScaffoldedHaplotype(hapID,start)? alleleFreq[start] : 1-alleleFreq[start]
                                       ,Info);
            else
                empError[start]+=Error[start];

        }
}



void MarkovModel::Impute(HaplotypeSet &tHap,vector<float> &rightFoldProb,
                         int &hapID,
                         vector<vector<float> > &Leftprob, vector<vector<float> > &leftNoRecomProb,
                         vector<float> &rightProb, vector<float> &rightNoRecomProb,
                         vector<float> &juncLeftprob,vector<float> &juncRightProb,
                         int start,int end,ReducedHaplotypeInfo &Info,int type,vector<double>  &alleleFreq)
{

    rightProb=rightFoldProb;
    vector<float> tempRightProb;
    vector<float> Constants(Info.uniqueCardinality.size(),0.0);
    noReducedStatesCurrent=(int)rightFoldProb.size();

    for(int i=0;i<refCount;i++)
            Constants[Info.uniqueIndexMap[i]]+=(juncLeftprob[i]*juncRightProb[i]);

    Impute(end,tHap.getScaffoldedHaplotype(hapID,end),tHap.getMissingScaffoldedHaplotype(hapID,end),Leftprob[end-start],rightProb,leftNoRecomProb[end-start],rightNoRecomProb,Leftprob[0],
               rightFoldProb,Constants,Info,alleleFreq);

    for (int markerPos=end-1; markerPos>start; markerPos--)
    {

         if (!missing[markerPos+1] && !tHap.getMissingScaffoldedHaplotype(hapID,markerPos+1))
        {

              Condition(markerPos+1,rightProb,rightNoRecomProb,
                      tHap.getScaffoldedHaplotype(hapID,markerPos+1),
                      Error[markerPos+1],
                       tHap.getScaffoldedHaplotype(hapID,markerPos+1)?
                          alleleFreq[markerPos+1] : 1-alleleFreq[markerPos+1],Info);
        }


        tempRightProb=rightProb;
        Transpose(tempRightProb,rightProb,rightNoRecomProb,Recom[markerPos],
                            Info.uniqueCardinality);
        Impute(markerPos,tHap.getScaffoldedHaplotype(hapID,markerPos),tHap.getMissingScaffoldedHaplotype(hapID,markerPos), Leftprob[markerPos-start],rightProb,leftNoRecomProb[markerPos-start],
                   rightNoRecomProb,Leftprob[0],rightFoldProb,Constants,Info,alleleFreq);
    }

    if (!missing[start+1] && !tHap.getMissingScaffoldedHaplotype(hapID,start+1))
        {

              Condition(start+1,rightProb,rightNoRecomProb,
                      tHap.getScaffoldedHaplotype(hapID,start+1),
                      Error[start+1],
                       tHap.getScaffoldedHaplotype(hapID,start+1)?
                          alleleFreq[start+1] : 1-alleleFreq[start+1],Info);
        }

    tempRightProb=rightProb;
    Transpose(tempRightProb,rightProb,rightNoRecomProb,Recom[start],Info.uniqueCardinality);


    if(start==0)
        Impute(start,tHap.getScaffoldedHaplotype(hapID,start),tHap.getMissingScaffoldedHaplotype(hapID,start),Leftprob[0],rightProb,leftNoRecomProb[0],
                       rightNoRecomProb,Leftprob[0],rightFoldProb,Constants,Info,alleleFreq);

}


void MarkovModel::CreatePosteriorProb(vector<float> &Leftprob,vector<float> &rightProb,
                         vector<float> &leftNoRecoProb,vector<float> &rightNoRecoProb,
                         vector<float> &leftEndProb,vector<float> &rightEndProb,
                         vector<float> &Constants,vector<float> &probHap,ReducedHaplotypeInfo &Info)
{
    probHap.clear();
    probHap.resize(Info.uniqueCardinality.size());

    double value=0.0;
    for(int i=0;i<noReducedStatesCurrent;i++)
    {

        value = Constants[i]*(leftNoRecoProb[i]*rightNoRecoProb[i]/(leftEndProb[i]*rightEndProb[i]))
            +(Leftprob[i]*rightProb[i]-leftNoRecoProb[i]*rightNoRecoProb[i])*(Info.InvuniqueCardinality[i]);

        probHap[i]=value;
    }
}


void MarkovModel::ImputeOld(int position, bool observed, bool observedMiss,
                         vector<float> &Leftprob,vector<float> &rightProb,
                         vector<float> &leftNoRecoProb,vector<float> &rightNoRecoProb,
                         vector<float> &leftEndProb,vector<float> &rightEndProb,
                         vector<float> &Constants,ReducedHaplotypeInfo &Info,
                         vector<double> &alleleFreq)
{

    float Pref=0.0,Palt=0.0;

    double value=0.0;

    for(int i=0;i<noReducedStatesCurrent;)
    {
        bool hp = Info.returnHapAtPosition(i,position);

        float LeftR=Leftprob[i]-leftNoRecoProb[i];
        float RightR=rightProb[i]-rightNoRecoProb[i];
        if(Info.uniqueCardinality[i]>0)
            value    =((Constants[i]/leftEndProb[i])
            *(leftNoRecoProb[i]/rightEndProb[i])*rightNoRecoProb[i])
                        +(((leftNoRecoProb[i]*RightR)
                           +(LeftR*rightNoRecoProb[i])
                           +(LeftR*RightR))/((double)Info.uniqueCardinality[i]));


        float pp=0.0;

        if(hp)
            pp = Palt + value;
        else
            pp = Pref + value;


		i++;
		while ((i < noReducedStatesCurrent) && (hp == Info.returnHapAtPosition(i,position)))
        {
            LeftR=Leftprob[i]-leftNoRecoProb[i];
            RightR=rightProb[i]-rightNoRecoProb[i];

            if(Info.uniqueCardinality[i]>0)
                value    =((Constants[i]/leftEndProb[i])
                        *(leftNoRecoProb[i]/rightEndProb[i])*rightNoRecoProb[i])
                        +(((leftNoRecoProb[i]*RightR)
                           +(LeftR*rightNoRecoProb[i])
                           +(LeftR*RightR))/((double)Info.uniqueCardinality[i]));
            pp += value;
            i++;
        }

        if(hp)
            Palt = pp ;
        else
            Pref = pp;
    }



    float ptotal =Pref+Palt;

    float pmajor = major[position]? Palt : Pref;


    bool mle = false;

    if(Pref<Palt)
    {
        mle=true;
    }



    imputedDose[position] += imputedHap[position]= (pmajor / ptotal);


    imputedAlleles[position] = mle;
    imputedAlleleNumber[position] = mle;

    double fmatch = 1.0 / (1. - Error[position] + Error[position] * ( major[position]? alleleFreq[position] : 1-alleleFreq[position]  ) + backgroundError );
    double fmismatch = 1.0 / (Error[position] * ( major[position]? alleleFreq[position] : 1-alleleFreq[position]  ) + backgroundError);




    if(!observedMiss)
    {
        if(observed)
        {
            Palt *= fmatch;
            Pref *= fmismatch;
        }
        else
        {
            Pref *= fmatch;
            Palt *= fmismatch;
        }
    }

    ptotal =Pref+Palt;
    pmajor = major[position]? Palt : Pref;
    leaveOneOut[position] = pmajor / ptotal;


}




void MarkovModel::Impute(int position, bool observed, bool observedMiss,
                         vector<float> &Leftprob,vector<float> &rightProb,
                         vector<float> &leftNoRecoProb,vector<float> &rightNoRecoProb,
                         vector<float> &leftEndProb,vector<float> &rightEndProb,
                         vector<float> &Constants,ReducedHaplotypeInfo &Info,
                         vector<double> &alleleFreq)
{



    float Pref=0.0,Palt=0.0;

    float *value = (float *)alloca(noReducedStatesCurrent*sizeof(float));
    for(int i=0; i<noReducedStatesCurrent; i++)
    {
        // careful: order of operations is important to avoid overflows
        value[i] = Constants[i]*(leftNoRecoProb[i]*rightNoRecoProb[i]/(leftEndProb[i]*rightEndProb[i]))
            +(Leftprob[i]*rightProb[i]-leftNoRecoProb[i]*rightNoRecoProb[i])*(Info.InvuniqueCardinality[i]);
    }

    for (int i=0; i<noReducedStatesCurrent;)
    {
        bool hp = Info.returnHapAtPosition(i,position);

        float pp=0.0;

        pp= value[i] + (hp? Palt:Pref);

        i++;
        while ((i < noReducedStatesCurrent) && (hp == Info.returnHapAtPosition(i,position)))
        {
            pp += value[i];
            i++;
        }

        if(hp)
            Palt = pp ;
        else
            Pref = pp;
    }


    float ptotal =Pref+Palt;

    float pmajor = major[position]? Palt : Pref;


    bool mle = false;

    if(Pref<Palt)
    {
        mle=true;
    }



    imputedDose[position] += imputedHap[position]= (pmajor / ptotal);
    imputedAlleles[position] = mle;
    imputedAlleleNumber[position] = mle;

    double fmatch = 1.0 / (1. - Error[position] + Error[position] * ( major[position]? alleleFreq[position] : 1-alleleFreq[position]  ) + backgroundError );
    double fmismatch = 1.0 / (Error[position] * ( major[position]? alleleFreq[position] : 1-alleleFreq[position]  ) + backgroundError);




    if(!observedMiss)
    {
        if(observed)
        {
            Palt *= fmatch;
            Pref *= fmismatch;
        }
        else
        {
            Pref *= fmatch;
            Palt *= fmismatch;
        }
    }

    ptotal =Pref+Palt;
    pmajor = major[position]? Palt : Pref;
    leaveOneOut[position] = pmajor / ptotal;


}



void MarkovModel::foldProbabilities(vector<float> &foldProb,int bridgeIndex,ReducedHaplotypeInfo &Info,int direction,int noReference) //0 - left; 1 - right
{

    foldProb.clear();
    foldProb.resize(Info.uniqueCardinality.size(),0.0);
    vector<int> *TempuniqueIndexMap=&Info.uniqueIndexMap;

    if(direction==0)
    {
        vector<float> *PrevjunctionLeftProb=&junctionLeftProb[bridgeIndex];
        for(int i=0;i<noReference;i++)
        {
            foldProb[(*TempuniqueIndexMap)[i]]+=(*PrevjunctionLeftProb)[i];
        }
    }
    else if(direction==1)
    {
        vector<float> *PrevjunctionRightProb=&junctionRightProb[bridgeIndex+1];
        for(int i=0;i<noReference;i++)
        {
            foldProb[(*TempuniqueIndexMap)[i]]+=(*PrevjunctionRightProb)[i];
        }
    }
}





void MarkovModel::unfoldProbabilities(int bridgeIndex,vector<float> &recomProb,
                                       vector<float> &noRecomProb,vector<float> &foldedProb,
                                     int direction,vector<ReducedHaplotypeInfo> &StructureInfo,
                                     int noReference)
{
    int N = StructureInfo[bridgeIndex].uniqueCardinality.size();

    float *adj_rec = (float *)alloca(N*sizeof(float));
    float *adj_norec = (float *)alloca(N*sizeof(float));
    ReducedHaplotypeInfo &thisInfo = StructureInfo[bridgeIndex];
    for (int i=0; i<N; i++)
    {
        adj_rec[i] = recomProb[i] * thisInfo.InvuniqueCardinality[i];
    }
    for (int i=0; i<N; i++)
    {
        adj_norec[i] = noRecomProb[i] / foldedProb[i];
    }
    vector<float> &prev = direction ? junctionRightProb[bridgeIndex+1] : junctionLeftProb[bridgeIndex];
    vector<float> &next = direction ? junctionRightProb[bridgeIndex] : junctionLeftProb[bridgeIndex+1];
    for (int i=0; i<noReference; i++)
    {
        int m = thisInfo.uniqueIndexMap[i];
        next[i] = adj_rec[m] + adj_norec[m]*prev[i];
    }

}



void MarkovModel::unfoldProbabilitiesOld(int bridgeIndex,vector<float> &recomProb,
                                       vector<float> &noRecomProb,vector<float> &foldedProb,
                                     int direction,vector<ReducedHaplotypeInfo> &StructureInfo,
                                     int noReference)
{

    ReducedHaplotypeInfo *TempStructureInfo=&StructureInfo[bridgeIndex];
    vector<int> *TempuniqueIndexMap=&TempStructureInfo->uniqueIndexMap;
    vector<int> *TempuniqueCardinality=&TempStructureInfo->uniqueCardinality;

    if(direction==0)
    {
        vector<float> *PrevjunctionLeftProb=&junctionLeftProb[bridgeIndex];
        vector<float> *NewjunctionLeftProb=&junctionLeftProb[bridgeIndex+1];
        for(int i=0;i<noReference;i++)
        {
            int tempIndex=(*TempuniqueIndexMap)[i];
            (*NewjunctionLeftProb)[i]=((recomProb[tempIndex])/(double)(*TempuniqueCardinality)[tempIndex])
                +(noRecomProb[tempIndex]*((*PrevjunctionLeftProb)[i]/foldedProb[tempIndex]));
        }
    }

    if(direction==1)
    {
        vector<float> *PrevjunctionRightProb=&junctionRightProb[bridgeIndex+1];
        vector<float> *NewjunctionRightProb=&junctionRightProb[bridgeIndex];
        for(int i=0;i<noReference;i++)
        {
            int tempIndex=(*TempuniqueIndexMap)[i];
            (*NewjunctionRightProb)[i]=((recomProb[tempIndex])/(double)(*TempuniqueCardinality)[tempIndex])+
            (noRecomProb[tempIndex]*((*PrevjunctionRightProb)[i]/foldedProb[tempIndex]));
        }
    }



}





void MarkovModel::WalkLeft(HaplotypeSet &tHap, int &hapID,
                           vector<vector<float> > &Leftprob,vector<vector<float> > &noRecomProb,
                           vector<float> &foldedProb,int start,int end,ReducedHaplotypeInfo &Info,
                           vector<double> &alleleFreq)
{

    Leftprob[0] = foldedProb;
    noReducedStatesCurrent=(int)foldedProb.size();

    for (int markerPos=start+1; markerPos<=end; markerPos++)
    {
        noRecomProb[markerPos-start]=noRecomProb[markerPos-start-1];
        PrecisionJump[markerPos]=Transpose(Leftprob[markerPos-start-1],
                  Leftprob[markerPos-start],noRecomProb[markerPos-start],
                  Recom[markerPos-1],Info.uniqueCardinality);


        if (!missing[markerPos] && !tHap.getMissingScaffoldedHaplotype(hapID,markerPos))
        {




           if(!tHap.getMissingScaffoldedHaplotype(hapID,markerPos))
                Condition(markerPos,Leftprob[markerPos-start],
                     noRecomProb[markerPos-start],
                     tHap.getScaffoldedHaplotype(hapID,markerPos),
                     Error[markerPos],
                     tHap.getScaffoldedHaplotype(hapID,markerPos)?
                          alleleFreq[markerPos] : 1-alleleFreq[markerPos],Info);

        }
    }

}



void MarkovModel::Condition(int markerPos,vector<float> &Prob,
                            vector<float> &noRecomProb, bool observed,double e,double freq,ReducedHaplotypeInfo &Info)
{

    double prandom = e*freq+backgroundError;
    double pmatch = (1.0 - e)+e*freq+backgroundError;

    for (int i = 0; i<noReducedStatesCurrent; i++)
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




bool MarkovModel::Transpose(vector<float> &from,
                            vector<float> &to, vector<float> &noRecomProb,
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

