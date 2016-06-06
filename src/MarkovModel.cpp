
#include "MarkovModel.h"




void MarkovModel::initializeMatrices(HaplotypeSet & rHap,HaplotypeSet &tHap)
{
    if(LowMemory)
    {

        junctionLeftFoldProb.clear();
        junctionLeftFoldProb.resize(rHap.NoBlocks);
        for(int i=0;i<rHap.NoBlocks;i++)
        {
            junctionLeftFoldProb[i].resize(rHap.maxRepSize);
        }

        ThisBlockLeftProb.clear();
        ThisBlockLeftProb.resize(rHap.maxBlockSize);
        for(int i=0;i<rHap.maxBlockSize;i++)
        {
            ThisBlockLeftProb[i].resize(rHap.maxRepSize);
        }

    }
    else
    {
        leftProb.resize(rHap.NoBlocks);
        for(int i=0;i<rHap.NoBlocks;i++)
        {
            ReducedHaplotypeInfo &TempBlock=rHap.ReducedStructureInfo[i];
            vector<vector<float> > &TempLeft=leftProb[i];
            TempLeft.resize(TempBlock.BlockSize);
            for(int j=0;j<TempBlock.BlockSize;j++)
            {
                TempLeft[j].resize(TempBlock.RepSize);
            }
        }
    }


    CurrentLeftNoRecoProb.clear();
    CurrentLeftNoRecoProb.resize(rHap.maxRepSize);

    ThisBlockLeftNoRecoProb.clear();
    ThisBlockLeftNoRecoProb.resize(rHap.maxBlockSize);
    for(int i=0;i<rHap.maxBlockSize;i++)
    {
        ThisBlockLeftNoRecoProb[i].resize(rHap.maxRepSize);
    }

    junctionLeftProb.clear();
    junctionLeftProb.resize(rHap.NoBlocks+1);
    PrevjunctionRightProb.clear();
    PrevjunctionRightProb.resize(rHap.numHaplotypes);

    for(int i=0;i<=rHap.NoBlocks;i++)
    {
        junctionLeftProb[i].resize(rHap.numHaplotypes);
    }

    probHap.resize(rHap.maxRepSize);
    Constants.resize(rHap.maxRepSize);
    tempRightProb.reserve(rHap.maxRepSize);



}

void MarkovModel::ReinitializeMatrices()
{

    NoPrecisionJumps=0;
    fill(junctionLeftProb[0].begin(), junctionLeftProb[0].end(), 1.0);
    fill(PrevjunctionRightProb.begin(), PrevjunctionRightProb.end(), 1.0);

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



void MarkovModel::ReCreateLeftNoRecoProb(HaplotypeSet &tHap, int &hapID,
                           int group, ReducedHaplotypeInfo &Info,
                           vector<double> &alleleFreq)
{

    int &Start=Info.startIndex;
    int &End=Info.endIndex;
    noReducedStatesCurrent=Info.RepSize;
    ThisBlockLeftNoRecoProb[0]=leftProb[group][0];

    for (int markerPos=Start+1; markerPos<=End; markerPos++)
    {
        vector<float> &NextnoRecomProb = ThisBlockLeftNoRecoProb[markerPos-Start];
        double complement = 1. - Recom[markerPos-1];
        NextnoRecomProb=ThisBlockLeftNoRecoProb[markerPos-Start-1];
        double freq=tHap.getScaffoldedHaplotype(hapID,markerPos)? alleleFreq[markerPos] : 1-alleleFreq[markerPos];
        double e=Error[markerPos];
        bool observed=tHap.getScaffoldedHaplotype(hapID,markerPos);
        double prandom = e*freq+backgroundError;
        double pmatch = (1.0 - e)+e*freq+backgroundError;

        for (int i = 0; i <noReducedStatesCurrent; i++)
        {
            NextnoRecomProb[i]*=(complement);
        }

        if (PrecisionJump[markerPos])
        {
            for (int i = 0; i <noReducedStatesCurrent; i++)
            {
                NextnoRecomProb[i]*=(1e15);
            }
        }

        if (!missing[markerPos] && !tHap.getMissingScaffoldedHaplotype(hapID,markerPos))
        {
            for (int i = 0; i<noReducedStatesCurrent; i++)
                {
                    NextnoRecomProb[i]*=(Info.returnHapAtPosition(i,markerPos)==observed)?pmatch:prandom;
                }
        }

    }
}





void MarkovModel::ReCreateBothLeftProb(HaplotypeSet &tHap, int &hapID,
                           int group, ReducedHaplotypeInfo &Info,
                           vector<double> &alleleFreq)
{

    vector<vector<float> > &Leftprob = ThisBlockLeftProb;
    Leftprob[0]=junctionLeftFoldProb[group];
    ThisBlockLeftNoRecoProb[0]=Leftprob[0];

    int &Start=Info.startIndex;
    int &End=Info.endIndex;

    noReducedStatesCurrent=Info.RepSize;

    for (int markerPos=Start+1; markerPos<=End; markerPos++)
    {
        ThisBlockLeftNoRecoProb[markerPos-Start]=ThisBlockLeftNoRecoProb[markerPos-Start-1];
        Transpose(Leftprob[markerPos-Start-1],
                  Leftprob[markerPos-Start],ThisBlockLeftNoRecoProb[markerPos-Start],
                  Recom[markerPos-1],Info.uniqueCardinality);

        if (!missing[markerPos] && !tHap.getMissingScaffoldedHaplotype(hapID,markerPos))
        {

           if(!tHap.getMissingScaffoldedHaplotype(hapID,markerPos))
                Condition(markerPos,Leftprob[markerPos-Start],
                     ThisBlockLeftNoRecoProb[markerPos-Start],
                     tHap.getScaffoldedHaplotype(hapID,markerPos),
                     Error[markerPos],
                     tHap.getScaffoldedHaplotype(hapID,markerPos)?
                          alleleFreq[markerPos] : 1-alleleFreq[markerPos],Info);

        }

    }

}




void MarkovModel::CountExpected(HaplotypeSet &tHap,int hapID,int group,
                                   vector<float> &PrevRightFoldedProb,
                                    vector<float> &CurrentRightProb, vector<float> &CurrentNoRecoRightProb,
                                    ReducedHaplotypeInfo &Info,vector<double> &alleleFreq)
{


    vector<float> &juncLeftprob = junctionLeftProb[group];
    vector<float> &juncRightProb = PrevjunctionRightProb;

    if(LowMemory)
        ReCreateBothLeftProb(tHap,hapID,group,Info,alleleFreq);
    else
        ReCreateLeftNoRecoProb(tHap,hapID,group,Info,alleleFreq);

    vector<vector<float> > &Leftprob = LowMemory? ThisBlockLeftProb: leftProb[group];
    vector<vector<float> > &leftNoRecomProb= ThisBlockLeftNoRecoProb;


    int &Start=Info.startIndex;
    int &End=Info.endIndex;

    CurrentRightProb=PrevRightFoldedProb;
    CurrentNoRecoRightProb=PrevRightFoldedProb;
    noReducedStatesCurrent=Info.RepSize;

    fill(Constants.begin(), Constants.end(), 0.0);
    for(int i=0;i<refCount;i++)
        Constants[Info.uniqueIndexMap[i]]+=(juncLeftprob[i]*juncRightProb[i]);


    CreatePosteriorProb(Leftprob[End-Start],CurrentRightProb,leftNoRecomProb[End-Start],
                        CurrentNoRecoRightProb,Leftprob[0],PrevRightFoldedProb,Constants,
                        probHap,Info);

    if(!missing[End] && !tHap.getMissingScaffoldedHaplotype(hapID,End))
         empError[End]+=CountErrors(probHap,End,tHap.getScaffoldedHaplotype(hapID,End),
                               Error[End],
                               tHap.getScaffoldedHaplotype(hapID,End)? alleleFreq[End] : 1-alleleFreq[End],Info);
    else
        empError[End]+=Error[End];


    for (int markerPos=End-1; markerPos>Start; markerPos--)
    {
        if (!missing[markerPos+1] && !tHap.getMissingScaffoldedHaplotype(hapID,markerPos+1))
        {

              Condition(markerPos+1,CurrentRightProb,CurrentNoRecoRightProb,
                      tHap.getScaffoldedHaplotype(hapID,markerPos+1),
                      Error[markerPos+1],
                       tHap.getScaffoldedHaplotype(hapID,markerPos+1)?
                          alleleFreq[markerPos+1] : 1-alleleFreq[markerPos+1],Info);
        }
        tempRightProb=CurrentRightProb;

        empRecom[markerPos]+=CountRecombinants(Leftprob[markerPos-Start],CurrentRightProb,probHap,
                                               Recom[markerPos]
                                               ,PrecisionJump[markerPos+1]);

        Transpose(tempRightProb,CurrentRightProb,CurrentNoRecoRightProb,Recom[markerPos],Info.uniqueCardinality);

        CreatePosteriorProb(Leftprob[markerPos-Start],CurrentRightProb,
                            leftNoRecomProb[markerPos-Start],CurrentNoRecoRightProb,
                            Leftprob[0],PrevRightFoldedProb,Constants,probHap,Info);


        if(!missing[markerPos] && !tHap.getMissingScaffoldedHaplotype(hapID,markerPos))
            empError[markerPos]+=CountErrors(probHap,markerPos,
                                             tHap.getScaffoldedHaplotype(hapID,markerPos),
                               Error[markerPos],
                               tHap.getScaffoldedHaplotype(hapID,markerPos)?
                                              alleleFreq[markerPos] : 1-alleleFreq[markerPos],Info);
        else
            empError[markerPos]+=Error[markerPos];

    }

    if (!missing[Start+1] && !tHap.getMissingScaffoldedHaplotype(hapID,Start+1))
        {

              Condition(Start+1,CurrentRightProb,CurrentNoRecoRightProb,
                      tHap.getScaffoldedHaplotype(hapID,Start+1),
                      Error[Start+1],
                       tHap.getScaffoldedHaplotype(hapID,Start+1)?
                          alleleFreq[Start+1] : 1-alleleFreq[Start+1],Info);
        }

    tempRightProb=CurrentRightProb;
    empRecom[Start]+=CountRecombinants(Leftprob[0],CurrentRightProb,probHap,Recom[Start],PrecisionJump[Start+1]);
    Transpose(tempRightProb,CurrentRightProb,CurrentNoRecoRightProb,Recom[Start],Info.uniqueCardinality);

    if(Start==0)
        {

            CreatePosteriorProb(Leftprob[0],CurrentRightProb,leftNoRecomProb[0],
                                CurrentNoRecoRightProb,Leftprob[0],PrevRightFoldedProb,
                                Constants,probHap,Info);

            if(!missing[Start] && !tHap.getMissingScaffoldedHaplotype(hapID,Start))
                 empError[Start]+=CountErrors(probHap,Start,tHap.getScaffoldedHaplotype(hapID,Start),
                                       Error[Start],
                                       tHap.getScaffoldedHaplotype(hapID,Start)? alleleFreq[Start] : 1-alleleFreq[Start]
                                       ,Info);
            else
                empError[Start]+=Error[Start];

        }
}



void MarkovModel::Impute(HaplotypeSet &tHap,int hapID,int group,
                                   vector<float> &PrevRightFoldedProb,
                                    vector<float> &CurrentRightProb, vector<float> &CurrentNoRecoRightProb,
                                    ReducedHaplotypeInfo &Info,vector<double> &alleleFreq)
{

    vector<float> &juncLeftprob = junctionLeftProb[group];
    vector<float> &juncRightProb = PrevjunctionRightProb;


    if(LowMemory)
        ReCreateBothLeftProb(tHap,hapID,group,Info,alleleFreq);
    else
        ReCreateLeftNoRecoProb(tHap,hapID,group,Info,alleleFreq);

    vector<vector<float> > &Leftprob = LowMemory? ThisBlockLeftProb: leftProb[group];
    vector<vector<float> > &leftNoRecomProb= ThisBlockLeftNoRecoProb;

    int &start=Info.startIndex;
    int &end=Info.endIndex;


    CurrentRightProb=PrevRightFoldedProb;
    CurrentNoRecoRightProb=PrevRightFoldedProb;


    fill(Constants.begin(), Constants.end(), 0.0);
    noReducedStatesCurrent=Info.RepSize;
    for(int i=0;i<refCount;i++)
            Constants[Info.uniqueIndexMap[i]]+=(juncLeftprob[i]*juncRightProb[i]);

    Impute(end,tHap.getScaffoldedHaplotype(hapID,end),tHap.getMissingScaffoldedHaplotype(hapID,end),Leftprob[end-start],
           CurrentRightProb,leftNoRecomProb[end-start],CurrentNoRecoRightProb,Leftprob[0],
               PrevRightFoldedProb,Constants,Info,alleleFreq);

    for (int markerPos=end-1; markerPos>start; markerPos--)
    {

         if (!missing[markerPos+1] && !tHap.getMissingScaffoldedHaplotype(hapID,markerPos+1))
        {

              Condition(markerPos+1,CurrentRightProb,CurrentNoRecoRightProb,
                      tHap.getScaffoldedHaplotype(hapID,markerPos+1),
                      Error[markerPos+1],
                       tHap.getScaffoldedHaplotype(hapID,markerPos+1)?
                          alleleFreq[markerPos+1] : 1-alleleFreq[markerPos+1],Info);
        }


        tempRightProb=CurrentRightProb;
        Transpose(tempRightProb,CurrentRightProb,CurrentNoRecoRightProb,Recom[markerPos],
                            Info.uniqueCardinality);
        Impute(markerPos,tHap.getScaffoldedHaplotype(hapID,markerPos),tHap.getMissingScaffoldedHaplotype(hapID,markerPos),
               Leftprob[markerPos-start],CurrentRightProb,leftNoRecomProb[markerPos-start],
                   CurrentNoRecoRightProb,Leftprob[0],PrevRightFoldedProb,Constants,Info,alleleFreq);
    }

    if (!missing[start+1] && !tHap.getMissingScaffoldedHaplotype(hapID,start+1))
        {

              Condition(start+1,CurrentRightProb,CurrentNoRecoRightProb,
                      tHap.getScaffoldedHaplotype(hapID,start+1),
                      Error[start+1],
                       tHap.getScaffoldedHaplotype(hapID,start+1)?
                          alleleFreq[start+1] : 1-alleleFreq[start+1],Info);
        }

    tempRightProb=CurrentRightProb;
    Transpose(tempRightProb,CurrentRightProb,CurrentNoRecoRightProb,Recom[start],Info.uniqueCardinality);


    if(start==0)
        Impute(start,tHap.getScaffoldedHaplotype(hapID,start),tHap.getMissingScaffoldedHaplotype(hapID,start),Leftprob[0],
               CurrentRightProb,leftNoRecomProb[0],
                       CurrentNoRecoRightProb,Leftprob[0],PrevRightFoldedProb,Constants,Info,alleleFreq);

}





void MarkovModel::CreatePosteriorProb(vector<float> &Leftprob,vector<float> &rightProb,
                         vector<float> &leftNoRecoProb,vector<float> &rightNoRecoProb,
                         vector<float> &leftEndProb,vector<float> &rightEndProb,
                         vector<float> &Constants,vector<float> &probHap,ReducedHaplotypeInfo &Info)
{

    double value=0.0;
    for(int i=0;i<noReducedStatesCurrent;i++)
    {

        value = Constants[i]*(leftNoRecoProb[i]*rightNoRecoProb[i]/(leftEndProb[i]*rightEndProb[i]))
            +(Leftprob[i]*rightProb[i]-leftNoRecoProb[i]*rightNoRecoProb[i])*(Info.InvuniqueCardinality[i]);

        probHap[i]=value;

    }
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
    bool mle = false;

    if(Pref<Palt)
    {
        mle=true;
    }



    imputedDose[position] += imputedHap[position]= (Palt / ptotal);
    imputedAlleleNumber[position] = mle;

    if(!observedMiss)
    {

        double fmatch = 1.0 / (1. - Error[position] + Error[position] * ( major[position]? alleleFreq[position] : 1-alleleFreq[position]  ) + backgroundError );
        double fmismatch = 1.0 / (Error[position] * ( major[position]? alleleFreq[position] : 1-alleleFreq[position]  ) + backgroundError);

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

        ptotal =Pref+Palt;
        leaveOneOut[position] = Palt / ptotal;
    }


}




void MarkovModel::foldProbabilities(vector<float> &foldProb,int bridgeIndex,ReducedHaplotypeInfo &Info,int direction,int noReference) //0 - left; 1 - right
{
    vector<int> *TempuniqueIndexMap=&Info.uniqueIndexMap;
    fill(foldProb.begin(), foldProb.end(), 0.0);
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

        for(int i=0;i<noReference;i++)
        {
            foldProb[(*TempuniqueIndexMap)[i]]+=PrevjunctionRightProb[i];
        }
    }
}




void MarkovModel::unfoldProbabilities(int bridgeIndex,vector<float> &recomProb,
                                       vector<float> &noRecomProb,vector<float> &PrevFoldedProb,
                                     int direction,vector<ReducedHaplotypeInfo> &StructureInfo,
                                     int noReference)
{
    ReducedHaplotypeInfo &thisInfo = StructureInfo[bridgeIndex];
    int N = thisInfo.RepSize;

    float *adj_rec = (float *)alloca(N*sizeof(float));
    float *adj_norec = (float *)alloca(N*sizeof(float));

    for (int i=0; i<N; i++)
    {
        adj_rec[i] = recomProb[i] * thisInfo.InvuniqueCardinality[i];
    }
    for (int i=0; i<N; i++)
    {
        adj_norec[i] = noRecomProb[i] / PrevFoldedProb[i];
    }
    vector<float> &prev = direction ? PrevjunctionRightProb : junctionLeftProb[bridgeIndex];
    vector<float> &next = direction ? PrevjunctionRightProb : junctionLeftProb[bridgeIndex+1];

    if(direction)
    {
        for (int i=0; i<noReference; i++)
        {
            int m = thisInfo.uniqueIndexMap[i];
            prev[i]*=adj_norec[m];
            prev[i]+=adj_rec[m];
        }
    }
    else
    {
         for (int i=0; i<noReference; i++)
        {
            int m = thisInfo.uniqueIndexMap[i];
            next[i] = adj_rec[m] + adj_norec[m]*prev[i];
        }
    }

}



void MarkovModel::WalkLeft(HaplotypeSet &tHap, int &hapID,
                           int group, ReducedHaplotypeInfo &Info,
                           vector<double> &alleleFreq)
{
    vector<vector<float> > &Leftprob = LowMemory ? ThisBlockLeftProb : leftProb[group];
    if(LowMemory)
        junctionLeftFoldProb[group]=Leftprob[0];

    int &Start=Info.startIndex;
    int &End=Info.endIndex;

    noReducedStatesCurrent=Info.RepSize;

    for (int markerPos=Start+1; markerPos<=End; markerPos++)
    {

        PrecisionJump[markerPos]=Transpose(Leftprob[markerPos-Start-1],
                  Leftprob[markerPos-Start],CurrentLeftNoRecoProb,
                  Recom[markerPos-1],Info.uniqueCardinality);

        if (!missing[markerPos] && !tHap.getMissingScaffoldedHaplotype(hapID,markerPos))
        {
           if(!tHap.getMissingScaffoldedHaplotype(hapID,markerPos))
                Condition(markerPos,Leftprob[markerPos-Start],
                     CurrentLeftNoRecoProb,
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




void MarkovModel::CheckSize(HaplotypeSet & rHap,HaplotypeSet &tHap)
{

    cout.precision(4);
    double LeftBlockSize=0.0,LeftNoRecoSum = 0.0, TransSize=0.0, HapData=0.0;


    for(int i=0;i<rHap.NoBlocks;i++)
    {
        ReducedHaplotypeInfo &TempBlock=rHap.ReducedStructureInfo[i];


        HapData+=TempBlock.uniqueIndexMap.capacity()*sizeof(int);
        HapData+=TempBlock.uniqueCardinality.capacity()*sizeof(int);
        HapData+=TempBlock.InvuniqueCardinality.capacity()*sizeof(float);

        for(int j=0;j<TempBlock.RepSize;j++)
        {
            HapData+=TempBlock.uniqueHaps[j].capacity()*sizeof(bool);
        }
    }

    cout<<" HAP_DATA = "<<HapData/(1024*1024*1024)<<" ";




    if(LowMemory)
    {

        for(int i=0;i<rHap.NoBlocks;i++)
        {
            LeftBlockSize+=junctionLeftFoldProb[i].capacity()*sizeof(float);
        }


        cout<<" R_LEFT_BLOCK = "<<LeftBlockSize/(1024*1024*1024)<<" ";

        double temp=0.0;


        for(int i=0;i<rHap.maxBlockSize;i++)
        {
            temp+=ThisBlockLeftProb[i].capacity()*sizeof(float);
        }
        cout<<" R_LEFT_SAVE = "<<temp/(1024*1024*1024)<<" ";
        LeftBlockSize+=temp;
    }
    else
    {
        for(int i=0;i<rHap.NoBlocks;i++)
        {
            ReducedHaplotypeInfo &TempBlock=rHap.ReducedStructureInfo[i];
            vector<vector<float> > &TempLeft=leftProb[i];
            for(int j=0;j<TempBlock.BlockSize;j++)
            {
//                cout<<TempLeft[j].capacity()<<"\t"<<sizeof(float)<<endl;
                LeftBlockSize+=TempLeft[j].capacity()*sizeof(float);
            }
        }
    }



    cout<<" R_LEFT = "<<LeftBlockSize/(1024*1024*1024)<<" ";
    LeftNoRecoSum+=CurrentLeftNoRecoProb.capacity()*sizeof(float);
    for(int i=0;i<rHap.maxBlockSize;i++)
    {
        LeftNoRecoSum+=ThisBlockLeftNoRecoProb[i].capacity()*sizeof(float);
    }
    cout<<" NOR_LEFT = "<<LeftNoRecoSum/(1024*1024*1024)<<" ";
    cout<<" R_RIGHT = "<<tempRightProb.capacity()*sizeof(float)/(1024*1024*1024)<<" ";
    cout<<" NOR_RIGHT = "<<tempRightProb.capacity()*sizeof(float)/(1024*1024*1024)<<" ";

    LeftBlockSize+=LeftNoRecoSum;
    LeftBlockSize+=(2*tempRightProb.capacity()*sizeof(float));

    cout<<" BLOCK_SUM = "<<LeftBlockSize/(1024*1024*1024)<<" ";



    for(int i=0;i<=rHap.NoBlocks;i++)
    {
        TransSize+=junctionLeftProb[i].capacity()*sizeof(float);
    }

    cout<<" J_LEFT = "<<TransSize/(1024*1024*1024)<<" ";
    cout<<" J_RIGHT = "<<PrevjunctionRightProb.capacity()*sizeof(float)/(1024*1024*1024)<<" ";
    TransSize+=PrevjunctionRightProb.capacity()*sizeof(float);
    cout<<" J_SUM = "<<TransSize/(1024*1024*1024)<<" ";

    cout<<endl;


}




