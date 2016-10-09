#include "HaplotypeSet.h"

void printErr(String filename)
{
    cout<<" Error in M3VCF File !!! "<<endl;
    cout<<" Please re-construct the following [.m3vcf] file using Minimac2 and try again ..."<<endl;
    cout<< " [ "<< filename<<" ] "<<endl;
    cout<<" Contact author if problem still persists : sayantan@umich.edu "<<endl;
    cout<<" Program Exiting ..."<<endl<<endl;
    abort();
}

void HaplotypeSet::CreateSummary()
{
    CalculateFreq();

    NoBlocks=(int)ReducedStructureInfo.size();
    maxBlockSize=0;
    maxRepSize=0;


    for(int i=0;i<NoBlocks;i++)
    {
        ReducedHaplotypeInfo &TempBlock=ReducedStructureInfo[i];

        TempBlock.BlockSize=TempBlock.endIndex - TempBlock.startIndex +1;
        TempBlock.RepSize=(int)TempBlock.uniqueCardinality.size();

        if(maxBlockSize<TempBlock.BlockSize)
            maxBlockSize=TempBlock.BlockSize;
        if(maxRepSize<TempBlock.RepSize)
            maxRepSize=TempBlock.RepSize;

        TempBlock.InvuniqueCardinality.resize(TempBlock.RepSize);

        for (int j = 0; j < TempBlock.RepSize; j++)
        {
            TempBlock.InvuniqueCardinality[j]=1.0/(float)TempBlock.uniqueCardinality[j];
        }
    }

}


void HaplotypeSet::PrintDosageForVcfOutputForIDMaleSamples(IFILE vcfdose,
                                                           int MarkerIndex,bool majorIsReference,char refAllele)
{

    bool colonIndex;
    int length=(int)Dosage.size();

    for(int hapId=0;hapId<length;hapId++)
        {
            bool a1=ImputedAlleles[hapId][MarkerIndex];
            colonIndex=false;
            ifprintf(vcfdose,"\t");
            if(GT)
            {
                int outAllele1=0;
                if(a1)
                    outAllele1=1;
                ifprintf(vcfdose,"%d",outAllele1);
                colonIndex=true;
            }

            double x=Dosage[hapId][MarkerIndex];

            if(DS)
            {
                if(colonIndex)
                    ifprintf(vcfdose,":");

                ifprintf(vcfdose,"%.3f",x);
                colonIndex=true;
            }
            if(GP)
            {
                if(colonIndex)
                    ifprintf(vcfdose,":");
                double p1,p3;
                p3=x;
                p1=(1-x);
                ifprintf(vcfdose,"%.3f,%.3f",p1,p3);
            }
        }
}


void HaplotypeSet::PrintDosageForVcfOutputForID(IFILE vcfdose, int MarkerIndex,bool majorIsReference,char refAllele)
{

    bool colonIndex;
    int length=(int)Dosage.size()/2;

    for(int hapId=0;hapId<length;hapId++)
        {
            bool a1=ImputedAlleles[2*hapId][MarkerIndex];
            bool a2=ImputedAlleles[2*hapId+1][MarkerIndex];
            colonIndex=false;
            ifprintf(vcfdose,"\t");
            if(GT)
            {
                int outAllele1=0,outAllele2=0;
                if(a1)
                    outAllele1=1;

                if(a2)
                    outAllele2=1;

                if(!unphasedOutput)
                    ifprintf(vcfdose,"%d|%d",outAllele1,outAllele2);
                else
                {
                    if((a1^a2)==1)
                        ifprintf(vcfdose,"0/1");
                    else if(a1 && a2)
                        ifprintf(vcfdose,"1/1");
                    else
                        ifprintf(vcfdose,"0/0");
                }
                colonIndex=true;
            }

            double x=Dosage[2*hapId][MarkerIndex];
            double y=Dosage[2*hapId+1][MarkerIndex];

            if(DS)
            {

                if(colonIndex)
                    ifprintf(vcfdose,":");
                ifprintf(vcfdose,"%.3f",x+ y);
                colonIndex=true;
            }

            if(GP)
            {

                if(colonIndex)
                    ifprintf(vcfdose,":");
                colonIndex=false;
                double p1,p2,p3;
                p3=x*y;
                p2=x*(1-y)+y*(1-x);
                p1=(1-x)*(1-y);
                ifprintf(vcfdose,"%.3f,%.3f,%.3f",p1,p2,p3);
            }
        }

}


void HaplotypeSet::PrintDosageGWASOnlyForVcfOutputForIDMaleSamples(HaplotypeSet &tHap,
                                                        IFILE vcfdose,
                                                        int MarkerIndex)
{

    bool colonIndex;
    int length=(int)GWASOnlySampleIndex.size();

    for(int index=0;index<length;index++)
    {

        int hapId=GWASOnlySampleIndex[index];
        bool a1=tHap.GWASOnlyMissingSampleUnscaffolded[hapId][MarkerIndex];
        ifprintf(vcfdose,"\t");

        double outAllele1=0.0;
        if(a1)
        {
            outAllele1=tHap.AlleleFreq[MarkerIndex];
            a1=round(outAllele1)==1?true:false;
        }
        else
        {
            a1=tHap.GWASOnlyhaplotypesUnscaffolded[hapId][MarkerIndex];
            if(a1)
                outAllele1=1.0;
        }
        colonIndex=false;
        if(GT)
        {

            ifprintf(vcfdose,"%d",(int)a1);
            colonIndex=true;
        }
        if(DS)
        {

            if(colonIndex)
                ifprintf(vcfdose,":");
            ifprintf(vcfdose,"%.3f",(float)outAllele1);
            colonIndex=true;
        }
        if(GP)
        {

            if(colonIndex)
                ifprintf(vcfdose,":");
            double p1,p3;

            p3=outAllele1;
            p1=(1-outAllele1);

            ifprintf(vcfdose,"%.3f,%.3f",p1,p3);
        }
    }
}


void HaplotypeSet::PrintDosageGWASOnlyForVcfOutputForID(HaplotypeSet &tHap,
                                                        IFILE vcfdose,
                                                        int MarkerIndex)
{

    bool colonIndex;
    int length=(int)GWASOnlySampleIndex.size();

    for(int index=0;index<length;index++)
    {

        int hapId=GWASOnlySampleIndex[index];

        bool a1=tHap.GWASOnlyMissingSampleUnscaffolded[2*hapId][MarkerIndex];
        bool a2=tHap.GWASOnlyMissingSampleUnscaffolded[2*hapId+1][MarkerIndex];

        ifprintf(vcfdose,"\t");

        double outAllele1=0.0,outAllele2=0.0;

        if(a1 || a2)
        {
            outAllele1=tHap.AlleleFreq[MarkerIndex];
            outAllele2=outAllele1;
            a1=round(outAllele1)==1?true:false;
            a2=a1;
        }
        else
        {
            a1=tHap.GWASOnlyhaplotypesUnscaffolded[2*hapId][MarkerIndex];
            a2=tHap.GWASOnlyhaplotypesUnscaffolded[2*hapId+1][MarkerIndex];
            if(a1)
                outAllele1=1.0;

            if(a2)
                outAllele2=1.0;
        }

        colonIndex=false;

        if(GT)
        {

            if(!unphasedOutput)
                ifprintf(vcfdose,"%d|%d",(int)a1,(int)a2);
            else
            {
                if((a1^a2)==1)
                    ifprintf(vcfdose,"0/1");
                else if(a1 && a2)
                    ifprintf(vcfdose,"1/1");
                else
                    ifprintf(vcfdose,"0/0");
            }
            colonIndex=true;
        }
        if(DS)
        {
            if(colonIndex)
                ifprintf(vcfdose,":");
            ifprintf(vcfdose,"%.3f",(float)outAllele1+ (float)outAllele2);
            colonIndex=true;
        }
        if(GP)
        {

            if(colonIndex)
                ifprintf(vcfdose,":");
            double p1,p2,p3;

            p3=outAllele1*outAllele2;
            p2=outAllele1*(1-outAllele2)+outAllele2*(1-outAllele1);
            p1=(1-outAllele1)*(1-outAllele2);

            ifprintf(vcfdose,"%.3f,%.3f,%.3f",p1,p2,p3);
        }
    }
}


void HaplotypeSet::InitializePartialDosageForVcfOutputMaleSamples(int NHaps,
                                                       int NMarkers,vector<bool> &Format)
{
    numHaplotypes = NHaps;
    numMarkers = NMarkers;
    Dosage.resize(numHaplotypes);

    if(TypedOnly)
    {
        GWASOnlySampleIndex.resize(numHaplotypes);
    }
    ImputedAlleles.resize(numHaplotypes);
    individualName.resize(numHaplotypes);
    for(int i=0;i<numHaplotypes;i++)
        {
            Dosage[i].resize(NMarkers,3.0);
            ImputedAlleles[i].resize(NMarkers,false);
        }
    GT=Format[0];
    DS=Format[1];
    GP=Format[2];
}

void HaplotypeSet::InitializePartialDosageForVcfOutput(int NHaps,
                                                       int NMarkers,vector<bool> &Format)
{
    numHaplotypes = NHaps;
    numMarkers = NMarkers;
    Dosage.resize(numHaplotypes);

    if(TypedOnly)
    {
        GWASOnlySampleIndex.resize(numHaplotypes/2);
    }

    ImputedAlleles.resize(numHaplotypes);
    individualName.resize(numHaplotypes/2);
    for(int i=0;i<numHaplotypes;i++)
        {
            Dosage[i].resize(NMarkers,3.0);
            ImputedAlleles[i].resize(NMarkers,false);
        }

    GT=Format[0];
    DS=Format[1];
    GP=Format[2];

}


void HaplotypeSet::SaveIndexForGWASOnlyForVcfOutput(int SamID,int HapId)

{
    GWASOnlySampleIndex[SamID]=HapId;
}

void HaplotypeSet::SaveDosageForVcfOutputSampleWise(int SamID,string &SampleName, vector<float> &dose1,vector<float> &dose2,
                                                    vector<bool> &impAlleles1,vector<bool> &impAlleles2)
{
    individualName[SamID]=SampleName;
    Dosage[2*SamID]=dose1;
    ImputedAlleles[2*SamID]=impAlleles1;
    Dosage[2*SamID+1]=dose2;
    ImputedAlleles[2*SamID+1]=impAlleles2;
}

void HaplotypeSet::SaveDosageForVcfOutputSampleWiseChrX(int SamID,string &SampleName, vector<float> &dose1,
                                                    vector<bool> &impAlleles1)
{
    individualName[SamID]=SampleName;
    Dosage[SamID]=dose1;
    ImputedAlleles[SamID]=impAlleles1;
}


void HaplotypeSet::reconstructHaplotype(vector<bool> &reHaplotypes,int &index)
{
    if((int)reHaplotypes.size()!=numMarkers)
    {
        cout<<" SIZE MISMATCH "<<endl;
        abort();
    }

    int markerIndex=0,k;
    bool checkAllele=false;

    for(int j=0;j<(int)ReducedStructureInfo.size();j++)
    {
        for(k=0;k<(ReducedStructureInfo[j].endIndex-ReducedStructureInfo[j].startIndex);k++)
        {
            reHaplotypes[markerIndex++]=ReducedStructureInfo[j].uniqueHaps[ReducedStructureInfo[j].uniqueIndexMap[index]][k];
            if(k==0 && j>1)
            {
                if(checkAllele!=(ReducedStructureInfo[j].uniqueHaps[ReducedStructureInfo[j].uniqueIndexMap[index]][k]))
                {
                    cout<<index<<"\t"<<j<<"\t"<<k<<"\t"<<(int)checkAllele<<"\t"<<(int)(ReducedStructureInfo[j].uniqueHaps[ReducedStructureInfo[j].uniqueIndexMap[index]][k])<<endl;
                    cout<<" CHECK ALLELEE  MISMATCH "<<endl;
                    abort();
                }

            }

        }

        if(j==((int)ReducedStructureInfo.size()-1))
        {

            reHaplotypes[markerIndex]=ReducedStructureInfo[j].uniqueHaps[ReducedStructureInfo[j].uniqueIndexMap[index]][k];
        }
        else
        {
            checkAllele=ReducedStructureInfo[j].uniqueHaps[ReducedStructureInfo[j].uniqueIndexMap[index]][k];
        }

    }
}


bool HaplotypeSet::readm3vcfFile(String m3vcfFile,String CHR,int START,int END,int WINDOW)
{
    string line,tempString,tempBlockPos,tempString2,tempString3,tempName,tempRsId,tempChr;
    variant tempVariant;
    variant tempVariant2;
    int InitialNMarkers=0,blockIndex,startIndexFlag=0,readIndex=0,writeBlockIndex=0,NoBlocks=0,tempPos,NoMarkersWritten=0,tempVarCount,tempRepCount;

    cout<<"\n Reading Reference Haplotype information from M3VCF files : "<<m3vcfFile<<endl<<endl;
    OrigStartPos=START;
    OrigEndPos=END;
    if (CHR!="" && WINDOW > 0)
    {
        if (START-WINDOW < 0)
            START = 0;
        else
            START -= WINDOW;

        END += WINDOW;
    }
    if(CHR!="")
    {
        std::cout << "\n Loading markers from chromosome " << CHR;
        if(END>0)
            std::cout << " from base position "<<START<<" to base position "<<END<<"."<< endl;
        else
            std::cout << " from base position "<<START<<" till end of M3VCF file."<< endl;
    }


    PrintStartIndex=99999999;
    PrintEndIndex=-1;


    IFILE m3vcfxStream = ifopen(m3vcfFile, "r");

    if(m3vcfxStream)
    {


        m3vcfxStream->readLine(line);
        if(line.compare("##fileformat=M3VCF")!=0 && line.compare("##fileformat=OPTM")!=0)
        {
            cout<<" Incorrect Header Information : "<<line<<endl;
            cout<<" Header line should be : ##fileformat=M3VCF "<<endl<<endl;
            printErr(m3vcfFile);
        }

        bool Header=true;
//        char * pch_split,* pch_split2,* pch_split3;
        char * pch;
        char *end_str1,*end_str2;

        while(Header)
        {
            line.clear();
            m3vcfxStream->readLine(line);
            if(line.substr(0,2).compare("##")==0)
                Header=true;
            else
                break;
            tempString = (string) strtok_r ((char*)line.c_str(),"=", &end_str1);
            pch = strtok_r (NULL, "=", &end_str1);

            if(tempString.compare("##n_blocks")==0)
            {
                NoBlocks=atoi(pch);
                continue;
            }
            else if(tempString.compare("##n_haps")==0)
            {
                numHaplotypes=atoi(pch);
                continue;
            }
            else if(tempString.compare("##n_markers")==0)
            {
                InitialNMarkers=atoi(pch);
//                abort();
                continue;

            }
            else if(tempString.compare("##chrxRegion")==0)
            {
                tempString = (string) pch;
                if(tempString.compare("NonPseudoAutosomal")==0)
                    PseudoAutosomal=false;
                else if(tempString.compare("PseudoAutosomal")==0)
                    PseudoAutosomal=true;
                else
                {
                    cout << "\n Inconsistent Tag for Chr X. "<<endl;
                    cout << " Please check the file properly..\n";
                    cout << " Program Aborting ... "<<endl;
                    return false;

                }
                continue;

            }

        }

        int colCount=0;
        pch = strtok_r ((char*)line.c_str(),"\t", &end_str2);


        while(pch!=NULL)
        {
            colCount++;
            if(colCount>9)
            {
                tempString2=string(pch);
                individualName.push_back(tempString2.substr(0,tempString2.size()-6));
            }
            pch = strtok_r (NULL,"\t", &end_str2);
        }


        cout<<" Reading  "<<numHaplotypes<< " haplotypes from data ..."<<endl<<endl;

        if((int)individualName.size()!=numHaplotypes)
        {
            cout<<endl<<" Error in Data consistency !!! "<<endl<<endl;
            cout<<" Number of Haplotypes should be : "<< numHaplotypes<<", but "<<individualName.size()<<" Haplotypes found in header row !!!"<< endl<<endl;
            printErr(m3vcfFile);
        }


        if(individualName.size()==0)
        {
            cout << "\n No haplotypes recorded from M3VCF Input File : "<<m3vcfFile<<endl;
            cout << " Please check the file properly..\n";
            cout << " Program Aborting ... "<<endl;
            return false;
        }



        for(blockIndex=0;blockIndex<NoBlocks;blockIndex++)
        {


            if (blockIndex % 1000 == 0)
			{
			    printf("  Loading Block %d out of %d blocks to be loaded... [%.1f%%] "
                , blockIndex + 1, NoBlocks, 100*(double)(blockIndex + 1)/(double)NoBlocks);
                cout<<endl;
            }

            int flag=0,blockEnterFlag=0,blocktoSave=0;
            line.clear();
            m3vcfxStream->readLine(line);
            char *end_str_new;

            pch = strtok_r ((char*)line.c_str(),"\t", &end_str_new);
            tempChr=string(pch);

            pch = strtok_r (NULL,"\t", &end_str_new);
            tempBlockPos=string(pch);

            char *end_str_new1;
            char *pch1;

            pch1 = strtok_r ((char*)tempBlockPos.c_str(),"-", &end_str_new1);
            int tempStartBlock=atoi(pch1);

            pch1 = strtok_r (NULL,"-", &end_str_new1);
            int tempEndBlock=atoi(pch1);

            if(CHR!="")
            {
                if(tempChr.compare(CHR.c_str())!=0)
                    flag=1;
                else
                {
                    if(END>0)
                    {
                        if(tempStartBlock>END)
                            flag=1;
                    }
                    if(tempEndBlock<START)
                        flag=1;
                }
            }


            pch = strtok_r (NULL,"\t", &end_str_new);
            string blockName=string(pch);

            pch = strtok_r (NULL,"\t", &end_str_new);
            pch = strtok_r (NULL,"\t", &end_str_new);
            pch = strtok_r (NULL,"\t", &end_str_new);
            pch = strtok_r (NULL,"\t", &end_str_new);

            pch = strtok_r (NULL,"\t", &end_str_new);
            tempString2=string(pch);

            char *pch2,*end_str_new2;
            pch2 = strtok_r ((char*)tempString2.c_str(),";", &end_str_new2);

            stringstream strs;
            strs<<(blockIndex+1);

            tempString3="B"+ (string)(strs.str());

            if(tempString3.compare((string)pch2)!=0)
            {
                cout<<endl<<" Error in INFO column (Block Identifier) for block : "<<blockName <<endl;
                cout<<" Block Identifier should be : "<< tempString3<<" but is : "<<pch2<<endl<<endl;
                printErr(m3vcfFile);
            }



            tempString3 = (string)strtok_r (NULL,";", &end_str_new2);
            char *pch3,*end_str_new3;
            pch3 = strtok_r ((char*)tempString3.c_str(),"=", &end_str_new3);
            pch3 = strtok_r (NULL,"=", &end_str_new3);
            tempVarCount=atoi(pch3);

            tempString3 = (string)strtok_r (NULL,";", &end_str_new2);
            end_str_new3=NULL;
            pch3 = strtok_r ((char*)tempString3.c_str(),"=", &end_str_new3);
            pch3 = strtok_r (NULL,"=", &end_str_new3);
            tempRepCount=atoi(pch3);

            ReducedHaplotypeInfo tempBlock;

            if(flag==1)
            {
                for(int tempIndex=0;tempIndex<tempVarCount;tempIndex++)
                m3vcfxStream->discardLine();
                    continue;
            }


            tempBlock.uniqueCardinality.resize(tempRepCount,0.0);
            tempBlock.InvuniqueCardinality.resize(tempRepCount,0.0);


            tempBlock.uniqueHaps.resize(tempRepCount);
            tempBlock.uniqueIndexMap.resize(numHaplotypes);

            pch = strtok_r (NULL,"\t", &end_str_new);
            pch = strtok_r (NULL,"\t", &end_str_new);



            int check=0;
            while(pch!=NULL)
            {
                int tempval=atoi(pch);
                tempBlock.uniqueIndexMap[check]=tempval;
                tempBlock.uniqueCardinality[tempval]++;

                check++;
                pch = strtok_r (NULL,"\t", &end_str_new);
            }

            for (int i = 0; i < tempRepCount; i++)
            {
                tempBlock.InvuniqueCardinality[i]=1.0/(float)tempBlock.uniqueCardinality[i];
            }
            if(check!=numHaplotypes)
            {
                cout<<endl<<" Error in Data consistency !!! "<<endl;
                cout<<" Number of Haplotypes should be : "<< numHaplotypes<<", but "<<check<<" indices found in block"<<blockName << endl<<endl;
                printErr(m3vcfFile);
            }
            for(int tempIndex=0;tempIndex<tempVarCount;tempIndex++)
            {
                flag=0;
                line.clear();
                m3vcfxStream->readLine(line);

                end_str_new3=NULL;


                pch3 = strtok_r ((char*)line.c_str(),"\t", &end_str_new3);
                tempChr=(string)pch3;

                pch3 = strtok_r (NULL,"\t", &end_str_new3);
                tempPos=atoi(pch3);

                stringstream strs3;
                strs3<<tempPos;


                tempName=tempChr+":"+(string)strs3.str();

                pch3 = strtok_r (NULL,"\t", &end_str_new3);
                tempRsId=(string)(pch3);

                if(CHR!="")
                {
                    if(tempChr.compare(CHR.c_str())!=0)
                        flag=1;
                    else
                    {
                        if(END>0)
                        {
                            if(tempPos>END || tempPos<START)
                                flag=1;
                        }
                        else
                            if(tempPos<START)
                                flag=1;
                         if(tempIndex==(tempVarCount-1) && writeBlockIndex==0)
                            flag=1;
                    }

                }
                readIndex++;

                if(flag==1)
                    continue;
                else
                {
                    if(blockEnterFlag==0)
                        tempBlock.startIndex=writeBlockIndex;
                    else
                        tempBlock.endIndex=writeBlockIndex;
                    if(tempIndex<(tempVarCount-1) && tempIndex>0) // to ensure a single marker from a block is NOT read.
                        blocktoSave=1;
                    blockEnterFlag=1;

                }

                if(tempPos>=OrigStartPos && startIndexFlag==0)
                    {
                        PrintStartIndex=writeBlockIndex;
                        startIndexFlag=1;
                    }

                if(CHR=="")
                {
                    PrintEndIndex=writeBlockIndex;
                }
                else
                {
                    if(OrigEndPos==0 || tempPos<=OrigEndPos)
                        PrintEndIndex=writeBlockIndex;
                }


                variant tempVariant;

                tempVariant.assignValues(tempName,tempChr,tempPos);
                tempVariant.rsid=tempRsId;

                double tempRecom=-3.0,tempError=0.0;

                pch3 = strtok_r (NULL,"\t", &end_str_new3);
                string tempString98=(string)(pch3);

                pch3 = strtok_r (NULL,"\t", &end_str_new3);
                string tempString99=(string)(pch3);

                tempVariant.assignRefAlt(tempString98,tempString99);

                pch3 = strtok_r (NULL,"\t", &end_str_new3);
                pch3 = strtok_r (NULL,"\t", &end_str_new3);
                pch3 = strtok_r (NULL,"\t", &end_str_new3);

                tempString=(string)(pch3);

                char *pch4, *end_str_new4;

                pch4=strtok_r ((char*)tempString.c_str(),";", &end_str_new4);
                stringstream strs1,strs2;
                strs1<<(blockIndex+1);
                strs2<<(tempIndex+1);
                tempString3="B"+ (string)(strs1.str())+ ".M"+ (string)(strs2.str());
                if(tempString3.compare((string)pch4)!=0)
                {
                    cout<<endl<<" Error in INFO column (Block Identifier) for variant : "<<tempName <<endl;
                    cout<<" Block Identifier should be : "<< tempString3<<" but is : "<<pch4<<endl<<endl;
                    printErr(m3vcfFile);
                }


                pch4=strtok_r (NULL,";", &end_str_new4);

                while(pch4!=NULL)
                {
                    tempString3=(string)(pch4);
                    char *pch5,*pch6,*end_str_new5;
                    pch5=strtok_r ((char*)tempString3.c_str(),"=", &end_str_new5);
                    pch6=strtok_r (NULL,"=", &end_str_new5);

                    if((string)pch5=="Err")
                    {
                        tempError=atof(pch6);
                    }
                    if((string)pch5=="Recom")
                    {
                        tempRecom=atof(pch6);
                    }
                    pch4=strtok_r (NULL,";", &end_str_new4);
                }

                char Rallele=0,Aallele=0;


                string refAllele=tempVariant.refAlleleString;
                string altAllele=tempVariant.altAlleleString;

                if (strlen(refAllele.c_str()) == 1
                     && strlen(altAllele.c_str()) == 1)
                {
                    switch (refAllele[0])
                    {
                        case 'A': case 'a': Rallele = 1; break;
                        case 'C': case 'c': Rallele = 2; break;
                        case 'G': case 'g': Rallele = 3; break;
                        case 'T': case 't': Rallele = 4; break;
                        case 'D': case 'd': Rallele = 5; break;
                        case 'I': case 'i': Rallele = 6; break;
                        case 'R': case 'r': Rallele = 7; break;
                        default:
                        {
                                   cout << "\n\n Data Inconsistency !!! \n";
                                   cout << " Error with reference allele for marker : " << tempVariant.rsid<< " in M3VCF File : " << m3vcfFile;
                                   cout << "\n VCF reference alleles for SNPs can only be A(a), C(c), G(g), or T(t).\n";
                                   cout << " " << tempVariant.rsid << " has " << refAllele << endl;
                                   cout << "\n Program Aborting ... \n\n";
                                   abort();
//                                       return false;
                        }
                    }

                    switch (altAllele[0])
                    {
                        case 'A': case 'a': Aallele = 1; break;
                        case 'C': case 'c': Aallele = 2; break;
                        case 'G': case 'g': Aallele = 3; break;
                        case 'T': case 't': Aallele = 4; break;
                        case 'D': case 'd': Aallele = 5; break;
                        case 'I': case 'i': Aallele = 6; break;
                        case 'R': case 'r': Aallele = 7; break;
                        default:
                        {
                                   cout << "\n\n Data Inconsistency !!! \n";
                                   cout << " Error with alternate allele for marker : " <<tempVariant.rsid<< " in M3VCF File : " << m3vcfFile;
                                   cout << "\n VCF alternate alleles for SNPs can only be A(a), C(c), G(g), or T(t).\n";
                                   cout << " " << tempVariant.rsid << " has " << altAllele << endl;
                                   cout << "\n Program Aborting ... \n\n";
                                   abort();
                        }
                    }
                }
                else
                {
                    Rallele = 7;
                    if(strlen(refAllele.c_str())<strlen(altAllele.c_str()))
                        Aallele=6;
                    else
                        Aallele=5;
                }
                tempVariant.refAllele=Rallele;
                tempVariant.altAllele=Aallele;



                if(tempIndex<(tempVarCount-1) || blockIndex==(NoBlocks-1))
                {

                    VariantList.push_back(tempVariant);
                    refAlleleList.push_back(tempVariant.refAllele);
                    markerName.push_back(tempName);
                    if(tempRecom!=-3.0)
                    {
                        Recom.push_back(tempRecom);
                        Error.push_back(tempError);
                    }
                    writeBlockIndex++;
                }

                pch3 = strtok_r (NULL,"\t", &end_str_new3);
                tempString=(string)(pch3);

                for(check=0;check<tempRepCount;check++)
                {

                    char t=tempString[check];

                    tempBlock.uniqueHaps[check].push_back(t=='0'? false:true);

                }

            }

            NoMarkersWritten+=(tempVarCount-1);
            if(blocktoSave==1)
                {
                optEndPoints.push_back(tempBlock.startIndex);

                ReducedStructureInfo.push_back(tempBlock);
                }
        }
        if(ReducedStructureInfo.size()>0)
            optEndPoints.push_back(ReducedStructureInfo[ReducedStructureInfo.size()-1].endIndex);

        if(Recom.size()>0)
            Recom.pop_back();

        finChromosome=tempChr;

    }
    else
    {
        cout<<" Following M3VCF File Not Available : "<<m3vcfFile<<endl;
        cout<<" Program Exiting ... "<<endl<<endl;
        return false;
    }


    numMarkers=writeBlockIndex;

    cout<<endl<<" Reference Haplotype information successfully recorded. "<<endl;


    if(finChromosome=="X")
    {
        cout<<"\n Chromosome X Detected !!! \n";
    }

	std::cout << "\n Number of Markers in File                           : " << InitialNMarkers << endl;
	std::cout << "\n Number of Markers Recorded                          : " << numMarkers << endl;
    std::cout << " Number of Haplotypes Recorded                       : " << numHaplotypes << endl;

    ifclose(m3vcfxStream);

   if(numMarkers<2)
    {
        cout << "\n None/Single marker left after filtering from Input File : "<<m3vcfFile<<endl;
		cout << " Please check the file or the filtering options properly ...\n";
		cout << " Program Aborting ... "<<endl;
		return false;
    }

    std::cout << "\n Haplotype Set successfully loaded from M3VCF File     : " << m3vcfFile << endl;

    return true;
}


bool HaplotypeSet::CheckValidChrom(string chr)
{
    bool result=false;

    if(MyChromosome!="" && chr==MyChromosome)
        return true;

    string temp[]={"1","2","3","4","5","6","7","8","9","10","11"
                    ,"12","13","14","15","16","17","18","19","20","21","22","X"};
    std::vector<string> ValidChromList (temp, temp + sizeof(temp) / sizeof(string) );

    for(int counter=0;counter<(int)ValidChromList.size();counter++)
        if(chr==ValidChromList[counter])
            result=true;

    return result;

}

void HaplotypeSet::writem3vcfFile(String filename,bool &gzip)
{


    IFILE m3vcffile = ifopen(filename + ".m3vcf" + (gzip ? ".gz" : ""), "wb",(gzip ? InputFile::BGZF : InputFile::UNCOMPRESSED));
    ifprintf(m3vcffile, "##fileformat=M3VCF\n");
    ifprintf(m3vcffile, "##version=1.2\n");
    ifprintf(m3vcffile, "##compression=block\n");
    ifprintf(m3vcffile, "##n_blocks=%d\n",ReducedStructureInfo.size());
    ifprintf(m3vcffile, "##n_haps=%d\n",numHaplotypes);
    ifprintf(m3vcffile, "##n_markers=%d\n",numMarkers);
    if(finChromosome=="X")
        ifprintf(m3vcffile, "##chrxRegion=%s\n",PseudoAutosomal?"PseudoAutosomal":"NonPseudoAutosomal");
    ifprintf(m3vcffile, "##<Note=This is NOT a VCF File and cannot be read by vcftools>\n");
    ifprintf(m3vcffile, "#CHROM\t");
    ifprintf(m3vcffile, "POS\t");
    ifprintf(m3vcffile, "ID\t");
    ifprintf(m3vcffile, "REF\t");
    ifprintf(m3vcffile, "ALT\t");
    ifprintf(m3vcffile, "QUAL\t");
    ifprintf(m3vcffile, "FILTER\t");
    ifprintf(m3vcffile, "INFO\t");
    ifprintf(m3vcffile, "FORMAT");
    int i,j,k;

    for(i=0;i<(int)individualName.size();i++)
    {
        ifprintf(m3vcffile, "\t%s_HAP_1",individualName[i].c_str());

        if(finChromosome!="X")
            ifprintf(m3vcffile, "\t%s_HAP_2",individualName[i].c_str());
        else if(SampleNoHaplotypes[i]==2)
            ifprintf(m3vcffile, "\t%s_HAP_2",individualName[i].c_str());
    }
    ifprintf(m3vcffile, "\n");

    int length=(int)ReducedStructureInfo.size();
    string cno;

    for(i=0;i<length;i++)
    {

        cno=VariantList[ReducedStructureInfo[i].startIndex].chr;
        int nvariants=ReducedStructureInfo[i].endIndex-ReducedStructureInfo[i].startIndex+1;
        int reps=ReducedStructureInfo[i].uniqueCardinality.size();


        ifprintf(m3vcffile, "%s\t",cno.c_str());
        ifprintf(m3vcffile, "%d-%d\t",VariantList[ReducedStructureInfo[i].startIndex].bp,VariantList[ReducedStructureInfo[i].endIndex].bp);
        ifprintf(m3vcffile, "<BLOCK:%d-%d>\t.\t.\t.\t.\t",ReducedStructureInfo[i].startIndex,ReducedStructureInfo[i].endIndex);

        ifprintf(m3vcffile, "B%d;VARIANTS=%d;REPS=%d\t.",i+1,nvariants,reps);


        for(j=0;j<numHaplotypes;j++)
            ifprintf(m3vcffile, "\t%d",ReducedStructureInfo[i].uniqueIndexMap[j]);

        ifprintf(m3vcffile, "\n");

        for(j=0;j<nvariants;j++)
        {
            ifprintf(m3vcffile, "%s\t",cno.c_str());
            ifprintf(m3vcffile, "%d\t",VariantList[j+ReducedStructureInfo[i].startIndex].bp);
            ifprintf(m3vcffile, "%s\t",RsId?VariantList[j+ReducedStructureInfo[i].startIndex].rsid.c_str():VariantList[j+ReducedStructureInfo[i].startIndex].name.c_str());
            ifprintf(m3vcffile, "%s\t%s\t.\t.\t",VariantList[j+ReducedStructureInfo[i].startIndex].refAlleleString.c_str(),VariantList[j+ReducedStructureInfo[i].startIndex].altAlleleString.c_str());

            ifprintf(m3vcffile, "B%d.M%d",i+1,j+1);
            if(Error.size()>0)
                ifprintf(m3vcffile, ";Err=%.5g;Recom=%.5g",
                     Error[j+ReducedStructureInfo[i].startIndex],(j+ReducedStructureInfo[i].startIndex)<(int)Recom.size()?Recom[j+ReducedStructureInfo[i].startIndex]:0);
            ifprintf(m3vcffile, "\t");

            for(k=0;k<reps;k++)
            {
                ifprintf(m3vcffile,"%d",!ReducedStructureInfo[i].uniqueHaps[k][j]? 0:1);
            }
            ifprintf(m3vcffile, "\n");
        }
    }

    std::cout << " Successfully written file ... "<<endl;
    ifclose(m3vcffile);

}

string HaplotypeSet::DetectTargetFileType(String filename)
{
    IFILE fileStream = ifopen(filename, "r");

    string line;
    if(fileStream)
    {
        fileStream->readLine(line);
        if(line.length()<1)
            return "Invalid";
        string tempString;
        tempString=(line.substr(0,17));

        char temp[tempString.length() + 1];
        std::strcpy(temp,tempString.c_str());
        for (char *iter = temp; *iter != '\0'; ++iter)
        {
           *iter = std::tolower(*iter);
        }

        if(((string)temp).compare("##fileformat=vcfv")==0)
            {
                return "vcf";
            }
       return "Invalid";

    }
    else
    {
        return "NA";
    }

    ifclose(fileStream);

    return "NA";

}

string HaplotypeSet::DetectReferenceFileType(String filename)
{
    IFILE fileStream = ifopen(filename, "r");
    string line;
    if(fileStream)
    {
        fileStream->readLine(line);
        if(line.length()<1)
            return "Invalid";
        string tempString;
        tempString=(line.substr(0,17));

        char temp[tempString.length() + 1];
        std::strcpy(temp,tempString.c_str());
        for (char *iter = temp; *iter != '\0'; ++iter)
        {
           *iter = std::tolower(*iter);
        }
        if(((string)temp).compare("##fileformat=m3vc")==0)
        {
            return "m3vcf";
        }
        else if(((string)temp).compare("##fileformat=vcfv")==0)
        {
            return "vcf";
        }
        else
            return "Invalid";

    }
    else
    {
        return "NA";
    }

    ifclose(fileStream);

    return "NA";
}


bool HaplotypeSet::FasterLoadHaplotypes(String filename, int maxIndiv, int maxMarker,String CHR,
                                      int START,int END,int WINDOW,bool rsid,bool compressOnly,bool filter,String &outfile, bool &gz)
{

    RsId=rsid;
    outFile=outfile;
    gzip=gz;

    string FileType=DetectReferenceFileType(filename);
    Filter=filter;

    if(FileType.compare("m3vcf")==0)
    {
        cout<<"\n Format = M3VCF (Minimac3 VCF File) "<<endl;
        return readm3vcfFile(filename,CHR,START,END,WINDOW);
    }
    else if(FileType.compare("NA")==0)
    {
        cout<<"\n Following File File Not Available : "<<filename<<endl;
        return false;
    }
    else if(FileType.compare("Invalid")==0)
    {

        cout << "\n Reference File provided by \"--refHaps\" must be a VCF or M3VCF file !!! \n";
        cout << " Please check the following file : "<<filename<<endl;
        return false;
    }
    cout<<"\n Format = VCF (Variant Call Format) "<<endl;
    vcfType=true;

    optEndPoints.clear();
	VcfFileReader inFile;
	VcfHeader header;
	VcfRecord record;

	if (!inFile.open(filename, header))
	{
		cout << "\n Program could NOT open file : " << filename << endl;
		return false;
	}

    std::cout << "\n Loading Reference Haplotype Set from VCF File       : " << filename << endl;

    OrigStartPos=START;
    OrigEndPos=END;



    if (CHR!="" && WINDOW > 0)
    {
        if (START-WINDOW < 0)
            START = 0;
        else
            START -= WINDOW;

        END += WINDOW;
    }


    stringstream strs;
    strs<<(END);

    if(CHR!="")
    {
        std::cout << "\n Region specified by user (including window = "<<WINDOW <<" bp) : chr" << CHR
        <<":"<<START <<"-"<< (END > 0 ? (string)(strs.str()) :"END") << endl;
    }



    PrintStartIndex=0;
    PrintEndIndex=0;
    int     numReadRecords = 0;
	int numtoBeWrittenRecords = 0;
	vector<bool> importIndex;
    int numSamplesRead = 0;
	int notBiallelic = 0;
	int failFilter = 0;
	int duplicates = 0;
	int insertions = 0;
	int deletions = 0;
	int inconsistent=0;
	string prevID="";
//	bool isInDel = false;
	inFile.setSiteOnly(true);
	numSamplesRead = header.getNumSamples();

	if(numSamplesRead==0)
    {
        cout << "\n No samples are available in the file : "<<filename<<endl;
		cout << " Please check the input VCF file properly ...\n";
		cout << " Program Aborting ... "<<endl;
		return false;
    }


    individualName.clear();
	for (int i = 0; i < numSamplesRead; i++)
	{
		string tempName(header.getSampleName(i));
		individualName.push_back(tempName);
	}


    if(individualName.size()==0)
    {
        cout << "\n No haplotypes recorded from VCF Input File : "<<filename<<endl;
		cout << " Please check the file properly..\n";
		cout << " Program Aborting ... "<<endl;
		return false;
    }



    string refAllele,altAllele,PrefrefAllele,PrevaltAllele,cno,fixCno,id;
    int bp;
	cout << "\n Reading VCF File to calculate number of records ... \n";
    cout<<endl;
	// pre-calculate number of samples read and number of markers read to allocate memory.
	while (inFile.readRecord(record))
	{

        int flag = 0;
        if (maxMarker != 0 && numtoBeWrittenRecords>=maxMarker)
            break;

		++numReadRecords;

		if (record.getNumAlts()>1)
		{
			notBiallelic++;
			flag = 1;
		}
		if (record.getFilter().getString(0).compare("PASS") != 0)
		{
			failFilter++;
			if(Filter)
                flag = 1;
		}
//		isInDel=false;

		refAllele = record.getRefStr();
		cno=record.getChromStr();
        bp=record.get1BasedPosition();
        id=record.getIDStr();
        altAllele = record.getAltStr();


        if(numReadRecords==1 && CHR=="")
        {
            if(!CheckValidChrom(cno))
            {
                cout << "\n Error !!! Reference VCF File contains chromosome : "<<cno<<endl;
                cout << " VCF File can only contain chromosomes 1-22 and X !!! "<<endl;
                cout << " Program Aborting ... "<<endl;
                return false;
            }

        }
        if(CHR=="" && fixCno!=cno && numReadRecords>1)
        {
            cout << "\n Error !!! Reference VCF File contains multiple chromosomes : "<<cno<<", "<<fixCno<<", ... "<<endl;
            cout << " Please use VCF file with single chromosome or specify chromosome using \"--chr\" option !!! "<<endl;
            cout << " Program Aborting ... "<<endl;
            return false;
        }


        fixCno=cno;

        string currID;

        stringstream strs1,strs2;
        strs1<<(cno);
        strs2<<(bp);


        currID=(string)strs1.str()+":"+(string)strs2.str();

        if(id==".")
            id=currID;

        if (strlen(refAllele.c_str()) == 1 && strlen(altAllele.c_str()) == 1)
        {
            switch (refAllele[0])
            {
                case 'A': case 'a': break;
                case 'C': case 'c': break;
                case 'G': case 'g': break;
                case 'T': case 't': break;
                case 'D': case 'd': break;
                case 'I': case 'i': break;
                case 'R': case 'r': break;
                default:
                {
                    cout << " WARNING !!! ";
                    cout << " Reference allele for marker : " << currID << " is : " <<refAllele<<".";
                    cout << " Variant will be ignored ... \n";
                    flag=1;
                    inconsistent++;
                }
            }

            if(flag==0)
                switch (altAllele[0])
                    {
                    case 'A': case 'a': break;
                    case 'C': case 'c': break;
                    case 'G': case 'g': break;
                    case 'T': case 't': break;
                    case 'D': case 'd': break;
                    case 'I': case 'i': break;
                    case 'R': case 'r': break;
                    default:
                    {
                        cout << " WARNING !!! ";
                        cout << " Alternate allele for marker : " <<currID << " is : " <<altAllele<<".";
                        cout << " Variant will be ignored ... \n";
                        flag=1;
                    }
                }
        }
        else if(strlen(refAllele.c_str())<strlen(altAllele.c_str()))
        {
//			isInDel = true;
			insertions++;
		}
		else
        {
//			isInDel = true;
			deletions++;
        }

        stringstream strs3,strs4;
        strs3<<(cno);
        strs4<<(bp);

        if(prevID==currID)
        {
            if(refAllele==PrefrefAllele && altAllele==PrevaltAllele)
            {

                cout <<  " WARNING !!!  Duplicate Variant found chr:"
                <<(string)strs3.str()+":"+
                (string)strs4.str()<<" with identical REF = "
                <<refAllele
                 <<" and ALT = "<<altAllele <<".";
                duplicates++;
                if(!Duplicates)
                {
                    cout << "\n              Repeated Instances will be ignored ... ";
                    flag=1;
                }
                cout<<endl;
            }

        }

        prevID=currID;
        PrefrefAllele=refAllele;
        PrevaltAllele=altAllele;


		if(CHR!="")
        {
            if(cno.compare(CHR.c_str())!=0)
                flag=1;
            else
            {
                 if(END>0)
                {
                    if(bp>END || bp<START)
                        flag=1;
                }
                else
                    if(bp<START)
                        flag=1;
            }

        }

		if (flag == 0)
		{
			if(bp<OrigStartPos)
                PrintStartIndex++;


            if(CHR=="")
            {
                PrintEndIndex=numtoBeWrittenRecords;
            }
            else
            {
                if(OrigEndPos==0 || bp<=OrigEndPos)
                    PrintEndIndex=numtoBeWrittenRecords;
            }

			++numtoBeWrittenRecords;
			markerName.push_back(currID);
            variant thisVariant(currID,cno,bp);

            thisVariant.rsid=id;

            VariantList.push_back(thisVariant);
			importIndex.push_back(true);

		}
		else
		{
			importIndex.push_back(false);
		}

	}


	inFile.close();


    if(CHR=="")
        finChromosome=cno;
    else
        finChromosome=CHR.c_str();




	std::cout << "\n Number of Markers read from VCF File                : " << numReadRecords << endl;
	std::cout << " Number of Markers with more than Two Alleles        : " << notBiallelic << endl;
	std::cout << " Number of Markers failing FILTER = PASS             : " << failFilter << endl;
    std::cout << " Number of Markers with inconsistent Ref/Alt Allele  : " << inconsistent << endl;
    std::cout << " Number of Markers with duplicate ID:POS:REF:ALT     : " << duplicates << endl;
    std::cout << " Number of Insertions                                : " << insertions << endl;
	std::cout << " Number of Deletions                                 : " << deletions << endl;


    if(numReadRecords==0)
    {
        cout << "\n No Markers recorded from VCF Input File : "<<filename<<endl;
		cout << " Please check the file properly..\n";
		cout << " Program Aborting ... "<<endl;
		return false;
    }


	if (maxIndiv == 0)
		maxIndiv = numSamplesRead;
	numMarkers = numtoBeWrittenRecords;
	numHaplotypes = (maxIndiv<numSamplesRead) ? (2 * maxIndiv) : (2 * numSamplesRead);
    numSamples=numHaplotypes/2;
    if(finChromosome!="X")
    {
        std::cout << "\n Number of Markers to be Recorded                    : " << numtoBeWrittenRecords << endl;
        std::cout << " Number of Haplotypes to be Recorded                 : " << (numHaplotypes) << endl;
    }
    else
    {
        cout<<"\n Chromosome X Detected !!! \n";
        std::cout << "\n Number of Markers to be Recorded                    : " << numtoBeWrittenRecords << endl;
        PseudoAutosomal=false;
        inFile.open(filename, header);
        inFile.setSiteOnly(false);
        int numWrittenRecords = 0;
//        int readIndex = -1;
        inFile.readRecord(record);
        int tempHapCount=0,MaleCount=0,FemaleCount=0;
        for (int i = 0; i<(numSamples); i++)
        {
            if(record.getNumGTs(i)==0)
            {
                std::cout << "\n Empty Value for Individual : " << individualName[i] << " at Marker : " << VariantList[numWrittenRecords].name << endl;
                std::cout << " Most probably a corrupted VCF file. Please check input VCF file !!! " << endl;
                return false;
            }
            else
            {
                if(record.getNumGTs(i)==1)
                    MaleCount++;
                else
                    FemaleCount++;
                SampleNoHaplotypes.push_back(record.getNumGTs(i));

                tempHapCount+=record.getNumGTs(i);
            }

        }
        inFile.close();
        std::cout << " Number of Samples (Haplotypes) to be Recorded       : " << numSamples << " ("<<tempHapCount <<") "<<endl;
        numHaplotypes=tempHapCount;

        if(MaleCount>0)
        {
        std::cout << " Number of MALE Samples (Haplotypes)                 : " << MaleCount << " ("<<MaleCount <<") "<<endl;
        std::cout << " Number of FEMALE Samples (Haplotypes)               : " << FemaleCount<< " ("<<FemaleCount*2 <<") "<<endl;
        }
        else
        {
            std::cout << "\n All " << numSamples<<" samples have two alleles on Chromosome X !!! "<< endl;
            cout<<" Cannot determine number of male/female samples from Pseudo-Autosomal Region ..."<<endl;
            PseudoAutosomal=true;
        }


    }


   if(numtoBeWrittenRecords<2)
    {
        cout << "\n None/Single marker left after filtering from Input File : "<<filename<<endl;
		cout << " Please check the file or the filtering options properly ...\n";
		cout << " Program Aborting ... "<<endl;
		return false;
    }

    std::cout <<"\n Starting to load data ...";

    std::cout <<"\n\n Number of Threads to be Used = "<<CPU<<endl;



    int    blockSize = 500;
    int    bufferSize = 10000;
       individualName.resize(numSamples);
    // open file again to start loading the data in the variable haplotype.

    cout<<endl;
    VcfFileReader inFileBuffer;
    VcfHeader headerBuffer;
    VcfRecord recordBuffer;
    int readIndex = -1;

    inFileBuffer.open(filename, headerBuffer);
    inFileBuffer.setSiteOnly(false);

    vector<vector<String> > Haplotypes(CPU);



    int currentPiece=0;

    for(currentPiece=0;currentPiece<CPU;currentPiece++)
    {
        Haplotypes[currentPiece].resize(numHaplotypes);
    }
    currentPiece=0;


    vector<int> BufferPosList(CPU);


    int NoMarkersWritten=1;
    int NoMarkersRead=0;

    BufferPosList[0]=1;
    for(int BufferNo=1;BufferNo<CPU;BufferNo++)
    {
        BufferPosList[BufferNo]= BufferPosList[BufferNo-1]+(bufferSize-1);

    }
    ReducedStructureInfoBuffer.clear();
    ReducedStructureInfoBuffer.resize(CPU);

    ReducedStructureInfo.clear();


    while (inFileBuffer.readRecord(recordBuffer) && NoMarkersRead<numMarkers)
	{
		// work only with bi-allelic markers
		readIndex++;


		if (importIndex[readIndex])
		{
            if (Haplotypes[currentPiece][0].Length()==1)
                {
                    printf("  Loading markers %d - %d  out of %d markers to be loaded... [%.1f%%] ",
                           NoMarkersWritten+currentPiece*(bufferSize-1),
                           min(NoMarkersWritten+(currentPiece+1)*(bufferSize-1),numMarkers),numMarkers,
                           100*(double)(NoMarkersWritten+currentPiece*(bufferSize-1))/numMarkers);
                    cout<<endl;
                }

                NoMarkersRead++;

                string refAllele = recordBuffer.getRefStr();
                string altAllele = recordBuffer.getAltStr();
                char Rallele,Aallele;

                if (strlen(refAllele.c_str()) == 1 && strlen(altAllele.c_str()) == 1)
                    {
                        switch (refAllele[0])
                        {
                            case 'A': case 'a': Rallele = 1; break;
                            case 'C': case 'c': Rallele = 2; break;
                            case 'G': case 'g': Rallele = 3; break;
                            case 'T': case 't': Rallele = 4; break;
                            case 'D': case 'd': Rallele = 5; break;
                            case 'I': case 'i': Rallele = 6; break;
                            case 'R': case 'r': Rallele = 7; break;
                            default:
                            {
                                       cout << "\n\n Data Inconsistency !!! \n";
                                       cout << " Error with reference allele for marker : " << recordBuffer.getIDStr() << " in VCF File : " << filename;
                                       cout << "\n VCF reference alleles for SNPs can only be A(a), C(c), G(g), or T(t).\n";
                                       cout << " " << recordBuffer.getIDStr() << " has " << refAllele << endl;
                                       cout << "\n Program Aborting ... \n\n";
                                       abort();
                            }
                        }

                        switch (altAllele[0])
                        {
                            case 'A': case 'a': Aallele = 1; break;
                            case 'C': case 'c': Aallele = 2; break;
                            case 'G': case 'g': Aallele = 3; break;
                            case 'T': case 't': Aallele = 4; break;
                            case 'D': case 'd': Aallele = 5; break;
                            case 'I': case 'i': Aallele = 6; break;
                            case 'R': case 'r': Aallele = 7; break;
                            default:
                            {
                                       cout << "\n\n Data Inconsistency !!! \n";
                                       cout << " Error with alternate allele for marker : " << recordBuffer.getIDStr() << " in VCF File : " << filename;
                                       cout << "\n VCF alternate alleles for SNPs can only be A(a), C(c), G(g), or T(t).\n";
                                       cout << " " << recordBuffer.getIDStr() << " has " << altAllele << endl;
                                       cout << "\n Program Aborting ... \n\n";
                                       abort();
                            }
                        }


                    }

                else
                    {
                        Rallele = 7;
                        if(strlen(refAllele.c_str())<strlen(altAllele.c_str()))
                            Aallele=6;
                        else
                            Aallele=5;
                    }

                refAlleleList.push_back(Rallele);
                VariantList[NoMarkersRead-1].refAllele=Rallele;
                VariantList[NoMarkersRead-1].altAllele=Aallele;
                VariantList[NoMarkersRead-1].refAlleleString=refAllele;
                VariantList[NoMarkersRead-1].altAlleleString=altAllele;
                int haplotype_index = 0;
                for (int i = 0; i<(numSamples); i++)
                {
                    int NoGt=recordBuffer.getNumGTs(i);

                    if(NoGt==0)
                    {
                        std::cout << "\n Empty Value for Individual : " << individualName[i] << " at Marker : " << VariantList[NoMarkersRead-1].name << endl;
                        std::cout << " Most probably a corrupted VCF file. Please check input VCF file !!! " << endl;
                        abort();
                    }
                    else if(NoGt==1 && finChromosome!="X")
                    {
                        std::cout << "\n Single Autosomal Haplotype for Individual : " << individualName[i] << " at Marker : " << VariantList[NoMarkersRead-1].name << endl;
                        std::cout << " Most probably a corrupted VCF file. Please check input VCF file !!! " << endl;
                        abort();
                    }
                    for (int j = 0; j<NoGt; j++)
                    {

                        int alleleIndex = recordBuffer.getGT(i, j);
                        if (alleleIndex<0)
                        {
                            std::cout << "\n Missing Value for Individual : " << individualName[i] << " at Marker : " << VariantList[NoMarkersRead-1].name << endl;
                            abort();
                        }

                        else
                        {
                            if(haplotype_index>=numHaplotypes && finChromosome=="X")
                            {
                                cout << "\n Error in Reference VCF File for Chromosome X !!! "<<endl;
                                cout << "\n Marker : " << VariantList[0].name << " has "<<numHaplotypes <<" haplotypes while ";
                                cout << "Marker : " << VariantList[NoMarkersRead-1].name << " has "<< haplotype_index+1<<" haplotypes."<<endl;
                                cout << " VCF file seems to have both Pseudo-Autosomal region (PAR) and non-PAR of chromosome X. \n";
                                cout << " Please use only either of the two regions ... \n";
                                cout << " See web-page for Minimac3 (Chromosome X Imputation) for more details ...\n";
                                cout << " Program Aborting ... "<<endl;
                                abort();
                            }
                            Haplotypes[currentPiece][haplotype_index]+= (char)('0'+alleleIndex);
                            haplotype_index++;
                        }

                    }


                }

                if(haplotype_index<(numHaplotypes-1) && finChromosome=="X")
                {
                    cout << "\n Error in Reference VCF File for Chromosome X !!! "<<endl;
                    cout << "\n Marker : " << VariantList[0].name << " has "<<numHaplotypes <<" haplotypes while ";
                    cout << "Marker : " << VariantList[NoMarkersRead-1].name << " has "<< haplotype_index+1<<" haplotypes."<<endl;
                    cout << " VCF file seems to have both Pseudo-Autosomal region (PAR) and non-PAR of chromosome X. \n";
                    cout << " Please use only either of the two regions ... \n";
                    cout << " See web-page for Minimac3 (Chromosome X Imputation) for more details ...\n";
                    cout << " Program Aborting ... "<<endl;
                    abort();
                }


                int NewPiece=CPU;
                if (Haplotypes[currentPiece][0].Length() == bufferSize)
                {

                    NewPiece=(currentPiece+1)%CPU;
                    vector<String> tempHaplotypes(numHaplotypes);
                    if(NewPiece!=0)
                    {
                        for (int i = 0; i < numHaplotypes; i++)
                        {
                            tempHaplotypes[i]=Haplotypes[currentPiece][i][Haplotypes[currentPiece][i].Length()-1];
                            Haplotypes[NewPiece][i]=tempHaplotypes[i];
                        }
                        currentPiece=NewPiece;
                    }
                }

                if(NewPiece==0 || NoMarkersRead==numMarkers)
                {

                    #pragma omp parallel for
                    for(int ThisPiece=0;ThisPiece<=currentPiece;ThisPiece++)
                    {

                        int LastflushPos=BufferPosList[ThisPiece]-1;
                        printf("     Processing Reference Chunk %d for M3VCF ... \n",ThisPiece+1);


                        vector<int> index(numHaplotypes),oldIndex;
                        vector<int> previousDifference(numHaplotypes);
                        vector<int> previousPredecessor(numHaplotypes);
                        vector<int> firstDifference(numHaplotypes-1,0);
                        vector<int> cost(bufferSize+1,0);
                        vector<int> bestSlice(bufferSize+1,0);
                        vector<int> bestComplexity(bufferSize+1,0);
                        vector<vector<int> > bestIndex(bufferSize+1);

                        ReducedStructureInfoBuffer[ThisPiece].clear();

                        findUnique RefUnique;

                        RefUnique.updateCoeffs(transFactor,cisFactor);
                        double blockedCost = 0.0;

                        for(int i=0;i<numHaplotypes;i++)
                            index[i]=i;


                        for(int length=1;length<=Haplotypes[ThisPiece][0].Length();length++)
                        {


                            vector<int> offsets(3,0);
                            for (int i = 0; i < numHaplotypes; i++)
                                offsets[Haplotypes[ThisPiece][i][length - 1] - '0' + 1]++;
                            offsets[2]+=offsets[1];


                            oldIndex = index;
                            for (int i = 0; i < numHaplotypes; i++)
                                {
                                    index[offsets[Haplotypes[ThisPiece][oldIndex[i]][length - 1] - '0']++] = oldIndex[i];
                                }

                            RefUnique.UpdateDeltaMatrix(Haplotypes[ThisPiece], index, firstDifference, length, blockSize,
                                   oldIndex, previousPredecessor, previousDifference);

                            RefUnique.AnalyzeBlocks(index, firstDifference, length, blockSize,
                               cost, bestSlice, bestComplexity, bestIndex);

                        }


                        if(Haplotypes[ThisPiece][0].Length()>1)
                            blockedCost += RefUnique.FlushBlocks(ReducedStructureInfoBuffer[ThisPiece],
                                                                VariantList,LastflushPos,
                                                                Haplotypes[ThisPiece], cost,
                                                                bestComplexity, bestSlice, bestIndex);

                    }


                    cout<<endl;
                    NoMarkersWritten+=(CPU*(bufferSize-1));


                    BufferPosList[0]=NoMarkersWritten;
                    for(int BufferNo=1;BufferNo<CPU;BufferNo++)
                    {
                        BufferPosList[BufferNo]= BufferPosList[BufferNo-1]+(bufferSize-1);

                    }

                    vector<String> tempHaplotypes(numHaplotypes);
                    for (int i = 0; i < numHaplotypes; i++)
                    {
                        tempHaplotypes[i]=Haplotypes[currentPiece][i][Haplotypes[currentPiece][i].Length()-1];
                        Haplotypes[0][i]=tempHaplotypes[i];
                    }

                    for(int ThisPiece=0;ThisPiece<CPU;ThisPiece++)
                    {
                        for(int jj=0;jj<(int)ReducedStructureInfoBuffer[ThisPiece].size();jj++)
                            {
                                ReducedStructureInfo.push_back(ReducedStructureInfoBuffer[ThisPiece][jj]);
                            }
                        ReducedStructureInfoBuffer[ThisPiece].clear();
                    }
                    currentPiece=0;
                }
        }
    }





    std::cout << "\n Number of Markers Recorded                          : " << markerName.size() << endl;
    std::cout << " Number of Haplotypes Recorded                       : " << numHaplotypes << endl;




    optEndPoints.clear();
    int i;
    for(i=0;i<(int)ReducedStructureInfo.size();i++)
        {
            optEndPoints.push_back(ReducedStructureInfo[i].startIndex);
        }
    optEndPoints.push_back(ReducedStructureInfo[i-1].endIndex);


    numMarkers = markerName.size();
    std::cout << "\n Haplotype Set successfully loaded from VCF File     : " << filename << endl;
	inFile.close();

    std::cout << "\n Writing draft [.m3vcf] file                         : " <<outFile+".draft"+".m3vcf" + (gzip ? ".gz" : "")<< endl<<endl;

    writem3vcfFile(outFile+".draft",gzip);

	return true;

}


bool HaplotypeSet::BasicCheckForTargetHaplotypes(String filename)
{

	VcfFileReader inFile;
	VcfHeader header;
	VcfRecord record;
	std::cout << "\n Performing basic file check on target/GWAS haplotype file : "<<filename << endl;

    std::cout << "\n Checking File ..." << endl;

	inFile.setSiteOnly(true);
    IFILE fileStream = ifopen(filename, "r");

    string line;
    if(!fileStream)
    {
		cout << "\n Program could NOT open file : " << filename << endl;
		return false;
	}
	else
    {
        std::cout << " File Exists ..." << endl;

        std::cout << "\n Checking File Type ..." << endl;

        fileStream->readLine(line);
        if(line.length()<1)
        {
            cout << "\n Target File provided by \"--haps\" must be a VCF file !!! \n";
            cout << " Please check the following file : "<<filename<<endl;
            return false;
        }

        string tempString;
        tempString=(line.substr(0,17));

        char temp[tempString.length() + 1];
        std::strcpy(temp,tempString.c_str());
        for (char *iter = temp; *iter != '\0'; ++iter)
        {
           *iter = std::tolower(*iter);
        }


        if(((string)temp).compare("##fileformat=vcfv")!=0)
        {
            cout << "\n Target File provided by \"--haps\" must be a VCF file !!! \n";
            cout << " Please check the following file : "<<filename<<endl<<endl;
            return false;
        }

    }

    ifclose(fileStream);

    std::cout << " VCF File Type Detected ..." << endl;

    std::cout <<"\n NOTE: Samples will be assumed to be phased irrespective "<<endl;
    cout<<       "       of GT delimiter (| or /) in VCF file !!! " << endl;


    std::cout << "\n Checking variant information ..." << endl;

    int numActualRecords=0,numReadRecords=0;
    int failFilter=0,duplicates=0,notBiallelic=0,inconsistent=0;
    string prevID,currID, refAllele,altAllele,PrefrefAllele,PrevaltAllele,cno,fixCno,id;
    inFile.open(filename, header);



    while (inFile.readRecord(record))
    {

        int flag=0;
        cno=record.getChromStr();
        int bp=record.get1BasedPosition();
        id=record.getIDStr();
        refAllele = record.getRefStr();
        altAllele = record.getAltStr();

        if(numActualRecords==0)
        {
            if(!CheckValidChrom(cno))
            {
                cout << "\n Error !!! Target VCF File contains chromosome : "<<cno<<endl;
                cout << " VCF File can only contain chromosomes 1-22 and X !!! "<<endl;
                cout << " Program Aborting ... "<<endl;
                return false;
            }

        }

        if (record.getNumAlts()>1)
		{
			notBiallelic++;
			flag = 1;
		}
		if (record.getFilter().getString(0).compare("PASS") != 0)
		{
			failFilter++;
			flag = 1;
		}

        stringstream strs3,strs4;
        strs3<<(cno);
        strs4<<(bp);

        currID=(string)strs3.str()+":"+(string)strs4.str();


        if(!CheckValidChrom(cno))
        {
            cout << "\n Error !!! Target VCF File contains chromosome : "<<cno<<endl;
            cout << " VCF File can only contain chromosomes 1-22 and X !!! "<<endl;
            cout << " Program Aborting ... "<<endl;
            return false;
        }
        if(fixCno!=cno && numReadRecords>0)
        {
            cout << "\n Error !!! Target VCF File contains multiple chromosomes : "<<cno<<", "<<fixCno<<", ... "<<endl;
            cout << " Please use VCF file with single chromosome !!! "<<endl;
            cout << " Program Aborting ... "<<endl;
            return false;
        }

        fixCno=cno;
        if (strlen(refAllele.c_str()) == 1 && strlen(altAllele.c_str()) == 1)
        {
            switch (refAllele[0])
            {
                case 'A': case 'a': ; break;
                case 'C': case 'c': ; break;
                case 'G': case 'g': ; break;
                case 'T': case 't': ; break;
                case 'D': case 'd': ; break;
                case 'I': case 'i': ; break;
                case 'R': case 'r': ; break;
                default:
                {
                        flag=1;
                        inconsistent++;
                }
            }
            if(flag==0)
                switch (altAllele[0])
                {
                    case '0':  ; break;
                    case 'A': case 'a': ; break;
                    case 'C': case 'c': ; break;
                    case 'G': case 'g': ; break;
                    case 'T': case 't': ; break;
                    case 'D': case 'd': ; break;
                    case 'I': case 'i': ; break;
                    case 'R': case 'r': ; break;
                    default:
                    {
                        flag=1;
                        inconsistent++;
                    }
                }
        }


        if(prevID==currID)
        {
            if(refAllele==PrefrefAllele && altAllele==PrevaltAllele)
            {
                duplicates++;
                cout << "\n Error !!! Duplicate Variant found chr:"<<cno<<":"<<bp<<" with identical REF = "<<refAllele <<" and ALT = "<<altAllele <<"\n";
                cout << " Program Aborting ... "<<endl;
                return false;
            }

        }

        prevID=currID;
        PrefrefAllele=refAllele;
        PrevaltAllele=altAllele;

        if(flag==0)
        {
            ++numReadRecords;
        }
        numActualRecords++;

    }



	std::cout << " "<<numActualRecords<<" variants in file with "<<numReadRecords<<" variants passing filters ..."<<endl<<endl;
    std::cout << " Checking sample information ..." << endl;

	inFile.close();

    inFile.open(filename, header);
    inFile.setSiteOnly(false);

	inFile.readRecord(record);

    std::cout << " "<<header.getNumSamples()<<" samples found in file ..."<<endl;

    if(header.getNumSamples()==0)
    {
        cout << "\n No haplotypes recorded from VCF Input File : "<<filename<<endl;
		cout << " Please check the file properly..\n";
		cout << " Program Aborting ... "<<endl;
		return false;
    }

    inFile.close();

    finChromosome=cno;

    if(finChromosome=="X")
    {
        cout<<"\n Checking Chromosome X information !!! \n";
        cout << " NOTE: All Samples in target VCF file must be either only MALE or only FEMALE ...  " << endl;

        inFile.open(filename, header);
        inFile.setSiteOnly(false);
        inFile.readRecord(record);
        int tempHapFlag=0;
        int tempHapCount=0;
        for (int i = 0; i<(header.getNumSamples()); i++)
        {
            if(record.getNumGTs(i)==0)
            {
                std::cout << "\n Empty Value for Individual at first variant !!! " << endl;
                std::cout << " Most probably a corrupted VCF file. Please check input VCF file !!! " << endl;
                return false;
            }
            else
            {
                if(tempHapFlag!=0 && tempHapFlag!=record.getNumGTs(i))
                {
                    cout << "\n ERROR: Both haploid and diploid samples found in target VCF file "<<filename<<endl;
                    cout << " For chromosome X imputation, Males and Females should be imputed separately.\n";
                    cout << " Program Aborting ... "<<endl;
                    return false;

                }
                tempHapFlag=record.getNumGTs(i);
                tempHapCount+=record.getNumGTs(i);
            }

        }
        inFile.close();
        std::cout << " "<<tempHapCount<<" haplotypes found in file "<<endl;
    }


    std::cout << "\n Initial basic file check on target/GWAS haplotype file successful !!!" << endl<<endl;
	return true;

}


bool HaplotypeSet::LoadTargetHaplotypes(String filename, String targetSnpfile, vector<string> &refSnpList,HaplotypeSet &rHap,bool typedOnly,bool passOnly)
{
	string FileType=DetectTargetFileType(filename);
    TypedOnly=typedOnly;
    Filter=passOnly;


    if(FileType.compare("vcf")==0)
    {
        cout<<"\n Format = VCF (Variant Call Format) "<<endl;
        return LoadVcfTargetHaplotypes(filename, targetSnpfile, refSnpList,rHap);
    }
    else if(FileType.compare("Invalid")==0)
    {

        cout << "\n Target File provided by \"--haps\" must be a VCF or MaCH file !!! \n";
        cout << " Please check the following file : "<<filename<<endl<<endl;
        return false;
    }
    else if(FileType.compare("NA")==0)
    {
        cout<<"\n Following File File Not Available : "<<filename<<endl<<endl;
        return false;
    }
    return false;
}

bool HaplotypeSet::LoadVcfTargetHaplotypes(String filename, String snpNames, vector<string> &refSnpList, HaplotypeSet &rHap)
{


	VcfFileReader inFile;
	VcfHeader header;
	VcfRecord record;
    string RefCno=rHap.VariantList[0].chr;


	inFile.setSiteOnly(true);
	if (!inFile.open(filename, header))
	{
		cout << "\n Program could NOT open file : " << filename << endl;
		return false;
	}
    vector<int> importIndexList;

    int numReadRecords = 0,numActualRecords=0;
    int failFilter=0,notBiallelic=0,inconsistent=0,otherChrom=0;
    std::cout << "\n Loading Target Haplotype SNP List from VCF File     : " << filename << endl<<endl;
    string cno;
    while (inFile.readRecord(record))
    {

        int flag=0;
        cno=record.getChromStr();
        int bp=record.get1BasedPosition();
        string id=record.getIDStr();
        string currID;
        string refAllele = record.getRefStr();
        string altAllele = record.getAltStr();
        char Rallele=0,Aallele=0;


        if(cno!=RefCno)
        {
            otherChrom++;
            flag=1;
        }

        if (record.getNumAlts()>1)
		{
			notBiallelic++;
			flag = 1;
		}
		if (record.getFilter().getString(0).compare("PASS") != 0)
		{
			failFilter++;
			if(Filter)
                flag = 1;
		}
        stringstream strs3,strs4;
        strs3<<(cno);
        strs4<<(bp);


        currID=(string)strs3.str()+":"+(string)strs4.str();


        if(id==".")
            id=currID;

        if (strlen(refAllele.c_str()) == 1 && strlen(altAllele.c_str()) == 1)
        {
            switch (refAllele[0])
            {
                case 'A': case 'a': Rallele = 1; break;
                case 'C': case 'c': Rallele = 2; break;
                case 'G': case 'g': Rallele = 3; break;
                case 'T': case 't': Rallele = 4; break;
                case 'D': case 'd': Rallele = 5; break;
                case 'I': case 'i': Rallele = 6; break;
                case 'R': case 'r': Rallele = 7; break;
                default:
                {
                        cout << " WARNING !!! Reference allele for SNP for "<<currID<<"is "<<refAllele<<". Will be ignored ..." <<endl;
                        flag=1;
                        inconsistent++;
                }
            }
            if(flag==0)
                switch (altAllele[0])
                {
                    case '0':  Aallele = 0; break;
                    case 'A': case 'a': Aallele = 1; break;
                    case 'C': case 'c': Aallele = 2; break;
                    case 'G': case 'g': Aallele = 3; break;
                    case 'T': case 't': Aallele = 4; break;
                    case 'D': case 'd': Aallele = 5; break;
                    case 'I': case 'i': Aallele = 6; break;
                    case 'R': case 'r': Aallele = 7; break;
                    default:
                    {
                        cout << " WARNING !!! Alternate allele for SNP for "<<currID<<"is "<<altAllele<<". Will be ignored ..." <<endl;
                        flag=1;
                        inconsistent++;
                    }
                }
        }
        else
        {
            Rallele = 7;
            if(strlen(refAllele.c_str())<strlen(altAllele.c_str()))
                Aallele=6;
            else
                Aallele=5;
        }


        variant thisVariant(currID,cno,bp);
        if(flag==0)
        {
            VariantList.push_back(thisVariant);

            VariantList[numReadRecords].refAllele=Rallele;
            VariantList[numReadRecords].altAllele=Aallele;
            VariantList[numReadRecords].refAlleleString=refAllele;
            VariantList[numReadRecords].altAlleleString=altAllele;
            VariantList[numReadRecords].rsid=id;
            markerName.push_back(currID);
            ++numReadRecords;
        }
        numActualRecords++;
        importIndexList.push_back(flag);

    }

	std::cout << "\n Number of Markers read from VCF File                : " << numActualRecords << endl;
	std::cout << " Number of Markers with more than One Allele         : " << notBiallelic << endl;
	std::cout << " Number of Markers failing FILTER = PASS             : " << failFilter << endl;
    std::cout << " Number of Markers with inconsistent Ref/Alt Allele  : " << inconsistent << endl;
    std::cout << " Number of Markers on other chromosomes (Non-Ref)    : " << otherChrom << endl;

    if(numActualRecords==0)
    {
        cout << "\n No Markers recorded from VCF Input File : "<<filename<<endl;
		cout << " Please check the file properly..\n";
		cout << " Program Aborting ... "<<endl;
		return false;
    }


    int refMarkerCount=(int)rHap.VariantList.size();

    vector<int> TargetMissing;
//    RefTypedIndex.resize(numReadRecords,0)


    finChromosome=cno;
	vector<int> knownPosition;
	knownPosition.clear();
	int counter = 0;
	GWASOnlycounter = 0;
	missing.resize(refMarkerCount, true);
	ScaffoldIndex.resize(refMarkerCount, -1);
    UnScaffoldIndex.resize(VariantList.size(), -1);
    AllMaleTarget=false;

    int flag;
	int markerIndex=0;
	vector<string> newMarkerName;
	newMarkerName.clear();
    vector<bool> RefAlleleSwap;

    RefAlleleSwap.clear();
    vector<bool> GWASOnlySkipIndex;
    int GWASCounterCount = 0 ;
	for (int j = 0; j<(int)VariantList.size(); j++)
	{

		int prevCounter = counter;
		flag=0;
		while(counter<refMarkerCount && flag==0 && rHap.VariantList[counter].bp<=VariantList[j].bp)
        {

            if(rHap.VariantList[counter].chr==VariantList[j].chr
             && rHap.VariantList[counter].bp==VariantList[j].bp)
            {
                prevCounter = counter;

                if(rHap.VariantList[counter].refAlleleString==VariantList[j].refAlleleString
                        && rHap.VariantList[counter].altAlleleString==VariantList[j].altAlleleString)
                    flag=1;
                else if(rHap.VariantList[counter].refAlleleString==VariantList[j].altAlleleString
                        && rHap.VariantList[counter].altAlleleString==VariantList[j].refAlleleString)
                    flag=1;
                else if (VariantList[j].refAlleleString==VariantList[j].altAlleleString
                        && rHap.VariantList[counter].refAlleleString==VariantList[j].refAlleleString)
                    flag=1;
                else if (VariantList[j].refAlleleString==VariantList[j].altAlleleString
                        && rHap.VariantList[counter].altAlleleString==VariantList[j].refAlleleString)
                    flag=1;
                else
                    counter++;
            }
            else
                counter++;



        }
        if(flag==1)
        {
            knownPosition.push_back(counter);
            newMarkerName.push_back(markerName[j]);

            if(rHap.VariantList[counter].refAlleleString==VariantList[j].refAlleleString)
                RefAlleleSwap.push_back(false);
            else
                {

                    RefAlleleSwap.push_back(true);
                    string tempa=VariantList[j].refAlleleString;
                    VariantList[j].refAlleleString=VariantList[j].altAlleleString;
                    VariantList[j].altAlleleString=tempa;
                }

            missing[counter] = false;
            UnScaffoldIndex[markerIndex]=counter;
            ScaffoldIndex[counter]=markerIndex++;

            counter++;
		}
		else
        {

            if(TypedOnly)
            {
                if(rHap.OrigEndPos>0)
                {
                    if(VariantList[j].bp>=rHap.OrigStartPos
                       && VariantList[j].bp<=rHap.OrigEndPos)

                    {
                        TargetMissing.push_back(counter-1);
                        GWASOnlycounter++;
                        GWASOnlySkipIndex.push_back(true);
                    }
                    else
                    {
                        GWASOnlySkipIndex.push_back(false);
                    }

                }
                else
                {

                    TargetMissing.push_back(counter-1);
                    GWASOnlycounter++;
                    GWASOnlySkipIndex.push_back(true);
                }


            }

            knownPosition.push_back(-1);
            counter = prevCounter;
            GWASCounterCount++;
        }

	}

    counter=0;
    rHap.RefTypedIndex.clear();
    int ThisIndex=0;

	while(counter<refMarkerCount && ThisIndex<(int)TargetMissing.size())
    {
        if(counter<=TargetMissing[ThisIndex])
        {
            rHap.RefTypedIndex.push_back(-1);
            counter++;
        }
        else
        {
            rHap.RefTypedIndex.push_back(ThisIndex);
            ThisIndex++;
        }
    }
    while(counter<refMarkerCount)
    {
        rHap.RefTypedIndex.push_back(-1);
        counter++;
    }
    while(ThisIndex<(int)TargetMissing.size())
    {
        rHap.RefTypedIndex.push_back(ThisIndex);
        ThisIndex++;
    }

    rHap.RefTypedTotalCount=GWASOnlycounter+rHap.numMarkers;


    if(rHap.RefTypedTotalCount!=(int)rHap.RefTypedIndex.size())
    {
        cout<<endl<<endl<<" Error in Code Construction [ERROR: 007] !!! "<<endl;
        cout<<" Please Contact author with ERROR number urgently : sayantan@umich.edu "<<endl;
        cout<<" Program Exiting ..."<<endl<<endl;
        abort();
    }



	numHaplotypes = 0;
    numMarkers = 0;

    cout<<endl;
	std::cout << " Number of Markers overlapping with Reference List   : " << (int)newMarkerName.size()  << endl;
    std::cout << " Number of Markers Only in Target Panel              : " << GWASCounterCount<< endl;
    std::cout << " Number of Markers to be Recorded                    : " << (int)newMarkerName.size() + GWASOnlycounter << endl<<endl;




	if (newMarkerName.size() == 0)
	{

		cout << "\n No overlap between Target and Reference markers !!!\n";
		cout << " Please check for consistent marker identifier in reference and target input files..\n";
		cout << " Program Aborting ... \n";
		return false;

	}

    TypedOnlyVariantList.clear();
//	vector<string> AllMarkerName=markerName;
	markerName = newMarkerName;
    numHaplotypes = 2 * header.getNumSamples();
    numSamples=header.getNumSamples();
    if(numHaplotypes==0)
    {
        cout << "\n No haplotypes recorded from VCF Input File : "<<filename<<endl;
		cout << " Please check the file properly..\n";
		cout << " Program Aborting ... "<<endl;
		return false;
    }


    for (int i = 0; i < numSamples; i++)
	{
		string tempName(header.getSampleName(i));
		individualName.push_back(tempName);
	}

    if(finChromosome=="X")
    {
        inFile.open(filename, header);
        inFile.setSiteOnly(false);
        inFile.readRecord(record);
        int tempHapCount=0;
        for (int i = 0; i<(numSamples); i++)
        {
            if(record.getNumGTs(i)==0)
            {
                std::cout << "\n Empty Value for Individual : " << individualName[i] << " at Marker : " << VariantList[0].name << endl;
                std::cout << " Most probably a corrupted VCF file. Please check input VCF file !!! " << endl;
                return false;
            }
            else
            {
                SampleNoHaplotypes.push_back(record.getNumGTs(i));
                tempHapCount+=record.getNumGTs(i);
            }
            if(i==0)
            {
                if(record.getNumGTs(i)==1)
                    AllMaleTarget=true;
                else
                    AllMaleTarget=false;
            }
        }
        inFile.close();
        numHaplotypes=tempHapCount;
    }



	inFile.close();

	std::cout << " Loading Target Haplotype Set from VCF File          : " << filename << endl<<endl;


	inFile.open(filename, header);
    inFile.setSiteOnly(false);

	inFile.readRecord(record);


	haplotypesUnscaffolded.resize(numHaplotypes);
	MissingSampleUnscaffolded.resize(numHaplotypes);


	if(TypedOnly)
    {
        GWASOnlyhaplotypesUnscaffolded.resize(numHaplotypes);
        GWASOnlyMissingSampleUnscaffolded.resize(numHaplotypes);
    }



	for (int i = 0; i<numHaplotypes; i++)
	{
		haplotypesUnscaffolded[i].resize(markerName.size(), false);
        MissingSampleUnscaffolded[i].resize(markerName.size(), false);
        if(TypedOnly)
        {
            GWASOnlyhaplotypesUnscaffolded[i].resize(GWASOnlycounter, false);
            GWASOnlyMissingSampleUnscaffolded[i].resize(GWASOnlycounter, false);
        }
	}


	// initialize this to 0 again to read each marker from each line of vcf file.
    int importIndex=-1;
    //ScaffoldIndex=knownPosition;
	int readIndex = -1;
	int numtoBeWrittenRecords = 0;
	int GWASnumtoBeWrittenRecords = 0;
	int GWASOnlySkipIndexpoint = 0;

	do
	{
		// work only with bi-allelic markers
		readIndex++;
        if (importIndexList[readIndex] ==0)
		{
		    importIndex++;

            if (knownPosition[importIndex] != -1)
            {
//                RefTypedIndex.push_back(numtoBeWrittenRecords);
                if (numtoBeWrittenRecords % 10000 == 0)
                       printf("  Loading markers %d out of %d markers to be loaded... [%.1f%%] \n",numtoBeWrittenRecords+1,(int)markerName.size(),100*(double)(numtoBeWrittenRecords + 1)/(int)markerName.size());

//                tempMarkerNames.push_back(record.getIDStr());
                string refAllele = record.getRefStr();
                string altAllele = record.getAltStr();
                char Rallele;
                if (strlen(refAllele.c_str()) == 1 && strlen(altAllele.c_str()) == 1)
                {
                    switch (refAllele[0])
                    {
                        case 'A': case 'a': Rallele = 1; break;
                        case 'C': case 'c': Rallele = 2; break;
                        case 'G': case 'g': Rallele = 3; break;
                        case 'T': case 't': Rallele = 4; break;
                        case 'D': case 'd': Rallele = 5; break;
                        case 'I': case 'i': Rallele = 6; break;
                        case 'R': case 'r': Rallele = 7; break;
                        default:
                        {
                                   cout << "\n Data Inconsistency !!! \n";
                                   cout << " Error with reference allele for marker : " << record.getIDStr() << " in VCF File : " << filename;
                                   cout << "\n VCF reference alleles for SNPs can only be A(a), C(c), G(g), or T(t).\n";
                                   cout << " " << record.getIDStr() << " has " << refAllele << endl;
                                   cout << " Program Aborting ... \n\n";
                                   return false;
                        }
                    }
                }
                else
                    Rallele = 7;

                refAlleleList.push_back(Rallele);

                int haplotype_index = 0;
                for (int i = 0; i<(numSamples); i++)
                {

                    if(record.getNumGTs(i)==0)
                    {
                        std::cout << "\n Empty Value for Individual : " << individualName[i] << " at Marker : " <<
                         VariantList[importIndex].name << endl;
                        std::cout << " Most probably a corrupted VCF file. Please check input VCF file !!! " << endl;
                        return false;
                    }
                    if(record.getNumGTs(i)==1 && finChromosome!="X")
                    {
                        std::cout << "\n Single Autosomal Haplotype for Individual : " << individualName[i] << " at Marker : "
                         << VariantList[importIndex].name << endl;
                        std::cout << " Most probably a corrupted VCF file. Please check input VCF file !!! " << endl;
                        return false;
                    }

                    if(rHap.PseudoAutosomal)
                    {
                        if(record.getNumGTs(i)!=2)
                        {
                            cout << "\n ERROR: Reference Panel on chromosome X seems to be on Pseudo-Autosomal Region "<<endl;
                            cout <<   "        since all samples have two alleles on Chromosome X. "<<endl;
                            cout <<   "        However, "<<VariantList[importIndex].name<<" in target panel has SINGLE haplotypes !!!\n";
                            cout << " Program Aborting ... "<<endl;
                            return false;
                        }
                    }
                    else
                    {
                        if(AllMaleTarget && record.getNumGTs(i)!=1)
                        {
                            cout << "\n ERROR: Reference Panel on chromosome X seems to be on NON-Pseudo-Autosomal Region "<<endl;
                            cout <<   "        since MALE samples had single alleles on Chromosome X. "<<endl;
                            cout <<   "        However, "<<VariantList[importIndex].name<<" in target panel (all MALE samples) \n";
                            cout <<   "        has DOUBLE haplotypes !!!\n";
                            cout << " Program Aborting ... "<<endl;
                            return false;
                        }
                        if(!AllMaleTarget && record.getNumGTs(i)!=2)
                        {
                            cout << "\n ERROR: Reference Panel on chromosome X seems to be on NON-Pseudo-Autosomal Region "<<endl;
                            cout <<   "        since MALE samples had single alleles on Chromosome X. "<<endl;
                            cout <<   "        However, "<<VariantList[importIndex].name<<" in target panel (all FEMALE samples) \n";
                            cout <<   "        has SINGLE haplotypes !!!\n";
                            cout << " Program Aborting ... "<<endl;
                            return false;
                        }



                    }

                    for (int j = 0; j<record.getNumGTs(i); j++)
                    {

                        int alleleIndex = record.getGT(i, j);
                        if (alleleIndex<0)
                        {
                            MissingSampleUnscaffolded[haplotype_index][numtoBeWrittenRecords] = true;
                        }
                        else
                        {

                            if(!RefAlleleSwap[numtoBeWrittenRecords])
                            {
                                if(alleleIndex==1)
                                {
                                    haplotypesUnscaffolded[haplotype_index][numtoBeWrittenRecords] = true;
                                }
                            }
                            else
                            {
                                if(alleleIndex==0)
                                {
                                    haplotypesUnscaffolded[haplotype_index][numtoBeWrittenRecords] = true;
                                }
                            }


                        }

                        haplotype_index++;



                    }
                }
            ++numtoBeWrittenRecords;
            }
            else if(TypedOnly)
            {

                if(GWASOnlySkipIndex[GWASOnlySkipIndexpoint])
                {
                    TypedOnlyVariantList.push_back(VariantList[importIndex]);
                    string refAllele = record.getRefStr();
                    string altAllele = record.getAltStr();
                    char Rallele;
                    if (strlen(refAllele.c_str()) == 1 && strlen(altAllele.c_str()) == 1)
                    {
                        switch (refAllele[0])
                        {
                            case 'A': case 'a': Rallele = 1; break;
                            case 'C': case 'c': Rallele = 2; break;
                            case 'G': case 'g': Rallele = 3; break;
                            case 'T': case 't': Rallele = 4; break;
                            case 'D': case 'd': Rallele = 5; break;
                            case 'I': case 'i': Rallele = 6; break;
                            case 'R': case 'r': Rallele = 7; break;
                            default:
                            {
                                       cout << "\n Data Inconsistency !!! \n";
                                       cout << " Error with reference allele for marker : " << record.getIDStr() << " in VCF File : " << filename;
                                       cout << "\n VCF reference alleles for SNPs can only be A(a), C(c), G(g), or T(t).\n";
                                       cout << " " << record.getIDStr() << " has " << refAllele << endl;
                                       cout << " Program Aborting ... \n\n";
                                       return false;
                            }
                        }
                    }
                    else
                        Rallele = 7;

                    refAlleleList.push_back(Rallele);

                    int haplotype_index = 0;
                    for (int i = 0; i<(numSamples); i++)
                    {

                        if(record.getNumGTs(i)==0)
                        {
                            std::cout << "\n Empty Value for Individual : " << individualName[i] << " at Marker : " <<
                             VariantList[importIndex].name << endl;
                            std::cout << " Most probably a corrupted VCF file. Please check input VCF file !!! " << endl;
                            return false;
                        }
                        if(record.getNumGTs(i)==1 && finChromosome!="X")
                        {
                            std::cout << "\n Single Autosomal Haplotype for Individual : " << individualName[i] << " at Marker : "
                             << VariantList[importIndex].name << endl;
                            std::cout << " Most probably a corrupted VCF file. Please check input VCF file !!! " << endl;
                            return false;
                        }

                        if(rHap.PseudoAutosomal)
                        {
                            if(record.getNumGTs(i)!=2)
                            {
                                cout << "\n ERROR: Reference Panel on chromosome X seems to be on Pseudo-Autosomal Region "<<endl;
                                cout <<   "        since all samples have two alleles on Chromosome X. "<<endl;
                                cout <<   "        However, "<<VariantList[importIndex].name<<" in target panel has SINGLE haplotypes !!!\n";
                                cout << " Program Aborting ... "<<endl;
                                return false;
                            }
                        }
                        else
                        {
                            if(AllMaleTarget && record.getNumGTs(i)!=1)
                            {
                                cout << "\n ERROR: Reference Panel on chromosome X seems to be on NON-Pseudo-Autosomal Region "<<endl;
                                cout <<   "        since MALE samples had single alleles on Chromosome X. "<<endl;
                                cout <<   "        However, "<<VariantList[importIndex].name<<" in target panel (all MALE samples) \n";
                                cout <<   "        has DOUBLE haplotypes !!!\n";
                                cout << " Program Aborting ... "<<endl;
                                return false;
                            }
                            if(!AllMaleTarget && record.getNumGTs(i)!=2)
                            {
                                cout << "\n ERROR: Reference Panel on chromosome X seems to be on NON-Pseudo-Autosomal Region "<<endl;
                                cout <<   "        since MALE samples had single alleles on Chromosome X. "<<endl;
                                cout <<   "        However, "<<VariantList[importIndex].name<<" in target panel (all FEMALE samples) \n";
                                cout <<   "        has SINGLE haplotypes !!!\n";
                                cout << " Program Aborting ... "<<endl;
                                return false;
                            }



                        }

                        for (int j = 0; j<record.getNumGTs(i); j++)
                        {


                            int alleleIndex = record.getGT(i, j);
                            if (alleleIndex<0)
                            {
                                GWASOnlyMissingSampleUnscaffolded[haplotype_index][GWASnumtoBeWrittenRecords] = true;
                            }
                            else
                            {
                                if(alleleIndex==1)
                                {
                                    GWASOnlyhaplotypesUnscaffolded[haplotype_index][GWASnumtoBeWrittenRecords] = true;
                                }
                            }

                            haplotype_index++;



                        }
                    }
                    ++GWASnumtoBeWrittenRecords;
                }
                ++GWASOnlySkipIndexpoint;
            }



		}
	} while (inFile.readRecord(record));



	std::cout << "\n\n Number of Haplotypes Recorded                       : " << (haplotypesUnscaffolded.size()) << endl;
    std::cout << " Number of Markers Recorded for Imputation           : " << markerName.size() << endl;

    if(TypedOnly)
        std::cout << " Number of Markers Only in GWAS Panel Recorded       : " <<GWASOnlycounter << endl;

    numMarkers = markerName.size();

    std::cout << "\n Haplotype Set successfully loaded from VCF File     : " << filename << endl;
	inFile.close();



	return true;

}

void HaplotypeSet::calculateFreq()
{
    alleleFreq.resize(8);
	for (int j = 0; j<8; j++)
		alleleFreq[j].resize(numMarkers, 0.0);

    int i,j,k;
    for(k=0;k<(int)ReducedStructureInfo.size();k++)
    {
        for (i = 0; i<(int)ReducedStructureInfo[k].uniqueCardinality.size(); i++)
        {
            for(j=ReducedStructureInfo[k].startIndex;j<ReducedStructureInfo[k].endIndex;j++)
            {
                alleleFreq[ReducedStructureInfo[k].uniqueHaps[i][j-ReducedStructureInfo[k].startIndex]][j]+=ReducedStructureInfo[k].uniqueCardinality[i];
            }

            if(k==(int)ReducedStructureInfo.size()-1)
                alleleFreq[ReducedStructureInfo[k].uniqueHaps[i][j-ReducedStructureInfo[k].startIndex]][j]+=ReducedStructureInfo[k].uniqueCardinality[i];
        }
    }
	major.resize(numMarkers, 0);
	minor.resize(numMarkers, 0);

	for (int i = 0; i<numMarkers; i++)
	{
		double max_freq = 0.0;
		for (int j = 0; j<8; j++)
		{
			if (max_freq<alleleFreq[j][i])
			{
				max_freq = alleleFreq[j][i];
				major[i] = j;
			}
			//cout<<alleleFreq[j][i]<<"\t";
			alleleFreq[j][i] /= (double)numHaplotypes;
		}
		//cout<<endl;
		if(major[i]==VariantList[i].refAllele)
            {
                minor[i]=VariantList[i].altAllele;
                VariantList[i].MinAlleleString=VariantList[i].altAlleleString;
                VariantList[i].MajAlleleString=VariantList[i].refAlleleString;
            }

        else
            {
                minor[i]=VariantList[i].refAllele;
                VariantList[i].MinAlleleString=VariantList[i].refAlleleString;
                VariantList[i].MajAlleleString=VariantList[i].altAlleleString;
            }


	}

}

void HaplotypeSet::CalculateFreq()
{

    AlleleFreq.resize(numMarkers, 0.0);

    int i,j,k;
    for(k=0;k<(int)ReducedStructureInfo.size();k++)
    {
        for (i = 0; i<(int)ReducedStructureInfo[k].uniqueCardinality.size(); i++)
        {
            for(j=ReducedStructureInfo[k].startIndex;j<ReducedStructureInfo[k].endIndex;j++)
            {
                if(ReducedStructureInfo[k].uniqueHaps[i][j-ReducedStructureInfo[k].startIndex])
                    {
//                         cout<<ReducedStructureInfo[k].uniqueHaps[i][j-ReducedStructureInfo[k].startIndex]<<"\t"<<AlleleFreq[j]<<endl;

                        AlleleFreq[j]+=ReducedStructureInfo[k].uniqueCardinality[i];

//                        cout<<" WELL ="<<ReducedStructureInfo[k].uniqueHaps[i][j-ReducedStructureInfo[k].startIndex]<<"\t"<<AlleleFreq[j]<<endl;

                    }
            }
//abort();

            if(k==(int)ReducedStructureInfo.size()-1)
                if(ReducedStructureInfo[k].uniqueHaps[i][j-ReducedStructureInfo[k].startIndex])
                    AlleleFreq[j]+=ReducedStructureInfo[k].uniqueCardinality[i];
        }
    }
	major.resize(numMarkers, false);
//	minor.resize(numMarkers, 0);

	for (int i = 0; i<numMarkers; i++)
	{

		AlleleFreq[i] /= (double)numHaplotypes;

		if (AlleleFreq[i]>0.5)
        {
            major[i] = true;
            VariantList[i].MinAlleleString=VariantList[i].refAlleleString;
            VariantList[i].MajAlleleString=VariantList[i].altAlleleString;
        }
        else
        {
            VariantList[i].MinAlleleString=VariantList[i].altAlleleString;
            VariantList[i].MajAlleleString=VariantList[i].refAlleleString;
        }
	}

}

void HaplotypeSet::CalculateGWASOnlyFreq()
{

    AlleleFreq.resize(GWASOnlycounter, 0.0);
    TotalSample.resize(GWASOnlycounter, 0);

    int i,k;
    for(k=0;k<numHaplotypes;k++)
    {
        for (i = 0; i<GWASOnlycounter; i++)
        {
            if(!GWASOnlyMissingSampleUnscaffolded[k][i])
            {
                TotalSample[i]++;
                if(GWASOnlyhaplotypesUnscaffolded[k][i])
                    AlleleFreq[i]++;

            }
        }
    }

    major.resize(GWASOnlycounter, false);

    for (int i = 0; i<GWASOnlycounter; i++)
	{

		AlleleFreq[i]/=(double)TotalSample[i];

		if (AlleleFreq[i]>0.5)
        {
            major[i] = true;
            TypedOnlyVariantList[i].MinAlleleString=TypedOnlyVariantList[i].refAlleleString;
            TypedOnlyVariantList[i].MajAlleleString=TypedOnlyVariantList[i].altAlleleString;
        }
        else
        {
            TypedOnlyVariantList[i].MinAlleleString=TypedOnlyVariantList[i].altAlleleString;
            TypedOnlyVariantList[i].MajAlleleString=TypedOnlyVariantList[i].refAlleleString;
        }
	}


}

bool HaplotypeSet::getScaffoldedHaplotype(int sample,int marker)
{

   if(missing[marker]==true)
        {
            return false;
        }
   else
        return haplotypesUnscaffolded[sample][ScaffoldIndex[marker]];
}


void HaplotypeSet::Create(int index, HaplotypeSet &rHap)
{


    vector<bool> padded(rHap.numMarkers);
    rHap.reconstructHaplotype(padded,index);

    numHaplotypes=1;
    numMarkers=(int)padded.size();
    missing.resize(numMarkers,false);
    ScaffoldIndex.resize(numMarkers);
    vector<bool> tempMissin(numMarkers,false);


    haplotypesUnscaffolded.clear();
    MissingSampleUnscaffolded.clear();

    haplotypesUnscaffolded.push_back(padded);
    MissingSampleUnscaffolded.push_back(tempMissin);


    for(int i=0;i<numMarkers;i++)
        ScaffoldIndex[i]=i;


}


bool HaplotypeSet::getMissingScaffoldedHaplotype(int sample,int marker)
{

   if(missing[marker]==true)
        {
            return true;
        }
   else
        return MissingSampleUnscaffolded[sample][ScaffoldIndex[marker]];
}

