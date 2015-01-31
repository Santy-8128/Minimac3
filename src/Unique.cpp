

#include "Unique.h"
double transFactor = 2.0;
double cisFactor = 3.0;



void findUnique::UpdateDeltaMatrix(vector<String> & haplotypes, vector<int> & index,
          vector<int> & firstDifference, int length, int blockSize,
          vector<int> & oldIndex,  vector<int> & previousPredecessor,  vector<int> & previousDifference)
   {
   int haps = index.size();

//    for (int i = 0; i < oldIndex.size(); i++)
//                {
////                    index[offsets[Haplotypes[oldIndex[i]][length - 1] - '0']++] = oldIndex[i];
//                    cout<<" OLD INDEX [ "<<i<<" = "<<oldIndex[i]<<endl;
//                }

   previousPredecessor[oldIndex[0]] = -1;
//    cout<<" NO = "<<oldIndex[0]<<endl;
   for (int i = 1; i < haps; i++)
      {
      previousPredecessor[oldIndex[i]] = oldIndex[i-1];
      previousDifference[oldIndex[i]]  = firstDifference[i-1];

//      cout<<" CHECK = "<<i<<" "<<oldIndex[i]<<"\t"<<previousPredecessor[oldIndex[i]]
//      <<"\t"<<previousPredecessor.size()<<endl;
      }

   for (int i = 1; i < haps; i++)
      {
      if (index[i-1] == previousPredecessor[index[i]])
         {
         firstDifference[i-1] =
            haplotypes[index[i]][length - 1] ==
            haplotypes[index[i-1]][length - 1] ?
               previousDifference[index[i]] + 1 : 0;
         continue;
         }

      int diff = 0;
      while (diff < length && diff < blockSize &&
               haplotypes[index[i]][length - diff - 1] ==
               haplotypes[index[i-1]][length - diff - 1])
         diff++;

      firstDifference[i - 1] = diff;
      }


   }


void findUnique::AnalyzeBlocks(
         vector<int> & index, vector<int> & firstDifference, int length, int blockSize,
         vector<int> & cost,
         vector<int> & bestSlice, vector<int> & bestComplexity, vector<vector<int> > &bestIndex)
   {
   int haps = index.size();

   // Try to figure out optimal block split
   for (int i = 1; i <= blockSize && i <= length; i++)
      {
      int distinctHaplos = 1;

      for (int j = 0; j < haps - 1; j++)
         if (i > firstDifference[j])
            distinctHaplos++;

      int currentCost=1;

      if(i>1)
        currentCost= transFactor * haps + (i) * distinctHaplos  * cisFactor + cost[length - i+1];

       // cout<<length<<" \t"<<i<<"\t"<<distinctHaplos<<endl;
      if (i==2)
         {
         cost[length] = currentCost;
         bestSlice[length] = 2;
         bestComplexity[length] = distinctHaplos;
         }
      else if (cost[length] > currentCost)
         {
         cost[length] = currentCost;
         bestSlice[length] = i;
         bestComplexity[length] = distinctHaplos;
         }
      else if (cost[length] + transFactor * haps < currentCost)
         break;
      }
    bestIndex[length-1] = index;


//
//    cout<<"\n\n\n CHECK INDEX FOR = "<<length<<endl;
//    for(int i=0;i<index.size();i++)
//        cout<<index[i]<<"\t";
//
//    cout<<endl;





   }


double findUnique::FlushBlocks(vector<int> &optEndPoints,vector<ReducedHaplotypeInfo> &HapInfo, vector<variant> &varList, int LastflushPos,vector<String> & haplotypes, vector<int> & cost,
                   vector<int> & bestComplexity, vector<int> & bestSlice, vector<vector<int> > &bestIndex)
{

    ReducedHaplotypeInfo tempInfo;
    int where   = haplotypes[0].Length();
    where--;
    int newCost = cost[where+1];
    int totalBlocks=0;
    double totalComplexity = 0;
    vector<int> examplars, labels;

    vector<int> blockBoundaries;
    blockBoundaries.clear();

    // Work out optimal set of block boundaries
    //cout<<where+LastflushPos<<endl;

    while (where != 0)
      {
      blockBoundaries.push_back(where);
      //  cout<<where<<"\t"<<bestSlice[where+1]<<"\t";
      totalBlocks++;
      totalComplexity += bestComplexity[where+1];
//      cout<<" WHERE 1 = "<<where<<endl;
      where = where - bestSlice[where+1]+1;
//      cout<<" WHERE 2 = "<<where<<endl;
      };



    int blockStart = 0;
//
    while (blockStart != (haplotypes[0].Length()-1))
    {

        int blockEnd = blockBoundaries.back();
        blockBoundaries.pop_back();
        examplars.clear();
        examplars.push_back(bestIndex[blockEnd][0]);

//
//    cout<<"BEST INDEX = "<<blockEnd<<endl;
//    for(int i=0;i<bestIndex[blockEnd].size();i++)
//        cout<<bestIndex[blockEnd][i]<<"\t";
//
//    cout<<endl;




        labels.resize(bestIndex[blockEnd].size());
        labels[bestIndex[blockEnd][0]] = 0;
        int countCard=0;

        for (int i = 1; i < (int) bestIndex[blockEnd].size(); i++)
        {
            int previous = bestIndex[blockEnd][i-1];
            int current  = bestIndex[blockEnd][i];
            countCard++;
            for (int j = blockStart; j <= blockEnd; j++)
            if (haplotypes[previous][j] != haplotypes[current][j])
               {
               //tempInfo.uniqueCardinality.push_back(countCard);
               //countCard=0;
               examplars.push_back(current);
               break;
               }

         labels[current] = examplars.size() - 1;

         //for(int j=blockStart+LastflushPos; j<blockEnd+LastflushPos; j++)


         //cout<<labels[current]<<"\t";
         }
         //tempInfo.uniqueCardinality.push_back(countCard);

         //cout<<labels.size()<<endl;
        tempInfo.uniqueIndexMap=labels;

        tempInfo.startIndex=blockStart+LastflushPos;
        tempInfo.endIndex=LastflushPos+blockEnd;
        //optEndPoints.push_back(tempInfo.startIndex);
        tempInfo.uniqueCardinality.clear();
        tempInfo.uniqueCardinality.resize(examplars.size(),0);
        for (int i = 0; i < (int)labels.size(); i++)
        {
            tempInfo.uniqueCardinality[labels[i]]++;

        }



        tempInfo.uniqueHaps.resize(examplars.size());
        for (int i = 0; i < (int)examplars.size(); i++)
            tempInfo.uniqueHaps[i].resize(blockEnd-blockStart+1);
        //cout<<tempInfo.startIndex<<"\t"<<tempInfo.endIndex<<endl;

//        for (int i = blockStart; i <= blockEnd; i++)
//            {
//
//
//                //cout<<i+LastflushPos<<"\t"<<varList[i+LastflushPos].name<<"\t";
//                for (int j = 0; j < (int)examplars.size(); j++)
//                    {
//                        tempInfo.uniqueHaps[j][i-blockStart]=haplotypes[examplars[j]][i]=='0'?varList[i+LastflushPos].refAllele:varList[i+LastflushPos].altAllele;
//cout<<(int)tempInfo.uniqueHaps[j][i-blockStart]<<" ";
//            assert((int)tempInfo.uniqueHaps[j][i-blockStart]>0 && (int)tempInfo.uniqueHaps[j][i-blockStart]<8);
//                    }
//            }

        for (int j = 0; j < (int)examplars.size(); j++)
       {

           //cout<<j<<" = "<<"\t";
            for (int i = blockStart; i <= blockEnd; i++)
            {
                tempInfo.uniqueHaps[j][i-blockStart]=haplotypes[examplars[j]][i]=='0'?varList[i+LastflushPos].refAllele:varList[i+LastflushPos].altAllele;

                //cout<<(int)tempInfo.uniqueHaps[j][i-blockStart]<<" ";
            }
            //cout<<endl;
    }

        HapInfo.push_back(tempInfo);
    //cout<<endl;


//  for (int i = blockStart; i <= blockEnd; i++)
//            {
//
//
//                //cout<<i+LastflushPos<<"\t"<<varList[i+LastflushPos].name<<"\t";
//                    {
//                        cout<<(int)tempInfo.uniqueHaps[j][i-blockStart]<<" ";
//            assert((int)tempInfo.uniqueHaps[j][i-blockStart]>0 && (int)tempInfo.uniqueHaps[j][i-blockStart]<8);
//                    }
//            }

         blockStart = blockEnd;

    }

    //optEndPoints.push_back(tempInfo.endIndex);

    //ReducedStructureInfo.push_back(tempInfo);

    return newCost;
    }

//
//
//int findUnique::findMaxBlockLength(int total)
//{
//
//    // formula based on worst case sequential comparison.
//    int check_number=4*(total+1),i;
//    for (i=1;i<total;i++)
//    {
//        if((i*(pow(2,i-1)))>check_number)
//            return (i-1);
//    }
//
//    return i;
//
//}
//
//void findUnique::swap(int &a,int &b)
//{
//    // function for swaping integer values.
//    int temp=a;
//    a=b;
//    b=temp;
//}
//
//
//void findUnique::unique(vector<vector<char> > &hap)
//{
//
//    cout<<"\n Constructing Optimal Allocation of Genomic Segments ... "<<endl;
//
//    cout<<"\n Finding Unique Number of Haplotypes between every pair of positions ... "<<endl<<endl;
//
//	// MTW
//	vector<int> new_permute_index;
//	vector<int> new_group_index;
//	vector<int> final_permute_index;
//
//    for (int endMarker=0;endMarker<(M);endMarker++)
//    {
//
//        if(endMarker%10000==0)
//            cout<<" Comparing Marker "<<endMarker+1<<" with previous Markers ... ["<< endMarker/M<<"\%] "<<endl;
//
//
//        vector<int> permute_index(N,0);
//        for(int indiv=0;indiv<N;indiv++)
//            permute_index[indiv]=indiv;
//
//
//        vector<int> group_index;
//        group_index.push_back(0);
//        group_index.push_back(N);
//        //cout<<endMarker<<endl;
//		// MTW
////		cout <<"START = " << endMarker << "   ";
////		cout << uniqueCount[endMarker].size() << endl;
//
//        uniqueGivenPrevious(hap,endMarker,endMarker,permute_index,group_index, new_permute_index, new_group_index, final_permute_index);
//		// MTW uniqueGivenPrevious(hap, endMarker, endMarker, permute_index, group_index);
////		cout << "END = " << endMarker << "   ";
////		cout << uniqueCount[endMarker].size() << endl;
//
//    }
//
//    return;
//}
//
//
//ReducedHaplotypeInfo findUnique::uniqueGivenPrevious(vector<vector<char> > &hap,int start,int end,vector<int> &prev_permute_index, vector<int> &prev_group_index,int posFirst)
//{
//
//    // use top 'if' statement if fixed block-length is being use, else use next 'if' statement.
//    //if((end-start)>=max_block || start<0 )
////    if(start<0 )
////        {
////            return;
////        }
//
//
//    vector<int> new_permute_index(N,0);
//    for(int indiv=0;indiv<N;indiv++)
//        new_permute_index[indiv]=indiv;
//
//    vector<int> new_group_index;
//    for(int group=1;group<(int)prev_group_index.size();group++)
//    {
//        new_group_index.push_back(prev_group_index[group-1]);
//        int first_alt_pos;
//        int indiv=prev_group_index[group-1]+1;
//
//        while(indiv<prev_group_index[group] && hap[prev_permute_index[indiv]][start]==hap[prev_permute_index[indiv-1]][start])
//            indiv++;
//
//
//        first_alt_pos=indiv;
//        indiv++;
//
//        while(indiv<prev_group_index[group])
//        {
//
//            if(hap[prev_permute_index[new_permute_index[indiv]]][start]!=hap[prev_permute_index[new_permute_index[indiv-1]]][start])
//                {
//                    swap(new_permute_index[indiv],new_permute_index[first_alt_pos]);
//                    first_alt_pos++;
//                }
//            indiv++;
//        }
//
//        if(first_alt_pos<prev_group_index[group])
//            new_group_index.push_back(first_alt_pos);
//
//    }
//
//    new_group_index.push_back(N);
//    vector<int> final_permute_index(N,0);
//
//    for(int indiv=0;indiv<N;indiv++)
//        final_permute_index[indiv]=prev_permute_index[new_permute_index[indiv]];
//
//
//
//    if(start==0 || start==posFirst)
//    {
//        ReducedHaplotypeInfo temp;
//
//        temp.startIndex=start;
//        temp.endIndex=end;
//        vector<int> tempRep;
//        tempRep.clear();
//
//        for(int i=0;i<(int)new_group_index.size()-1;i++)
//        {
//            tempRep.push_back(final_permute_index[new_group_index[i]]);
//
//
//            temp.uniqueCardinality.push_back(new_group_index[i+1]-new_group_index[i]);
//        }
//
//
//        temp.uniqueIndexMap.resize(N);
//        for(int i=1;i<(int)new_group_index.size();i++)
//        {
//            for(int j=new_group_index[i-1];j<new_group_index[i];j++)
//                temp.uniqueIndexMap[final_permute_index[j]]=i-1;
//        }
////
//        temp.uniqueHaps.resize(tempRep.size());
//        for(int i=0;i<(int)tempRep.size();i++)
//            {
//
//                temp.uniqueHaps[i].resize(end-start+1);
//                for(int j=start;j<=end;j++)
//                    temp.uniqueHaps[i][j-start]=hap[tempRep[i]][j];
//            }
//
//
//        return temp;
//    }
//    else
//        return uniqueGivenPrevious(hap,start-1,end,final_permute_index,new_group_index,posFirst);
//}
//
//
//
//
//
//void findUnique::uniqueGivenPrevious(vector<vector<char> > &hap, int start, int end, vector<int> &prev_permute_index, vector<int> &prev_group_index,
//	vector<int>& new_permute_index, vector<int>& new_group_index, vector<int>& final_permute_index)
//{
//
//    // use top 'if' statement if fixed block-length is being use, else use next 'if' statement.
//    if((end-start)>=max_block || start<0 )
//    //if(start<0 )
//    {
//        return;
//    }
//    if(prev_group_index.size()==(hap.size()+1))
//    {
//        return;
//
//    }
//
//
//	new_permute_index.resize(N, 0);
//// MTW    vector<int> new_permute_index(N,0);
//	for (int Indiv = 0; Indiv<N; Indiv++)
//		new_permute_index[Indiv] = Indiv;
//
//	new_group_index.clear();
//// MTW    vector<int> new_group_index;
//	//int indiv;
//    for(int group=1;group<(int)prev_group_index.size();group++)
//    {
//        new_group_index.push_back(prev_group_index[group-1]);
//        int first_alt_pos;
//		int indiv=prev_group_index[group-1]+1;
//		//cout << indiv << "\t";
//
//        while(indiv<prev_group_index[group] && hap[prev_permute_index[indiv]][start]==hap[prev_permute_index[indiv-1]][start])
//            indiv++;
//
//
//        first_alt_pos=indiv;
//        indiv++;
//
//        while(indiv<prev_group_index[group])
//        {
//
//            if(hap[prev_permute_index[new_permute_index[indiv]]][start]!=hap[prev_permute_index[new_permute_index[indiv-1]]][start])
//                {
//                    swap(new_permute_index[indiv],new_permute_index[first_alt_pos]);
//                    first_alt_pos++;
//                }
//            indiv++;
//        }
//
//        if(first_alt_pos<prev_group_index[group])
//            new_group_index.push_back(first_alt_pos);
//
//    }
//
//    new_group_index.push_back(N);
//	final_permute_index.resize(N, 0);
//	// MTW    vector<int> final_permute_index(N,0);
//
//    for(int indiv=0;indiv<N;indiv++)
//        final_permute_index[indiv]=prev_permute_index[new_permute_index[indiv]];
//
//
//    // use top formula if NOT pre-allocated (slower), if pre-allocated use next form (faster)
//    //uniqueCount[end].push_back((int)new_group_index.size()-1);
//    uniqueCount[end][end-start]=(int)new_group_index.size()-1;
//
//    int mid_point;
//
//
//    // recursive function to calculate the minimum cost and minimum allocation
//    if(start<end)
//    {
//        if(start==(end-1))
//            minCost[end]=minCost[start]+(2*N)+(2*uniqueCount[end][1]);
//        else if (minCost[end]>(minCost[start]+(2*N)+((end-start+1)*uniqueCount[end][end-start])))
//        {
//            minCost[end]=minCost[start]+(2*N)+((end-start+1)*uniqueCount[end][end-start]);
//            minAllocation[end]=start;
//        }
//
//    // check sequential comparison for variable block length selection
//
//    mid_point=(start+end)/2;
//
//    if((2*N)<( ( (end-start+1) * (uniqueCount[end][end-start]) ) - ( 2 * ( (end-mid_point+1) * (uniqueCount[end][end-mid_point] ) ) ) ) )
//    //if((2*N)<(((end-start)*(uniqueCount[end][end-start]-uniqueCount[end][end-start-1]))-(2*uniqueCount[end][1])+(uniqueCount[end][end-start])))
//         return;
//    }
//
//
//
//
//
//
//
//    // recursive function to go back a step previous.
//	uniqueGivenPrevious(hap, start - 1, end, final_permute_index, new_group_index, prev_permute_index, prev_group_index, new_permute_index);
//	// MTW uniqueGivenPrevious(hap, start - 1, end, final_permute_index, new_group_index);
//    return;
//}
//
//
//void findUnique::printAllocation(int position)
//{
//
//    if(position==0)
//    {
//        optimalAllocation.push_back(position);
//    }
//    else
//    {
//        printAllocation(minAllocation[position]);
//        optimalAllocation.push_back(position);
//    }
//
//    return;
//
//
//}
//
//
//void findUnique::createOutput(String filename)
//{
//
//    printAllocation(M-1);
//
//
//
//
//
//
//    ofstream myfile2;
//    myfile2.open(filename+".UniqueCount");
//
//    for(int i=0;i<(int)uniqueCount.size();i++)
//    {
//        for(int j=0;j<(int)uniqueCount[i].size();j++)
//            myfile2<<uniqueCount[i][j]<<"\t";
//        myfile2<<"\n";
//    }
//
//}
//
//
//
//
