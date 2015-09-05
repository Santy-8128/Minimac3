#ifndef HAPLOTYPESET_H_INCLUDED
#define HAPLOTYPESET_H_INCLUDED

#include "StringBasics.h"
#include "VcfFileReader.h"
#include "Unique.h"
#include <fstream>
#include <sstream>
#include <algorithm>


using namespace std;



class HaplotypeSet
{

	public:
		int         numHaplotypes,numSamples;
		int         numMarkers;
		bool        unphasedOutput;
		vector<int>        optEndPoints;
		vector<int>        ScaffoldIndex;
		vector<int>        UnScaffoldIndex;
		vector<ReducedHaplotypeInfo> ReducedStructureInfo;
		vector<vector<ReducedHaplotypeInfo> > ReducedStructureInfoBuffer;

		vector<vector<bool> >     haplotypesUnscaffolded;
		vector<vector<bool> >     MissingSampleUnscaffolded;
        vector<vector<bool> >     GWASOnlyhaplotypesUnscaffolded;
		vector<vector<bool> >     GWASOnlyMissingSampleUnscaffolded;



		vector<vector<double> > alleleFreq;
		vector<double> AlleleFreq;
        vector<int> TotalSample;
		vector<vector<float> > Dosage;
		vector<vector<bool> > ImputedAlleles;
		vector<vector<bool> > GWASOnlyAlleles;
		vector<vector<bool> > GWASOnlyMissing;
		vector<int> GWASOnlySampleIndex;
		int PrintStartIndex,PrintEndIndex;
		bool RsId;
		bool GWASPanelOnly;
        bool GT,DS,GP;
        bool onlyCompress,filter;
        bool EstimateNcompress;
        bool PseudoAutosomal;
        bool AllMaleTarget;
        bool Duplicates;
        bool TypedOnly;
        int GWASOnlycounter;
        vector<int> RefTypedIndex;
        int RefTypedTotalCount;
        bool Filter;
        int transFactor, cisFactor;
        string MyChromosome;
        int CPU;

        vector<double>      Recom,Error;
        vector<string> markerName;
		vector<string> individualName;
		vector<int> SampleNoHaplotypes;
		vector<char> refAlleleList;
		vector<bool> major, minor;

        int OrigStartPos;
        int OrigEndPos;


        vector<char> GWASOnlyrefAlleleList;


		vector<bool> missing, MarkerIndices;
		vector<variant> VariantList;
		vector<variant> TypedOnlyVariantList;



        string finChromosome;
		bool allowMissing, vcfType,m3vcfxType,machType;


        HaplotypeSet()
		{
			numHaplotypes = 0;
			numMarkers = 0;
			ReducedStructureInfo.clear();
			alleleFreq.clear();
			individualName.clear();
			SampleNoHaplotypes.clear();


			markerName.clear();
			refAlleleList.clear();
			GWASOnlyrefAlleleList.clear();
			major.clear();
			minor.clear();
			missing.clear();
			MarkerIndices.clear();
			allowMissing = true;
			vcfType = false;
			m3vcfxType=false;
			machType=false;
			PrintStartIndex=0;
			PrintEndIndex=0;
            Recom.clear();
            Error.clear();

		}
        bool    getScaffoldedHaplotype                      (int sample,int marker);
        bool    getMissingScaffoldedHaplotype                      (int sample,int marker);


        void    Create                                      (vector<bool> &tempHaplotype);
        void    calculateFreq                               ();
		void    CalculateFreq                               ();
		void    CalculateGWASOnlyFreq                       ();
        bool    LoadSnpList                                 (String filename);
		bool    LoadMachHaplotypes                          (String filename, String targetSnpfile, vector<string> &refSnpList);
		bool    LoadMachHaplotypes                          (String filename, String targetSnpfile);
        char    convertAlleles                              (string markerId, string indivId, const char *alleles, string refAlleleString, string altAlleleString);
		bool    LoadHaplotypes                              (String filename, String snpNames, int maxIndiv, int maxMarker,bool optmx,String CNO,int START,int END);
		bool    FastLoadHaplotypes                          (String filename, int maxIndiv, int maxMarker,String CNO,int START,int END,int WINDOW,bool rsid,bool compressOnly,bool filter);
		bool    FasterLoadHaplotypes                         (String filename, int maxIndiv, int maxMarker,String CNO,int START,int END,int WINDOW,bool rsid,bool compressOnly,bool filter);
//		bool    LoadTargetHaplotypes                        (String filename, String targetSnpfile, vector<string> &refSnpList, HaplotypeSet &rHap);

        bool    LoadTargetHaplotypes                        (String filename, String targetSnpfile, vector<string> &refSnpList,HaplotypeSet &rHap,bool typedOnly,bool passOnly);
		bool    LoadVcfTargetHaplotypes                     (String filename, String snpNames, vector<string> &refSnpList,HaplotypeSet &rHap);
        void    PrintDosageGWASOnlyForVcfOutputForID        (HaplotypeSet &tHap,
                                                            IFILE vcfdose,
                                                            int MarkerIndex);
        void    PrintDosageGWASOnlyForVcfOutputForIDMaleSamples        (HaplotypeSet &tHap,
                                                            IFILE vcfdose,
                                                            int MarkerIndex);

        void    SaveIndexForGWASOnlyForVcfOutput            (int SamID,int HapId);

        void    convertToReducedStructure                   (vector<int> &optStructure);
        void    writem3vcfFile                              (String &filename,bool &gzip);
        bool    readm3vcfFile                               (String m3vcfFile,String CHR,int START,int END,int WINDOW);
        void    reconstructHaplotype                        (vector<bool> &reHaplotypes,int &index);
        void    SaveDosageForVcfOutput                      (int hapID,vector<float> dose,vector<bool> impAlleles);
        void    SaveDosageForVcfOutputSampleWise            (int SamID,string &SampleName, vector<float> &dose1,vector<float> &dose2,vector<bool> &impAlleles1,vector<bool> &impAlleles2);
        void    InitializeDosageForVcfOutput                (int NHaps,int NMarkers);
        void    InitializePartialDosageForVcfOutput         (int NHaps,int NMarkers, vector<bool> &Format);
        void    InitializePartialDosageForVcfOutputMaleSamples    (int NHaps,int NMarkers, vector<bool> &Format);
        void    PrintDosageForVcfOutputForID                (IFILE vcfdose, int MarkerIndex,bool majorIsReference,char refAllele);
        void    PrintPartialDosageForVcfOutputForID         (IFILE vcfdose, int MarkerIndex,bool majorIsReference,char refAllele);
        bool    BasicCheckForTargetHaplotypes               (String filename);
        string  DetectTargetFileType                        (String filename);
        string  DetectReferenceFileType                     (String filename);
        void    SaveDosageForVcfOutputSampleWiseChrX        (int SamID,string &SampleName, vector<float> &dose1,
                                                            vector<bool> &impAlleles1);
        bool    CheckValidChrom                             (string chr);

        void PrintDosageForVcfOutputForIDMaleSamples(IFILE vcfdose, int MarkerIndex,bool majorIsReference,char refAllele);

        void updateCoeffs(int trans,int cis)
        {
            transFactor = trans;
            cisFactor = cis;

        }


};


#endif // HAPLOTYPESET_H_INCLUDED
