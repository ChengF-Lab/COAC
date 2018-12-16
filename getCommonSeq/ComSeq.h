#define ARMA_64BIT_WORD
#define ARMA_BLAS_LONG_LONG
#include <iterator>
#include <iostream>
#include <armadillo>
#include <stdlib.h>
#include <map>
#include <vector>
#include <list>
#include <set>
#include "assert.h"
#include <algorithm>
#include "string.h"
#include <float.h>
// #include <boost/math/distributions/fisher_f.hpp>
 

using namespace std;
using namespace arma;

 typedef struct PDFS {
	unsigned int ID;	// The ID of the original input sequence
	unsigned int subID; // The location  of the subsequence  in the same sequence
}PDFS;

 typedef struct ResultEntry{
  unsigned int ID;
  unsigned int subID;
  unsigned int length;
 }ResultEntry;
 typedef struct filteKey{
  unsigned int ID;
  unsigned int endsubID;
 }filteKey;

void svd_test_();

void ClusterRead( FILE *ip);

typedef struct DFS {
    unsigned  short label;// Label in the sequence
    unsigned int depth;
    unsigned int support;
    std::vector <PDFS> Projected;
}DFS;

typedef struct DFSR {
    unsigned int support;
    unsigned int index;
    unsigned int depth;
    std::vector <int> Projected;
}DFSR;
typedef struct DFSC{
    unsigned int Rank;
    unsigned int index;
    unsigned int depth;
    std::vector <int> Projected;
}DFSC;



typedef struct Merged_Result
{
    std::vector <int> PatternIDs;
    std::vector <int> Items;

}Merged_Result;

typedef std::map<unsigned  short, std::vector <PDFS> > Projected_map;



typedef struct clusterEntry
{
   int Rank;
   std::vector <int> IDArray;
   std::vector <int> PatternIDs;
}clusterEntry;




typedef struct MotifElem
{
    std::list<std::string>  Motif;
    int Freq;
}MotifElem;


 static const  int ResultEntryCompare(ResultEntry a,ResultEntry b)
{

	if(a.ID<b.ID)
	{
	    return 1;

	}else if(a.ID>b.ID)
	{
	    return  0;

	}else if(a.subID<b.subID)
	{
		return 1;

	}else
	{
	    return  0;
	 }

}
 static const bool merged_asc( const std::vector <ResultEntry>&a, const std::vector <ResultEntry>&b) {
  if(a.size()!=b.size())
  {
          return false;
  }
  for(int i=0;i<a.size();i++)
  {
      if(a[i].ID != b[i].ID||(a[i].subID>(b[i].subID-a[i].length)))
      {
          return false;
      }
  }
         return true;
       }

 static const bool merged_asc_fuzz( const std::vector <ResultEntry>&a, const std::vector <ResultEntry>&b) {

 double cross = 0.0;
 int i =0;
 int j =0;
  for(;i<a.size()&&j<b.size();i++,j++)
  {
      if(a[i].ID == b[j].ID)
      {
          cross ++;

      if (a[i].subID > b[j].subID||(a[i].subID>(b[j].subID-a[i].length)))
      {
          return false;
      }

      }
      if(a[i].ID > b[j].ID)
      {
         j++;
      }
      if(a[i].ID < b[j].ID)
      {
         i++;
      }

  }

        if((cross*cross)/((a.size()*b.size())*1.0)>0.9*0.9)
        {
           return true;
        }else
         {
           return false;
         }
}


struct Compare
{
bool operator()(const std::vector<int> &a, const std::vector<int> &b) const
{

      if(a.size()<b.size())
      {
          return false;
      }
        if(a.size()>b.size())
      {
          return true;
      }
    for(int i=0;i<a.size();i++)
    {
              if (a[i] > b[i])
              {
                 return  true;
              }
              if(a[i] < b[i])
              {
                 return  false;
              }


    }
          return   false;
}
};



struct compare
{

bool operator()( const std::vector <filteKey>&a, const std::vector <filteKey>&b) const
{

      if(a.size()<b.size())
      {
          return false;
      }
        if(a.size()>b.size())
      {
          return true;
      }
    for(int i=0;i<a.size();i++)
    {
              if (a[i].ID > b[i].ID)
              {
                 return  true;
              }
              if(a[i].ID < b[i].ID)
              {
                 return  false;
              }
              if(a[i].ID == b[i].ID)
              {
                  if(a[i].endsubID < b[i].endsubID)
                  {
                      return  true;
                  }
              }
    }
          return   false;
}
};

template <typename T>
std::vector<size_t> ordered(std::vector<T> const& values,bool ascending) {
    std::vector<size_t> indices(values.size());
    std::iota(begin(indices), end(indices), static_cast<size_t>(0));

    if(ascending)
    {
        std::sort(
            begin(indices), end(indices),
            [&](size_t a, size_t b) { return values[a] < values[b]; }
         );
    }else
    {
        std::sort(
            begin(indices), end(indices),
            [&](size_t a, size_t b) { return values[a] > values[b]; }
         );
    };
    return indices;
}

struct DFSPath: public std::vector <DFS> {
public:
	void push (int label,std::vector <PDFS>&Projected,int depth)
	{
		resize (size() + 1);
		DFS &d = (*this)[size()-1];
		d.label = label;
		d.depth = depth;
		d.Projected = Projected;
	}
	void pop () { resize (size()-1); }
};//The  type  of deep first search



struct DFSCPath: public std::vector <DFSC> {
public:
	void push (std::vector <int>&Projected,int index)
	{
		resize (size() + 1);
        DFSC &d = (*this)[size()-1];
        d.Projected = Projected;
		d.index = index;
	}
	void pop () { resize (size()-1); }
};//The  type  of deep first search



struct DFSRPath: public std::vector <DFSR> {
public:
	void push (std::vector <int>&Projected,int Index)
	{
		resize (size() + 1);
		DFSR &d = (*this)[size()-1];
		d.Projected = Projected;
        d.index     = Index;

  }
	void pop () { resize (size()-1); }
};



class ClusterToFeature{

  public:
  std::vector < std::vector <std::string> > FeatureNameArray;
  std::vector < mat> FeatureComponentArray;
  std::map < std::vector <int>,int > ComponentTable;
  std::vector < mat> FeatureM; 
  std::map < std::string ,int > NetworkName;
  std::map < std::vector <int>,std::string > NetworkDictionary;
  mat Sample;
  std::vector<std::string> SampleName;
  std::set<std::string> FilterSet;
  void SampleDecipher(FILE *Op);
  void ClusterRead( FILE *ip);
  void OutPutFeature(FILE *Op);
  void OutPutFeature();
  void SampleRead(FILE *ip);
  void ReadFilteFile(FILE *ip);
  void FilteFeature();
  void ReadFeature(FILE *Op);
  void OutPutFeatureS(FILE *Op);
  void NetWorkFilteFeature();
  void ReadNetWork(FILE *ip);

  mat U;
  mat V;
  vec s;
  double tolerance;
  double CorThreshold=0;
  double CompRThreshold=0;
  void ReadSpareM( std::vector < std::vector < unsigned int > > & Index,std::vector < std::vector < double > > & Value,FILE *Op);
  void StatisticSpareM( std::vector < std::vector < unsigned int > > &Index,std::vector < std::vector < double > > &Value, std::vector < double > &SumValue,std::vector < double > &SquareValue);
 void StatisticSpareMGroup( std::vector < std::vector < unsigned int > > &Index,std::vector < std::vector < double > > &Value, std::vector<unsigned int> &GroupFlag,std::vector < double > &SumValue,std::vector < double > &SquareValue);
void StatisticSpareCor( std::vector < std::vector < unsigned int > > &Index,std::vector < std::vector < double > > &Value, std::vector<double> & CorFlag,std::vector < double > & CorValue);

};








class ComSeq{
private:
    std::vector <unsigned short>  SubSeq;
    DFSPath   DFSpath;
    DFSRPath   DFSRpath;
    DFSCPath  DFSCpath;
    Projected_map  new_projected;
    void  support(DFS &Node);
    int   support(std::vector<ResultEntry> Node);
    void report(std::vector <PDFS>&reportset,int length);

public:
    unsigned int minsupport ;
    int minlength = 4;
    int Threshold = 5;
    int LimitedRank = 1024;
    int MaxDepth = 1024;
    double tolerance;
    double ScannedTolerance=0.1;
    double range;
    double NoZeroRatio=0;
    double Log_Pvalue;
    bool rankSkip = false;
    unsigned int PC_Sup;
    mat Sample;
    mat Restriction;
    bool Centering;
    std::vector<std::string> SampleName;
    std::vector<std::string> RestrictionName;
    mat Correlation_Matrix;
    mat Retrict_Matrix;
    mat U;
    mat V;
    vec s;
    mat KCPath;
    mat cor;
    int Rank =0;
    void RestrictionReadWithName( FILE *ip);
    void RestrictionRead( FILE *ip);
    void ReadFilteFile(FILE *ip);
    void OutCluster(FILE * op);
    void filter();
    std::vector<clusterEntry> clusterResult;
    std::vector<std::vector<int> >CentralStack;
    std::vector<std::vector<double> >CentralVStack;
    double * ReportPathMem;
    double   PC_CutOff;
    double   CutThresHold=0;
    std::map <std::vector <filteKey>,std::vector<ResultEntry>,compare> FiltedResult;
    std::map <int, std::vector< std::vector<int> > >rankMap;
    bool rankFilter(int Rank,std::vector <int> ReportPath,std::vector<int> parentSet);
    void rankInsert(int Rank, std::vector <int> ReportPath);
    std::vector<std::vector<ResultEntry> >TranArray;
    std::vector<std::vector<unsigned short> >Results;
    std::vector<std::vector<unsigned short> >sample;
    std::set < std::vector<int>,Compare>filterset;
    std::vector< std::vector<ResultEntry> >ResultEntryArray;
    std::vector<Merged_Result>Merged_Result_Array;
    std::vector< std::vector<int> >ResultIDArray;
    std::vector< std::vector<int> >Simple_Nozero;
    void matrixRead(FILE*op);
    void SetCutOff(double CutOffInit,float PCSupport);
    int RobustRank(mat matPath);
    int RobustRank(int index,mat matPath);
    void AttributeCluster();
    void fastaread(FILE * ip);
    void read(FILE *ip);
    void output(FILE *op);
    void OutClusterWithName(FILE * op);
    void OutClusterWithNameCombact(FILE *op);
    void run_intern ();
    void KernalCluster();
    void PreAttributeCluster();
    void matrixReadWithName(FILE *ip);
    std::set<std::string> FilterSet;
    std::set<int> FilterSetID;

    std::vector <int> Common_Set(std::vector <int>a,std::vector <int>b);

};
