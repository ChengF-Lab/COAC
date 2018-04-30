#include "ComSeq.h"
#include "stdio.h"
#include "math.h"
#include <fstream>

using namespace std;
using namespace arma;




void ComSeq::RestrictionRead( FILE *ip)
{
      std::vector<double>Temp;
      char * seqName[1024];
      unsigned int entry = 1;
      int N =0;
      int M =0;

       do{

            double value;
			fscanf(ip,"%[\b|\t]*",seqName);
			fscanf(ip,"%s",seqName);
			N++;
            fscanf(ip,"%[\b|\t]*",seqName);
	} while(!fscanf(ip,"%[\n]",seqName));
    rewind(ip);

      while(!feof(ip))
	  {
		  entry = fgetc(ip);
		  if(entry=='\n')
		  {
		      M++;
		  }
	   }
      rewind(ip);
      Restriction.set_size(M,N);
      for(int i=0;i<M;i++)
	  {
         for(int j=0;j<N;j++)
		 {
			    double value;
              	fscanf(ip,"{\b|\t|\n}*",seqName);
			    fscanf(ip,"%lf",&value);
				Restriction(i,j)=value;
		 }
	 }
}

void ComSeq::RestrictionReadWithName( FILE *ip)
{

      std::vector<double>Temp;
      char * seqName[1024];
      char buffer[1024];
      unsigned int entry = 1;
      int N =0;
      int M =0;

       do{

            double value;
			fscanf(ip,"%[\b|\t]*",seqName);
			fscanf(ip,"%s",seqName);
			N++;
            fscanf(ip,"%[\b|\t]*",seqName);
	} while(!fscanf(ip,"%[\n]",seqName));
    rewind(ip);

      while(!feof(ip))
	  {
		  entry = fgetc(ip);
		  if(entry=='\n')
		  {
		      M++;
		  }
	   }
      rewind(ip);
      Restriction.set_size(M-1,N);
      RestrictionName.resize(N);

      for(int j=0;j<N;j++)
      {
         fscanf(ip,"{\b|\t|\n}*",buffer);
         fscanf(ip,"%s",buffer);
         RestrictionName[j].assign(buffer);
      }

      for(int i=0;i<M-1;i++)
	  {
         for(int j=0;j<N;j++)
		 {
			    double value;
              	fscanf(ip,"{\b|\t|\n}*",seqName);
			    fscanf(ip,"%lf",&value);
				Restriction(i,j)=value;
		 }
	 }


}



void ComSeq::matrixRead(FILE *ip)
{
  std::vector<double>Temp;
  char * seqName[1024];
  unsigned int entry = 1;
  int N =0;
  int M =0;

    do{

            double value;
			fscanf(ip,"%[\b|\t]*",seqName);
			fscanf(ip,"%lf",&value);
			N++;
            fscanf(ip,"%[\b|\t]*",seqName);
	} while(!fscanf(ip,"%[\n]",seqName));
rewind(ip);
     while(!feof(ip))
	{
		  entry = fgetc(ip);
		  if(entry=='\n')
		  {
		      M++;
		  }
	 }
   rewind(ip);
   Sample.set_size(M,N);
	 for(int i=0;i<M;i++)
	 {


         for(int j=0;j<N;j++)
		 {
			    double value;
              	fscanf(ip,"{\b|\t|\n}*",seqName);
			    fscanf(ip,"%lf",&value);
				Sample(i,j)=value;
		 }
	 }
     if(Centering)
     {
        for(int j=0;j<N;j++)
        {
            double mean = 0.0;

            for(int i=0;i<M;i++)
            {
                mean = mean +Sample(i,j);
            }
            mean = mean/M;
         for(int i=0;i<M;i++)
         {
            Sample(i,j) = Sample(i,j) - mean;
         }
        }
     }
}




void ComSeq::matrixReadWithName(FILE *ip)
{
  std::vector<double>Temp;
  char buffer[1024];
  unsigned int entry = 1;
  int N =0;
  int M =0;

   std::string tempword;

     do{


	    fscanf(ip,"%[\b|\t]*",buffer);
	    fscanf(ip,"%s",&buffer);
            N++;
            fscanf(ip,"%[\b|\t]*",buffer);
	} while(!fscanf(ip,"%[\n]",buffer));

     rewind(ip);
     while(!feof(ip))
	{
		  entry = fgetc(ip);
		  if(entry=='\n')
		  {
		      M++;
		  }
	 }
   rewind(ip);
   Sample.set_size(M-1,N);
   SampleName.resize(N);


    for(int j=0;j<N;j++)
    {
              	fscanf(ip,"{\b|\t|\n}*",buffer);
			    fscanf(ip,"%s",buffer);
			    SampleName[j].assign(buffer);
    }
	 for(int i=0;i<M-1;i++)
	 {	 for(int j=0;j<N;j++)
		 {
			    double value;
              	fscanf(ip,"{\b|\t|\n}*",buffer);
			    fscanf(ip,"%lf",&value);
				Sample(i,j)=value;
		 }
	 }

   if(Centering)
 {
        for(int j=0;j<N;j++)
        {
            double mean = 0.0;

            for(int i=0;i<M-1;i++)
            {
                mean = mean +Sample(i,j);
            }

                mean = mean/(M-1);

            for(int i=0;i<M-1;i++)
            {
                Sample(i,j) = Sample(i,j) - mean;
            }

   }

}

}




 std::vector<std::string>Sample;

void ComSeq::fastaread(FILE * ip)
{
     char * seqName[1024];
     std::vector<unsigned short int>*Temp;
     unsigned int entry = 1;

     while(!feof(ip))
	{

        sample.resize(sample.size()+1);
        Temp = &sample[sample.size()-1];

		fscanf(ip,"%[\b|\t]*",seqName);
	    fscanf(ip,"%s",seqName);
        printf("%s\n",seqName);
		fscanf(ip,"{\b|\t|\n}*",seqName);
        fscanf(ip,"%[\b|\t]*",seqName);
        entry = fgetc(ip);

      do{
          Temp->resize(Temp->size()+1);
          (*Temp)[Temp->size()-1] = entry;
          fscanf(ip,"{\b|\t}*",seqName);
          entry = fgetc(ip);
      } while(entry!='\n'&&entry!=EOF);

      Sample.resize(Sample.size()+1);


	}
}

void ComSeq::read(FILE * ip)
{
    std::vector<unsigned short int>*Temp;
    char buffer[1024];
    unsigned int entry = 1;

    while(entry!=EOF)
    {
        sample.resize(sample.size()+1);
        Temp = &sample[sample.size()-1];
        fscanf(ip,"%[\b|\t|\n]*",buffer);
        entry = fgetc(ip);
     do{
          Temp->resize(Temp->size()+1);
          (*Temp)[Temp->size()-1] = entry;
          fscanf(ip,"{\b|\t}*",buffer);
          entry = fgetc(ip);
      } while(entry!='\n'&&entry!=EOF);
    };
}
int ComSeq::support(std::vector<ResultEntry> Node)
{    unsigned int  oldID = 0xffffffff;
     int support =0;

     for (unsigned int n = 0;n<Node.size();n++){
          ResultEntry *cur = &Node[n];
          if(cur->ID!=oldID){
               support++;
          }
          oldID = cur->ID;
     }
     return support;

}

void ComSeq::support(DFS &Node)
{    unsigned int  oldID = 0xffffffff;

     Node.support =0;

     for (unsigned int n = 0;n<Node.Projected.size();n++){
          PDFS *cur = &Node.Projected[n];
          if(cur->ID!=oldID){
             Node.support++;
          }
          oldID = cur->ID;
     }
}

 void printfvector(std::vector <unsigned short>&SubSeq)
 {
   for(int i=0;i<SubSeq.size();i++)
   {
        printf("%d ",SubSeq[i]);
   }
     printf("\n");
 }

 void printfvector( std::vector<PDFS> &SubSeq)
 {
   for(int i=0;i<SubSeq.size();i++)
   {
        printf("%d-",SubSeq[i].ID);
        printf("%d ",SubSeq[i].subID);
   }
     printf("\n");
 }
void ComSeq::report(std::vector <PDFS>&reportset,int length)
{
   std::vector<ResultEntry>subset;
   unsigned int  oldID = 0xffffffff;

   for (std::vector <PDFS>::iterator cur = reportset.begin(); cur != reportset.end(); ++cur)
   {

            subset.resize(subset.size()+1);
            subset[subset.size()-1].ID = (*cur).ID;
            subset[subset.size()-1].subID = (*cur).subID;
            subset[subset.size()-1].length = SubSeq.size();
   }

       std::sort(subset.begin(),subset.end(),ResultEntryCompare);

       TranArray.resize(TranArray.size()+1);
       TranArray[TranArray.size()-1]=subset;

}

std::vector <int> ComSeq::Common_Set(std::vector <int>a,std::vector <int>b)
{

        std::vector <int> CoVector;

        for(int i=0,j=0;i<a.size()&&j<b.size();)
        {
            if(a[i]==b[j])
            {

               CoVector.resize(CoVector.size()+1);
               CoVector[CoVector.size()-1]= a[i];
               i++;
               j++;
               continue;
            }
            if(a[i]>b[j])
            {
                j++;
            }else
            {
                i++;
            }

        }
        return  CoVector;

    }



void ComSeq::merged_report()
{

    std::vector <ResultEntry> old_one;
    ResultEntryArray.resize(FiltedResult.size());
    ResultIDArray.resize(FiltedResult.size());

     int i =0;
     for(std::map <std::vector <filteKey>,std::vector<ResultEntry>,compare>::iterator cur = FiltedResult.begin();cur != FiltedResult.end(); ++cur)
    {
          ResultEntryArray[i]=cur->second;
          i++;
     }
    for(int i=0;i<ResultEntryArray.size();i++)
    {
        if(ResultEntryArray[i][0].length<2)
            continue;
        int Num = support(ResultEntryArray[i]);
        ResultIDArray[i].resize(Num);
      {
            unsigned int  oldID = 0xffffffff;
            int j =0;
            for (unsigned int n = 0;n<ResultEntryArray[i].size();n++){
               ResultEntry *cur = &ResultEntryArray[i][n];
               if(cur->ID!=oldID){
                    ResultIDArray[i][j] = cur->ID;
                    j++;
                  }
           oldID = cur->ID;
       }
    }
    }
    associte();
    std::list<std::string> Motif;
    std::vector <unsigned short>  merged_result;
    std::vector <MotifElem>Motiflist;
    std::vector <ResultEntry> Temp;
    std::string MotifEntry;

    for(int i=0;i<Merged_Result_Array.size();i++)
    {

        Motif.resize(0);
        Motiflist.resize(0);

        for(int j=0;j<Merged_Result_Array[i].Items.size();j++)
        {
            int SampleID = Merged_Result_Array[i].Items[j];
            Temp.resize(0);
             for(int k=0;k<Merged_Result_Array[i].PatternIDs.size();k++)
             {
                  int Pindex = Merged_Result_Array[i].PatternIDs[k];
                  for(int l1 = 0;l1<ResultEntryArray[Pindex].size();l1++)
                  {
                      if (ResultEntryArray[Pindex][l1].ID==SampleID&&ResultEntryArray[Pindex][l1].length>1)
                      {
                            Temp.push_back(ResultEntryArray[Pindex][l1]);
                      }

                  }
                 // printf("%d ",Pindex);
             }
                 // printf("\n");
               if(Temp.size()>1)
               {
                   int lower_bound_;
                   int up_bound_;
                   MotifEntry.resize(0);
                   Motif.resize(0);
                   std:sort(Temp.begin(),Temp.end(),ResultEntryCompare);
                   up_bound_ = lower_bound_ =-1;
                   for(int i=0;i<Temp.size();i++)
                   {
                       if( up_bound_ == -1)
                       {
                           lower_bound_ = Temp[i].subID;
                           up_bound_    = Temp[i].subID + Temp[i].length-1;
                           continue;
                       }
                       if( Temp[i].subID <= up_bound_+1)
                       {
                            int EndIndex = Temp[i].subID + Temp[i].length-1;
                            if(EndIndex>up_bound_)
                            {
                                up_bound_ = EndIndex;
                            }
                            continue;
                       }
                        const char *debug;
                       if( Temp[i].subID > up_bound_+1)
                       {

                           for(int j=lower_bound_;j<=up_bound_;j++)
                            {

                                MotifEntry.append((const char*)&sample[SampleID][j]);

                            }
                              const char a = '\0';
                              MotifEntry.append(&a);

                              printf("%s .",MotifEntry.c_str());
                              Motif.push_back(MotifEntry);
                              MotifEntry.resize(0);

                           lower_bound_ = Temp[i].subID;
                           up_bound_    = Temp[i].subID + Temp[i].length-1;
                           continue;
                       }
                   }

                        {
                            for(int j=lower_bound_;j<=up_bound_;j++)
                            {
                                MotifEntry.append((const char*)&sample[SampleID][j]);
                                }
                                const char a = '\0';
                                MotifEntry.append(&a);
                                printf("%s .",MotifEntry.c_str());
                                Motif.push_back(MotifEntry);
                                MotifEntry.resize(0);
                           }

                  bool Insert = true;
                  int freq =1;


                  for(int i=0;i<Motiflist.size();i++)
                {

                    std::list<std::string>::iterator it1;
                    std::list<std::string>::iterator it2;

                    bool Mismatch = false;
                    for (it1 =  Motiflist[i].Motif.begin(),it2 = Motif.begin(); it1 != Motiflist[i].Motif.end()&&it2!= Motif.end(); it1++,it2++)
                    {       std::string temp1;
                            std::string temp2;
                            temp1 = (*it1);
                            temp2 = (*it2);
                           if(!(temp1.compare(temp2)==0))
                           {
                               Mismatch = true;
                               break;
                           }
                    }
                    int aa1 = Motiflist[i].Motif.size();
                    int bb1 = Motif.size();


                    if(it1 != Motiflist[i].Motif.end()||it2!= Motif.end())
                    {
                       Mismatch = true;
                    }
                    if(Mismatch)
                    {
                            std::string temp2;
                            std::string temp1;

                            it2 = Motif.begin();
                            temp2 = (*it2);
                            const  char * str = temp2.c_str();
                            for (it1 = Motiflist[i].Motif.begin(); it1 != Motiflist[i].Motif.end()&&it2!=Motif.end(); it1++)
                            {
                                      temp1 = (*it1);
                                      while(it2!= Motif.end())
                                     {

                                         str = strstr(str,temp1.c_str());

                                         if(str)
                                         {

                                            str = str + strlen(temp1.c_str());

                                            break;
                                         }else
                                         {
                                            it2++;
                                            if(it2!= Motif.end())
                                           {
                                               temp2 = (*it2);
                                               str = temp2.c_str();
                                           }

                                         }
                                     }
                                  if(it2 == Motif.end())
                                 {
                                    break;
                                 }

                            }
                            if(it1 == Motiflist[i].Motif.end())
                            {
                                Motiflist[i].Freq = Motiflist[i].Freq+1;
                            }



                            it2 = Motiflist[i].Motif.begin();
                            temp2 = (*it2);
                            str = temp2.c_str();
                            for (it1 = Motif.begin(); it1 != Motif.end()&&it2!=Motiflist[i].Motif.end(); it1++)
                            {
                                      temp1 = (*it1);
                                      while(it2!= Motiflist[i].Motif.end())
                                     {

                                           str = strstr(str,temp1.c_str());



                                         if(str)
                                         {
                                            str = str + strlen(temp1.c_str());

                                            break;
                                         }else
                                         {
                                           it2++;
                                           if(it2!= Motiflist[i].Motif.end())
                                           {
                                               temp2 = (*it2);
                                               str = temp2.c_str();
                                           }

                                         }
                                     }
                                 if(it2 == Motiflist[i].Motif.end())
                                 {
                                    break;
                                 }

                            }
                            if(it1 == Motif.end())
                            {
                                 freq = freq+1;
                            }

                    }else{

                          Motiflist[i].Freq = Motiflist[i].Freq+1;
                          Insert = false;
                    }
                  }
                  if(Insert)
                  {
                      MotifElem Entry;
                      Entry.Motif = Motif;
                      int bb1 = Motif.size();
                      Entry.Freq  = freq;
                      Motiflist.push_back(Entry);
                  }
                  printf("\n");
                }
            }


      for(int i=0;i<Motiflist.size();i++)
      {
        MotifElem One;
        One = Motiflist[i];
        merged_result.resize(0);
        if(One.Freq > Threshold)
        {
           std::list<std::string>::iterator it;
           for (it = One.Motif.begin();it!=One.Motif.end();it++)
           {
               std::string temp;
               temp = (*it);
               const char * str = temp.c_str();
               for(int i=0;i<temp.length();i++)
              {
                 merged_result.resize(merged_result.size()+1);
                 merged_result[merged_result.size()-1] = str[i];
               }
                 merged_result.resize(merged_result.size()+1);
                 merged_result[merged_result.size()-1]='.';
           }
            if(std::find(Results.begin(),Results.end(),merged_result)==Results.end())
            {
                Results.resize(Results.size()+1);
                Results[Results.size()-1]=merged_result;
            }
        }
      }
      }
    }


bool ComSeq::rankFilter(int Rank, std::vector <int> ReportPath,std::vector<int> parentSet)
{

     if(rankMap.find(Rank)==rankMap.end())
   {
         return false;
   }else
   {
             if(filterset.find(parentSet) == filterset.end())
            {
                return false;
            }
            std::vector < std::vector <int> > ReportPathArray;
            ReportPathArray = rankMap[Rank];
            for(int i=0;i<ReportPathArray.size();i++)
            {
                int m = 0;
                for( int n=0;n<ReportPathArray[i].size()&&m<ReportPath.size();)
                {
                    if(ReportPath[m] == ReportPathArray[i][n])
                    {
                            n++;
                            m++;
                            continue;
                    }
                    if(ReportPath[m] < ReportPathArray[i][n])
                    {
                            n++;
                            continue;
                    }
                 if(ReportPath[m] >ReportPathArray[i][n])
                 {
                     break;
                 }
              }
              if(m ==ReportPath.size())
              {
                 return true;
              }
          }
   }
                return false;
}


void ComSeq::rankInsert(int Rank, std::vector <int> ReportPath)
{


   if(rankMap.find(Rank)==rankMap.end())
   {
        std::vector < std::vector <int> > ReportPathArray;
        ReportPathArray.push_back(ReportPath);
        rankMap.insert(make_pair(Rank,ReportPathArray));
   }else
   {
        rankMap[Rank].push_back(ReportPath);
   }

}
void ComSeq::SetCutOff(double CutOffInit,float PCSupport)
{
    PC_CutOff =   CutOffInit*(1/sqrt((double)Sample.n_rows));
    PC_Sup     =   (unsigned int)(PCSupport*((float)Sample.n_rows));


}

int ComSeq::RobustRank(mat matPath)
{
     double PC =0;// principal components
     double AC =0;// all  components
     int Rank =0;


     vec s = eig_sym(matPath);
     for(int i=0;i<s.n_elem;i++)
     {
       AC = AC+s[i];
     }

      for(int i=s.n_elem-1;i>=0;i--)
     {
        PC = PC+s[i];
        Rank ++;

        if(PC>((1.0-tolerance)*(1.0-tolerance)*AC))
        {
           break;
        }
    }

    double n1 = sqrt(matPath(matPath.n_rows-1,matPath.n_rows-1));
    double n2;
    double product;


    for(int i=0;i<matPath.n_rows-1;i++)
    {
       n2 = sqrt(matPath(i,i));
       product = matPath(i,matPath.n_rows-1);

      if(abs(product/(n1*n2))<range)
      {
               return -1;
      }

    }
  if( LimitedRank < Rank )
        return -1;
    else
        return Rank;

}

int  ComSeq::RobustRank(int index,mat matPath)
{

  double PC =0;// principal components
  double AC =0;// all the components
  int Rank =0;
  ofstream  U_d;
  ofstream  s_d;
  ofstream  V_d;
  ofstream  M_d;

   //M_d.open("M.txt");
   //V_d.open("V.txt");
  // s_d.open("S.txt");
  // U_d.open("U.txt");



  matPath.resize(matPath.n_rows,matPath.n_cols+1);
  matPath.col(matPath.n_cols-1) = Sample.col(index);
  svd_econ(U,s,V,matPath);
   //  M_d<<matPath;
  //  U_d<<U;
  //  s_d<<s;
  //  V_d<<V;
// M_d.close();
 // U_d.close();
 // s_d.close();
 // V_d.close();


  assert(V.n_rows==s.n_elem);

	for(int i=0;i<s.n_elem;i++)
    {
       AC = AC+s[i]*s[i];
    }
    for(int i=0;i<s.n_elem;i++)
    {
       PC = PC+s[i]*s[i];
       Rank ++;

       if(PC>((1.0-tolerance)*(1.0-tolerance)*AC))
       {
           break;
       }
    }
    int RankAdjusted = Rank;
    for(int i=0;i<Rank;i++)
    {
         assert(U.n_rows == matPath.n_rows);
         int Below_CutOff_Count = 0;
         for(int j=0;j<U.n_rows;j++)
         {
               if(abs(U(j,i))<PC_CutOff)
                Below_CutOff_Count++;
         }
         if(Below_CutOff_Count>PC_Sup)
         {
             RankAdjusted--;
         }
    }



     double AverageCorrelation =0;

     for(int j=0;j<V.n_rows-1;j++)
     {
        double product =0;
        double n1 = 0;
        double n2 = 0;

           for(int i =0;i<Rank;i++)
          {
              product = product +V(j,i)*V(V.n_rows-1,i)*s[i]*s[i];
              n1 = n1 + V(j,i)*V(j,i)*s[i]*s[i];
              n2 = n2 + V(V.n_rows-1,i)*V(V.n_rows-1,i)*s[i]*s[i];
          }

          AverageCorrelation  = product/sqrt(n1*n2);
         if(abs(AverageCorrelation)<range)
         {
              matPath.resize(matPath.n_rows,matPath.n_cols-1);
             return -1;
         }
     }
       matPath.resize(matPath.n_rows,matPath.n_cols-1);
      {
          if( LimitedRank < RankAdjusted )
            return -1;
          else
            return RankAdjusted;

      };
}

void ComSeq::KernalCluster()
{

   PreAttributeCluster();


   for(int i=0;i<Correlation_Matrix.n_cols;i++)
   {
          DFSCpath.push(ResultIDArray[i],i);
          DFSCpath[DFSCpath.size()-1].Rank  = 1;
          DFSCpath[DFSCpath.size()-1].depth = 0;
   }
   {
          std::vector <int>reportPath;
          uvec PathIndices;
          reportPath.resize(0);
          while(DFSCpath.size())
         {
            DFSC * DFSset = &DFSCpath[DFSCpath.size()-1];
            int order = DFSset->index;
            int old_rank = DFSset->Rank;
            int old_depth = DFSset->depth;
            std::vector<int> parentSet;
            std::vector<int> CoVector;
            reportPath.resize(DFSset->depth+1);
            reportPath[DFSset->depth]  = DFSset->index;
            parentSet.resize(DFSset->Projected.size());
            std::copy(DFSset->Projected.begin(),DFSset->Projected.end(),parentSet.begin());
            DFSCpath.pop();
            if(rankFilter(old_rank,reportPath,parentSet))
            {
                 continue;
             }
            int MinRank = old_rank;
            bool Isminimum = true;
            int max_support = 0;

            for(int i=0;i<order;i++)
            {
                  CoVector = Common_Set(ResultIDArray[i],parentSet);
                  int Rank;

                  PathIndices.resize(reportPath.size()+1);
                  int k;
                  for( k =0;k<reportPath.size();k++)
                  {

                      PathIndices(k) = reportPath[k];

                  }
                  PathIndices(k) = i;

                  if(rankSkip)
                  {     int k;
                        for(k=0;k<reportPath.size();k++)
                        {
                            if(Correlation_Matrix(reportPath[k],i)<range)
                            {
                                Rank = -1;
                                break;
                            }
                        }
                        if(k==reportPath.size())
                        {
                            Rank = 1;
                        }
                  }else{

                        KCPath =  Correlation_Matrix.submat(PathIndices,PathIndices);
                        Rank = RobustRank(KCPath);
                  }

                  if(max_support<CoVector.size())
                  {
                     max_support=CoVector.size();
                  }
                  if(Rank>0)
                  {
                        DFSCpath.push(CoVector,i);
                        DFSCpath[DFSCpath.size()-1].depth = old_depth+1;
                        DFSCpath[DFSCpath.size()-1].Rank = Rank;
                        if(MinRank >= Rank)
                        {
                           Isminimum = false;
                        }
                  }
            }

            if(max_support<parentSet.size()||order==0)
            {
                filterset.insert(parentSet);
            }



            if(Isminimum&&(!(rankFilter(old_rank,reportPath,parentSet))))
            {
                    rankInsert(old_rank,reportPath);
                    clusterEntry  Temp;
                    Temp.IDArray = reportPath;
                    Temp.Rank = old_rank;
                    clusterResult.push_back(Temp);
            }

         }
   }


}



void ComSeq::PreAttributeCluster()
{

    double PC =0;// Principal components
    double AC =0;// All the components
    int Rank =0;
    mat U;
    vec s;
    mat V;
    mat cor;


    svd_econ(U,s,V,Sample);

    for(int i=0;i<s.n_elem;i++)
    {
       AC = AC+s[i]*s[i];
    }

    for(int i=0;i<s.n_elem;i++)
    {
       PC = PC+s[i]*s[i];
       Rank ++;
       if(PC>((1.0-tolerance)*(1.0-tolerance)*AC))
       {
           break;
       }
    }

    U.resize(U.n_rows,Rank);
    cor= Sample.t()*U;


    for(int j=0;j< cor.n_cols;j++)
    {
      for(int i=0;i< cor.n_rows;i++)
      {
            cor(i,j) = cor(i,j)/sqrt(Correlation_Matrix(i,i));
      }
   }



     ResultIDArray.resize(0);
     ResultIDArray.resize(Sample.n_cols);


     for(int i=0;i< ResultIDArray.size();i++)
     {
        for(int j=0; j< Rank;j++)
          {
                if(cor(i,j)>= (1.0/sqrt(s.n_elem))*DimRatio)
                {

                     ResultIDArray[i].resize(ResultIDArray[i].size()+1);
                     ResultIDArray[i][ResultIDArray[i].size()-1] = j;
                }
                if(cor(i,j) <= -(1.0/sqrt(s.n_elem))*DimRatio)
                {
                     ResultIDArray[i].resize(ResultIDArray[i].size()+1);
                     ResultIDArray[i][ResultIDArray[i].size()-1] = j+Rank;
                }
                std:sort(ResultIDArray[i].begin(),ResultIDArray[i].end());
          }
     }
}

void ComSeq::AttributeCluster()
{
      DFSCpath.resize(0);


      Correlation_Matrix =  Sample.t()*Sample;
      assert(Correlation_Matrix.n_rows==Correlation_Matrix.n_cols);
      assert(Correlation_Matrix.n_rows==Sample.n_cols);

      for(int i=0;i<Correlation_Matrix.n_rows;i++)
      {
         for(int j=0;j<Correlation_Matrix.n_cols;j++)
         {
               if(i==j)
                 continue;
                Correlation_Matrix(i,j) = Correlation_Matrix(i,j)/sqrt(Correlation_Matrix(i,i)*Correlation_Matrix(j,j));
          }
      }

        if(Restriction.n_rows!=0&&Restriction.n_cols!=0)
       {
           Retrict_Matrix =  Sample.t()*Restriction;
       }


        for(int j=0;j< Retrict_Matrix.n_cols;j++)
        {
             mat  Rn = Restriction.col(j).t()*Restriction.col(j);


             for(int i=0;i< Retrict_Matrix.n_rows;i++)
             {
               Retrict_Matrix(i,j) = Retrict_Matrix(i,j)/sqrt(Correlation_Matrix(i,i)*(Rn(0,0)));
             }
        }

   //      Retrict_Matrix.print("EE:");
         PreAttributeCluster();

         for(int i=0;i<Correlation_Matrix.n_rows;i++)
        {

            Correlation_Matrix(i,i) = 1.0;

        }


      for(int i=0;i<Sample.n_cols;i++)
     {
            DFSCpath.push(ResultIDArray[i],i);
            DFSCpath[DFSCpath.size()-1].Rank  = 1;
            DFSCpath[DFSCpath.size()-1].depth = 0;
     }

     {
         std::vector <int>reportPath;
         mat matPath;//(ReportPathMem,Sample.n_rows,120,false,true);
         reportPath.resize(0);
         int Rank;

         while(DFSCpath.size())
        {
               DFSC * DFSset = &DFSCpath[DFSCpath.size()-1];
               std::vector<int> parentSet;
               std::vector<int> CoVector;
               std::vector<int>  Central;
               std::vector<double> CentralValue;

               int order = DFSset->index;
               int old_rank = DFSset->Rank;
               int old_depth = DFSset->depth;

               parentSet.resize(DFSset->Projected.size());
               std::copy(DFSset->Projected.begin(),DFSset->Projected.end(),parentSet.begin());
               reportPath.resize(DFSset->depth+1);
               matPath.resize(Sample.n_rows,DFSset->depth+1);
               reportPath[DFSset->depth]  = DFSset->index;
               matPath.col(DFSset->depth) = Sample.col(DFSset->index);
               DFSCpath.pop();
               /*Filter Code*/


               if(rankFilter(old_rank,reportPath,parentSet))
                {

                   continue;
                }


               int MinRank = old_rank;
               bool Isminimum = true;
               int max_support = 0;
               int j;
               Central.resize(0);

             if(Retrict_Matrix.n_cols != 0)
            {

                for(j=0;j<Retrict_Matrix.n_cols;j++)
                {
                   bool Is_filter = false;
                   double minCorrelation = 1;

                      for(int k=0;k<reportPath.size();k++)
                      {
                          if( abs(Retrict_Matrix(reportPath[k],j))< abs(minCorrelation))
                          {
                                minCorrelation = Retrict_Matrix(reportPath[k],j);
                          }
                      }
                      if( abs(minCorrelation) >= range)
                      {
                            Central.resize(Central.size()+1);
                            CentralValue.resize(CentralValue.size()+1);
                            Central[Central.size()-1] = j;
                            CentralValue[CentralValue.size()-1] = minCorrelation;
                      }
                }

          }
               for(int i=0;i<order;i++)
               {


                  CoVector = Common_Set(ResultIDArray[i],parentSet);
                  if(rankSkip)
                  {

                        /*    int k;
                        for(k=0;k<reportPath.size();k++)
                        {
                            if(abs(Correlation_Matrix(reportPath[k],i))<range)
                            {
                                Rank = -1;
                                break;
                            }
                        } */

                     /*  if(k==reportPath.size())
                      {
                            Rank = 1;
                        }*/
                   if(CoVector.size() == 0 )
                   {
                         Rank = -1;

                    }else
                    {
                         Rank = 1;
                       }
                  }else{
                    Rank = RobustRank(i,matPath);
                  }

                  // CoVector = Common_Set(ResultIDArray[i],parentSet);

                  if(max_support<CoVector.size())
                  {
                     max_support=CoVector.size();
                  }


                  if(Rank>0)
                  {
                        DFSCpath.push(CoVector,i);
                        DFSCpath[DFSCpath.size()-1].depth = old_depth+1;
                        DFSCpath[DFSCpath.size()-1].Rank = Rank;
                        if(MinRank >= Rank)
                        {
                           Isminimum = false;
                        }
                  }
               }


                if(max_support<parentSet.size()||order==0||Isminimum)
                {
                    filterset.insert(parentSet);
                    rankInsert(old_rank,reportPath);
                    clusterEntry  Temp;
                    Temp.IDArray = reportPath;
                    Temp.Rank = old_rank;
                    clusterResult.push_back(Temp);
                    CentralStack.push_back(Central);
                    CentralVStack.push_back(CentralValue);
               }
        }
    }
}

void ComSeq::OutCluster(FILE * op)
{
    for(int i =0;i<clusterResult.size();i++)
    {
            fprintf(op,"ID %d ",i);
            fprintf(op," Rank: %d IDArray:",clusterResult[i].Rank);
             for(int j=0;j< clusterResult[i].IDArray.size();j++)
            {
                        fprintf(op," %d",clusterResult[i].IDArray[j]+1);
            }
            fprintf(op,"\n");
    }
}

void ComSeq::OutClusterWithName(FILE * op)
{

     std::vector <string> labels;

    for(int i =0;i<clusterResult.size();i++)
    {




             if( CentralStack[i].size() == 0 )
             {
                continue;
             }


            fprintf(op,"ID  %d ",i);
            fprintf(op," Rank:  %d  NameArray:",clusterResult[i].Rank);
            labels.resize(0);



             for(int j=0;j< clusterResult[i].IDArray.size();j++)
             {
                labels.push_back(SampleName[clusterResult[i].IDArray[j]]);
             }

            std::sort(labels.begin(),labels.end());

            for(int k=0;k<labels.size();k++)
            {
                fprintf(op," %s",labels[k].c_str());
            }

            fprintf(op," Central: ");
           // labels.resize(0);




           /* for(int j=0;j< CentralStack[i].size();j++)
            {
                labels.push_back(RestrictionName[CentralStack[i][j]]);
            }

            for(int k=0;k<labels.size();k++)
            {
                fprintf(op," %s",labels[k].c_str());
            }*/

            for(int j=0;j< CentralStack[i].size();j++)
            {
                fprintf(op," %s ",RestrictionName[CentralStack[i][j]].c_str());
                fprintf(op," %lf",CentralVStack[i][j]);
            }

            fprintf(op,"\n");

    }

}



void ComSeq::associte()
{

    for(int i=0;i<ResultIDArray.size();i++)
    {
        DFSRpath.push(ResultIDArray[i],i);
        DFSRpath[DFSRpath.size()-1].support = ResultIDArray[i].size();
        DFSRpath[DFSRpath.size()-1].depth = 0;
    }
{
    std::vector <int>reportPath;
    reportPath.resize(0);

    while(DFSRpath.size())
    {
         DFSR * DFSset = &DFSRpath[DFSRpath.size()-1];
         int order = DFSset->index;
         std::vector<int> CoVector;
         std::vector<int> parentSet;
         int max_support = 0;

         parentSet.resize(DFSset->Projected.size());
         int old_support = DFSset->Projected.size();
         int old_depth = DFSset->depth;
         std::copy(DFSset->Projected.begin(),DFSset->Projected.end(),parentSet.begin());
         reportPath.resize(DFSset->depth+1);
         reportPath[DFSset->depth]= DFSset->index;
         DFSRpath.pop();

         if(filterset.find(parentSet)!= filterset.end())
         {
                continue;
         }




         for(int i=0;i<order;i++)
        {
            CoVector = Common_Set(ResultIDArray[i],parentSet);
            if(max_support<CoVector.size())
            {
               max_support=CoVector.size();
            }
            if(CoVector.size()<=minsupport)
            {
               continue;
            }
            if(filterset.find(CoVector)!= filterset.end())
            {
               continue;
            }
            DFSRpath.push(CoVector,i);
            DFSRpath[DFSRpath.size()-1].depth = old_depth+1;
          }

            if((max_support<old_support)||order==0)
          {
            float Changed_Ratio = abs(max_support-old_support)*1.0/old_support*1.0;
            printf("%d\n",reportPath.size());
            filterset.insert(parentSet);
            Merged_Result Entry;
            Entry.Items = parentSet;
            Entry.PatternIDs = reportPath;
            Merged_Result_Array.push_back(Entry);
          }


    }

}
}




void ComSeq::output( FILE * op)
{
   int N =0;

   for(int i=0;i<Results.size();i++)
   {
      if(Results[i].size()<minlength)
      {
          continue;
      }else
       {
          N++;
       }
      for(int j =0;j<Results[i].size();j++)
      {
         fprintf(op,"%c",Results[i][j]);
      }
       fprintf(op,"%c",'\n');
   }
      fclose(op);
      printf("size :%d",N);
}
void ComSeq::filter()
{

  std::vector <filteKey> key;


  for(int i=0;i<TranArray.size();i++)
  {
         key.resize(0);
      for(int j=0;j<TranArray[i].size();j++)
      {
         key.resize(key.size()+1);
         key[key.size()-1].ID = TranArray[i][j].ID;
         key[key.size()-1].endsubID = TranArray[i][j].subID + TranArray[i][j].length-1;
      }
      if(FiltedResult.find(key)!=FiltedResult.end())
      {
          if(FiltedResult[key][0].length<TranArray[i][0].length)
          {
             FiltedResult[key] = TranArray[i];
          }

      }else{
          FiltedResult[key] = TranArray[i];
      }
  }
}
std::vector< std::list<std::string> >DataSet;
void featureread( FILE * ip)
{
    std::string tempword;
    std::list<std::string> tempsentence;
    char charbuffer[1024];

    while(!feof(ip))
   {
        fscanf(ip,"%[\b|\t]*",charbuffer);
       if(fscanf(ip,"%[^.]s",charbuffer)!=0)
       {
           tempword.assign(charbuffer);
           tempsentence.push_back(tempword);
       }
        fscanf(ip,"%[.]s",charbuffer);
        fscanf(ip,"%[\b|\t]*",charbuffer);
       if(fscanf(ip,"%[\n]*",charbuffer)!=0)
       {
           DataSet.push_back(tempsentence);
           tempsentence.resize(0);
       }

   }
}

void featurestat(FILE *filename)
{

          for(int i=0;i<DataSet.size();i++)
        {
              int count =0;
              bool matched = true;

                for(int j=0;j<Sample.size();j++)
                {
                        std::list<std::string>::iterator it;

                      const  char * str = Sample[j].c_str();


                    for (it = DataSet[i].begin(); it != DataSet[i].end(); it++)
                    {
                            std::string temp;
                            temp = (*it);
                            str = strstr(str,temp.c_str());
                        if(!str)
                        {
                                matched =false;
                                break;
                        }else
                        {
                            printf("matchesting:%s\n",temp.c_str());
                        }
                    }
            if (matched)
            {
				count++;
            }
        }
            fprintf(filename,"%d\n",count);
        }

}



void ComSeq::run_intern()
{

       unsigned int old_depth = 0;
       unsigned int max_support = 0;
       unsigned int old_support = 0;
       std::vector <PDFS> reportset;

       FiltedResult.clear();
       TranArray.clear();
       new_projected.clear();
       SubSeq.resize(1);

       for(int i=0;i<sample.size();i++){
           for(int j=0;j<sample[i].size();j++){
             PDFS cur;
             cur.ID = i;
             cur.subID = j;
             new_projected[sample[i][j]].push_back(cur);
           }
       }

       for(Projected_map::iterator cur =new_projected.begin();cur != new_projected.end(); ++cur){
            DFSpath.push(cur->first,cur->second,0);
            support(DFSpath[DFSpath.size()-1]);
            if(DFSpath[DFSpath.size()-1].support<minsupport)
            {
               DFSpath.pop();
            }

        }
       //inite


       while(DFSpath.size())
      {

        DFS * DFSset = &DFSpath[DFSpath.size()-1];
        new_projected.clear();

        for (unsigned int n = 0; n < DFSset->Projected.size(); ++n){
            PDFS *cur = &DFSset->Projected[n];
                if(cur->subID+DFSset->depth+1<sample[cur->ID].size())
                    new_projected[sample[cur->ID][cur->subID+DFSset->depth+1]].push_back(*cur);
                }

        SubSeq.resize(DFSset->depth+1);
        SubSeq[DFSset->depth] = DFSset->label;

        old_depth = DFSset->depth;
        old_support = DFSset->support;
        reportset = DFSset->Projected;
        DFSpath.pop();
        max_support = 0;

        for(Projected_map::iterator cur =new_projected.begin();cur != new_projected.end(); ++cur)
        {
            DFSpath.push(cur->first,cur->second,old_depth+1);
            support(DFSpath[DFSpath.size()-1]);
            if(DFSpath[DFSpath.size()-1].support<minsupport)
            {
                DFSpath.pop();
                continue;
            }//
            if(DFSpath[DFSpath.size()-1].support>max_support)
            {
                max_support = DFSpath[DFSpath.size()-1].support;
            }
        }
        if(max_support<old_support)
        {
             report(reportset,old_depth);
        }
         reportset.clear();
    }
     filter();
     merged_report();
}


