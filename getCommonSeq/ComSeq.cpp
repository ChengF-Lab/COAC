#include "ComSeq.h"
#include "stdio.h"
#include "math.h"
#include <fstream>
#include <boost/math/distributions/fisher_f.hpp>

using namespace std;
using namespace arma;
using namespace boost::math;






void ComSeq::ReadFilteFile(FILE *ip)
{
    char buffer[1024];

    while((!feof(ip)))
    {
         fscanf(ip,"%[\b|\t|\n]*",buffer);
         if( fscanf(ip,"%s",buffer)!=0)
          FilterSet.insert(buffer);
    }

}

void ComSeq::RestrictionRead( FILE *ip)
{
      std::vector<double>Temp;
      char seqName[1024];
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
      char seqName[1024];
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
  char seqName[1024];
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
	    fscanf(ip,"%s",buffer);
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
     char seqName[1024];
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



bool ComSeq::rankFilter(int Rank, std::vector <int> ReportPath,std::vector<int> parentSet)
{




       if((filterset.find(parentSet) == filterset.end()))
       {
           if(FilterSetID.size()!=0)
           {
              for(int i=0;i<parentSet.size();i++)
             {
                if( FilterSetID.find(parentSet[i])!= FilterSetID.end() )
                {
                   return false;
                }
             }
                  return true;
           }

            return false;
       }else
         {
            return true;
         }



    /* if(rankMap.find(Rank)==rankMap.end())
   {
         return false;
   }else
   {
          /*
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
          */
    //      return true;
   //}
            //    return false;
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

void ComSeq::filter()
{
     mat Temp_Sample;
     std::vector<std::string> TempSampleName;
     Temp_Sample.resize(Sample.n_rows,Sample.n_cols);

     int colSize=0;

     for(int j=0;j<Sample.n_cols;j++)
     {
          int NozeroCount = 0;
          for(int i=0;i<Sample.n_rows;i++)
          {
               if(Sample(i,j)>0)
               {
                   NozeroCount++;
               }
          }
          if(NozeroCount>Sample.n_rows*NoZeroRatio)
          {
               for(int i=0;i<Sample.n_rows;i++)
               {
                   Temp_Sample(i,colSize)=Sample(i,j);
               }
               colSize++;
               TempSampleName.resize(colSize);
               TempSampleName[TempSampleName.size()-1]=SampleName[j];
          }
     }

     Temp_Sample.resize(Sample.n_rows,colSize);
     Sample.resize(0,0);
     Sample=Temp_Sample;
     Temp_Sample.resize(0,0);
     SampleName.resize(0);
     SampleName = TempSampleName;

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
  //ofstream  U_d;
 // ofstream  s_d;
 // ofstream  V_d;
 // ofstream  M_d;

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
    mat U;
    vec s;
    mat V;
    std::vector<double> ComponentDistribution;





    svd_econ(U,s,V,Sample);

    for(int i=0;i<s.n_elem;i++)
    {
       AC = AC+s[i]*s[i];
    }

    for(int i=0;i<s.n_elem;i++)
    {
       PC = PC+s[i]*s[i];
       Rank ++;
       if(PC>((1.0-tolerance)*AC))
       {
           break;
       }
    }
    printf("Rank: %d\n",Rank);

  //  U.resize(U.n_rows,Rank);
    //cor= Sample.t()*U;
    cor.resize(Sample.n_cols,Rank);


    for(int j=0;j< cor.n_cols;j++)
    {
           for(int i=0;i< cor.n_rows;i++)
           {
              cor(i,j) = V(i,j)*s[j];
           }
    }

    for(int j=0;j< cor.n_cols;j++)
    {
      for(int i=0;i< cor.n_rows;i++)
      {
             if(Correlation_Matrix[1,i]==0)
             {
                cor(i,j) =0;
             }
             else
             {
                cor(i,j) = cor(i,j)/sqrt(Correlation_Matrix[1,i]);
             }
      }
    }


    ResultIDArray.resize(2*Rank);
    ComponentDistribution.resize(cor.n_rows);

    for(int j=0;j<cor.n_cols;j++)
    {
        int NegtiveSize =0;
        double  NegativeGroupAverage = 0;
        double  PositiveGroupAverage = 0;
        double  NegativeGroupSquare  = 0;
        double  PositiveGroupSquare  = 0;
        double  allAverage  = 0;
        double  allSquare   = 0;
        bool UNFILTER = false;





        for(int i=0;i<cor.n_rows;i++)
        {
            ComponentDistribution[i] = cor(i,j);
            if(ComponentDistribution[i]>0)
            {
                PositiveGroupAverage = PositiveGroupAverage + ComponentDistribution[i];
                PositiveGroupSquare  = PositiveGroupSquare + ComponentDistribution[i]*ComponentDistribution[i];
            }else
            {
                NegativeGroupAverage = NegativeGroupAverage + ComponentDistribution[i];
                NegativeGroupSquare  = NegativeGroupSquare + ComponentDistribution[i]*ComponentDistribution[i];
            }
        }

       allAverage = PositiveGroupAverage + NegativeGroupAverage;
       allSquare  = PositiveGroupSquare + NegativeGroupSquare;
       double sd = allSquare/((cor.n_rows-1)*1.0);
       double mean = allAverage/((cor.n_rows)*1.0);
       sd = sd - mean*mean;
       sd = sqrt(sd);


      std::vector<size_t>Index = ordered(ComponentDistribution,true);


        for(int i=0;i< ComponentDistribution.size();i++)
       {
           if(ComponentDistribution[Index[i]]>0)
           {
               break;
           }else
           {
               NegtiveSize++;
           }
       }

         if( NegtiveSize<ComponentDistribution.size() )
         {

             int PositiveSize = ComponentDistribution.size()-NegtiveSize;

                mat data(1,PositiveSize);
                for(int k=0;k<PositiveSize;k++)
                {
                    data(0,k) = ComponentDistribution[Index[ComponentDistribution.size()-k-1]];
                }
                if(data.n_cols>2)
                {
                     double  Group1Average = 0;
                     double  Group1Square  = 0;
                     double  Group2Average = 0;
                     double  Group2Square  = 0;
                     double  MaxGroup1Average = 0;
                     double  MaxGroup1Square = 0;
                     double  MinSqrtAverage = -1;
                     double  SqrtAverage;
                     double  Max_F;
                     int MinvarID = 0;
                     int GroupSize = data.n_cols;
                     double F;
                     boost::math::fisher_f_distribution<double>s(1,cor.n_rows-1);
                     double  scalar_likelihood;


                     //std::fisher_f_distribution<double> distribution(1.0,data.n_cols-2);



                for(int step = 0;step<16;step++)
               {


                    MinSqrtAverage = (0xFFFFFFFF)*1.0;
                    Group1Average  = 0;
                    Group1Square   = 0;



                    for(int k=0;k<GroupSize;k++)
                   {

                     Group1Average = Group1Average + data(0,k);
                     Group1Square  = Group1Square + data(0,k)*data(0,k);
                     Group2Average = PositiveGroupAverage - Group1Average;
                     Group2Square  = PositiveGroupSquare - Group1Square;


                     SqrtAverage =  (Group1Square-Group1Average*Group1Average/((k+1)*1.0))*((k+1)*1.0) + (Group2Square-(Group2Average*Group2Average/((GroupSize-k)*1.0)))*(GroupSize-k);

                       if(SqrtAverage < MinSqrtAverage)
                      {
                         MinSqrtAverage = SqrtAverage;
                         MaxGroup1Square = Group1Square;
                         MaxGroup1Average = Group1Average;
                         MinvarID = k;
                      }
                    }

                    {

                       double x = data(0,MinvarID);

                       scalar_likelihood = 1-cdf(s,((x-mean)/sd)*((x-mean)/sd));
                       scalar_likelihood = log(scalar_likelihood);
                    }


                     PositiveGroupSquare = MaxGroup1Square;
                     PositiveGroupAverage = MaxGroup1Average;
                     GroupSize = MinvarID;
                      if( scalar_likelihood < Log_Pvalue )
                      {
                        std::vector < double > CorSum;
                        double MinCorAverageRatio = (0xFFFFFFFF)*1.0;
                        int MinID=0;
                        double CorSumSquare=0;
                        double SelectGroup=0;
                        double SelectGroupAverage=0;

                          CorSum.resize(cor.n_cols);

                          for(int j=0;j<cor.n_cols;j++)
                          {

                                CorSum[j]=0;

                                for(int i=0;i<MinvarID+1;i++)
                                {
                                   CorSum[j]=CorSum[j]+cor(Index[i],j);
                                };
                                CorSumSquare =CorSumSquare + CorSum[j]*CorSum[j];
                          }

                           std::vector < double >SubCorSum;

                           SubCorSum.resize(cor.n_cols);


                           for(int i=0;i<cor.n_cols;i++)
                           {
                              SubCorSum[i]=0;
                           }

                            for(int i=0;i<MinvarID;i++)
                           {

                              for(int j=0;j<cor.n_cols;j++)
                              {
                                 SubCorSum[j]=SubCorSum[j]+cor(Index[i],j);
                              };

                             double GroupBetween = CorSumSquare;
                              SelectGroup = 0;

                             for(int j=0;j<cor.n_cols;j++)
                             {
                               GroupBetween = GroupBetween-((CorSum[j]-SubCorSum[j])*(CorSum[j]-SubCorSum[j])+SubCorSum[j]*SubCorSum[j]);
                               SelectGroup = SelectGroup + SubCorSum[j]*SubCorSum[j];
                             }

                             if(i>minsupport&i>1)
                             {

                                 double CorAverageRatio = GroupBetween/((MinvarID+1)*(MinvarID+1)-(i+1)*(i+1)-(MinvarID-i)*(MinvarID-i));


                                 if(MinCorAverageRatio>CorAverageRatio)
                                 {
                                    MinCorAverageRatio = CorAverageRatio;
                                    MinID =i;
                                    SelectGroupAverage = abs((SelectGroup-i-1)/((i+1)*(i)*1.0));
                                 }
                             }
                            }
                            if( abs(MinCorAverageRatio/SelectGroupAverage)<CutThresHold )
                            {
                               MinvarID = MinID;
                            }

                         break;
                      }


                 }

                 if(scalar_likelihood < Log_Pvalue)
                 {
                    for(int k=0;k<MinvarID+1;k++)
                   {
                     ResultIDArray[j].resize(ResultIDArray[j].size()+1);
                     ResultIDArray[j][ResultIDArray[j].size()-1] = Index[ComponentDistribution.size()-k-1];

                     if(FilterSet.size()!=0)
                     {
                       if( FilterSet.find( SampleName[Index[k]].c_str()) != FilterSet.end() )
                       {
                         UNFILTER = true;
                       }
                     }


                   }

                 }
                }
         }
        if(NegtiveSize!=0)
        {



             mat data(1,NegtiveSize);

             for(int k=0;k<NegtiveSize;k++)
             {
                data(0,k) = ComponentDistribution[Index[k]];
             }

              if(data.n_cols>2)
             {

                     double  Group1Average = 0;
                     double  Group1Square  = 0;
                     double  Group2Average = 0;
                     double  Group2Square  = 0;
                     double  MaxGroup1Average = 0;
                     double  MaxGroup1Square = 0;
                     double  MinSqrtAverage = -1;
                     double  SqrtAverage;
                     double  scalar_likelihood;
                     double  Max_F;
                     int MinvarID = 0;
                     int GroupSize = data.n_cols;
                     double F;
                     boost::math::fisher_f_distribution <double> s(1,cor.n_rows-1);


               //  std::fisher_f_distribution<double> distribution(1.0,data.n_cols-2);
                for(int step = 0;step<16;step++)
               {


                    MinSqrtAverage = (0xFFFFFFFF)*1.0;
                    Group1Average  = 0;
                    Group1Square   = 0;


                 for(int k=0;k<GroupSize;k++)
                {



                     Group1Average = Group1Average + data(0,k);
                     Group1Square  = Group1Square + data(0,k)*data(0,k);
                     Group2Average = NegativeGroupAverage - Group1Average;
                     Group2Square  = NegativeGroupSquare - Group1Square;
                     SqrtAverage =  (Group1Square-Group1Average*Group1Average/((k+1)*1.0))*((k+1)*1.0) + (Group2Square-Group2Average*Group2Average/((GroupSize-k)*1.0));

                      if(SqrtAverage < MinSqrtAverage)
                      {
                         MinSqrtAverage = SqrtAverage;
                         MaxGroup1Square = Group1Square;
                         MaxGroup1Average = Group1Average;
                         MinvarID = k;
                      }
                    // scalar_likelihood = log(distribution(F));
                 };



                    {

                       double x = data(MinvarID);

                       scalar_likelihood = 1- cdf(s,((x-mean)/sd)*((x-mean)/sd));
                       scalar_likelihood = log(scalar_likelihood);

                    }
                     NegativeGroupSquare = MaxGroup1Square;
                     NegativeGroupAverage = MaxGroup1Average;
                     GroupSize = MinvarID;
                     if( scalar_likelihood < Log_Pvalue )
                      {
                        std::vector < double > CorSum;
                        double MinCorAverageRatio = (0xFFFFFFFF)*1.0;
                        int MinID=0;
                        double CorSumSquare=0;
                        double SelectGroup=0;
                        double SelectGroupAverage=0;

                          CorSum.resize(cor.n_cols);

                          for(int j=0;j<cor.n_cols;j++)
                          {

                                CorSum[j]=0;

                                for(int i=0;i<MinvarID+1;i++)
                                {
                                   CorSum[j]=CorSum[j]+cor(Index[i],j);
                                };
                                CorSumSquare =CorSumSquare + CorSum[j]*CorSum[j];
                          }

                           std::vector < double >SubCorSum;

                           SubCorSum.resize(cor.n_cols);


                           for(int i=0;i<cor.n_cols;i++)
                           {
                              SubCorSum[i]=0;
                           }

                            for(int i=0;i<MinvarID;i++)
                           {

                              for(int j=0;j<cor.n_cols;j++)
                              {
                                 SubCorSum[j]=SubCorSum[j]+cor(Index[i],j);
                              };

                             double GroupBetween = CorSumSquare;

                             SelectGroup = 0;
                             for(int j=0;j<cor.n_cols;j++)
                             {
                               GroupBetween = GroupBetween-((CorSum[j]-SubCorSum[j])*(CorSum[j]-SubCorSum[j])+SubCorSum[j]*SubCorSum[j]);
                               SelectGroup = SelectGroup + SubCorSum[j]*SubCorSum[j];
                             }

                             if(i>minsupport&i>1)
                             {
                                 double CorAverageRatio = GroupBetween/((MinvarID+1)*(MinvarID+1)-(i+1)*(i+1)-(MinvarID-i)*(MinvarID-i));


                                 if(MinCorAverageRatio>CorAverageRatio)
                                 {
                                    SelectGroupAverage = abs((SelectGroup-i-1)/((i+1)*(i)*1.0));
                                    MinCorAverageRatio = CorAverageRatio;
                                    MinID =i;
                                 }
                             }
                         }
                      if( abs(MinCorAverageRatio/SelectGroupAverage)<CutThresHold )
                       {
                           MinvarID = MinID;
                       }

                         break;
                      }
               }

                if(scalar_likelihood < Log_Pvalue) {

                  for(int k=0;k<MinvarID+1;k++)
                  {


                     ResultIDArray[j+Rank].resize(ResultIDArray[j+Rank].size()+1);
                     ResultIDArray[j+Rank][ResultIDArray[j+Rank].size()-1] = Index[k];
                    if(FilterSet.size()!=0)
                    {
                       if( FilterSet.find( SampleName[Index[k]].c_str()) != FilterSet.end() )
                      {
                        UNFILTER = true;
                      }
                    }
                  };
                }
             }
         }
       if(FilterSet.size()!=0)
       {
         if(UNFILTER==false)
         {
            ResultIDArray[j+Rank].resize(0);
            ResultIDArray[j].resize(0);
         }
       }
    }


     FilterSetID.clear();

     for(int i=0;i<SampleName.size();i++)
     {
        if( FilterSet.find( SampleName[i].c_str()) != FilterSet.end() )
        {
            FilterSetID.insert(i);
        }
     }

     for(int i=0;i< ResultIDArray.size();i++)
     {
              /*if(cor(i,j)>= (1.0/sqrt(s.n_elem))*DimRatio)
                {

                     ResultIDArray[i].resize(ResultIDArray[i].size()+1);
                     ResultIDArray[i][ResultIDArray[i].size()-1] = j;
                }
                if(cor(i,j) <= -(1.0/sqrt(s.n_elem))*DimRatio)
                {
                     ResultIDArray[i].resize(ResultIDArray[i].size()+1);
                     ResultIDArray[i][ResultIDArray[i].size()-1] = j+Rank;
                }*/
                std:sort(ResultIDArray[i].begin(),ResultIDArray[i].end());
     }





}

void ComSeq::AttributeCluster()
{
      DFSCpath.resize(0);


     // Correlation_Matrix =  Sample.t()*Sample;
    //  assert(Correlation_Matrix.n_rows==Correlation_Matrix.n_cols);
      //assert(Correlation_Matrix.n_rows==Sample.n_cols);

    /*  for(int i=0;i<Correlation_Matrix.n_rows;i++)
      {
         for(int j=0;j<Correlation_Matrix.n_cols;j++)
         {
               if(i==j)
                 continue;
                Correlation_Matrix(i,j) = Correlation_Matrix(i,j)/sqrt(Correlation_Matrix(i,i)*Correlation_Matrix(j,j));
          }
      }*/

        Correlation_Matrix.set_size(1,Sample.n_cols);

        for( int i=0;i<Sample.n_cols;i++)
        {
             mat  Rn = Sample.col(i).t()*Sample.col(i);
             Correlation_Matrix[1,i] = Rn(0,0);
        }

        if(Restriction.n_rows!=0&&Restriction.n_cols!=0)
        {
           Retrict_Matrix =  Sample.t()*Restriction;

            for(int j=0;j< Retrict_Matrix.n_cols;j++)
           {
              mat  Rn = Restriction.col(j).t()*Restriction.col(j);


              for(int i=0;i< Retrict_Matrix.n_rows;i++)
              {
                if(Correlation_Matrix[1,i]==0||Rn(0,0)==0)
                {
                     Retrict_Matrix(i,j) =0;
                }else{
                     Retrict_Matrix(i,j) =Retrict_Matrix(i,j)/sqrt(Correlation_Matrix[1,i]*(Rn(0,0)));
                  }
              }
           }
       }


   //      Retrict_Matrix.print("EE:");
         PreAttributeCluster();

     /*    for(int i=0;i<Correlation_Matrix.n_rows;i++)
        {

            Correlation_Matrix(i,i) = 1.0;

        }*/


      for(int i=0;i<ResultIDArray.size();i++)
     {
            DFSCpath.push(ResultIDArray[i],i);
            DFSCpath[DFSCpath.size()-1].Rank  = 1;
            DFSCpath[DFSCpath.size()-1].depth = 0;
     }

     {
         std::vector <int>reportPath;
         reportPath.resize(0);
         int Rank;

         while(DFSCpath.size())
        {
               DFSC * DFSset = &DFSCpath[DFSCpath.size()-1];
               std::vector<int> parentSet;
               std::vector<int> CoVector;


               int order = DFSset->index;
               int old_rank = DFSset->Rank;
               int old_depth = DFSset->depth;

               parentSet.resize(DFSset->Projected.size());
               std::copy(DFSset->Projected.begin(),DFSset->Projected.end(),parentSet.begin());
               reportPath.resize(DFSset->depth+1);
               reportPath[DFSset->depth]  = DFSset->index;
               DFSCpath.pop();
               /*Filter Code*/


               if(rankFilter(old_rank,reportPath,parentSet))
                {

                   continue;
                }

                if(reportPath.size()> MaxDepth)
                {
                    continue;
                }
               int MinRank = old_rank;
               bool Isminimum = true;
               int max_support = 0;
               int j;








         /*    if(Retrict_Matrix.n_cols != 0)
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

          }*/
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
                   if(CoVector.size() < minsupport)
                   {
                         Rank = -1;

                    }else
                    {
                         Rank = 1;
                       }
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




                if((max_support<parentSet.size()||order==0||Isminimum)&(!(rankFilter(old_rank,reportPath,parentSet))))
                {
                    filterset.insert(parentSet);
                    // rankInsert(old_rank,reportPath);

                  if(parentSet.size() >= minsupport)
                 {
                     clusterEntry  Temp;

                     Temp.IDArray = reportPath;
                     Temp.Rank = old_rank;
                     Temp.PatternIDs = parentSet;
                     clusterResult.push_back(Temp);
                    //CentralStack.push_back(Central);
                    //CentralVStack.push_back(CentralValue);
                  }
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

    fprintf(op,"Total size:%d\n", (int)clusterResult.size());


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

//            std::sort(labels.begin(),labels.end());

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
            std::vector<size_t>CentralIndex = ordered(CentralVStack[i],false);

            for(int j=0;j<CentralIndex.size();j++)
            {
                fprintf(op," %s ",RestrictionName[CentralStack[i][(int)CentralIndex[j]]].c_str());
                fprintf(op," %lf",CentralVStack[i][(int)CentralIndex[j]]);
            }

            fprintf(op,"\n");



             for(int i1=0;i1< clusterResult[i].IDArray.size();i1++)
            {
                    int SampleID = clusterResult[i].IDArray[i1];
                    double norm = 0.0;

                    for(int j=0;j< clusterResult[i].PatternIDs.size();j++)
                    {

                        int ComponentID = clusterResult[i].PatternIDs[j];
                        if(ComponentID>=Rank)
                        {
                            ComponentID = ComponentID - Rank;
                        }

                        norm = norm + cor(SampleID,ComponentID)*cor(SampleID,ComponentID)*(Correlation_Matrix[1,SampleID]);

                    }
                   norm = sqrt(norm);

                   fprintf(op," %lf",norm);
            }
                   fprintf(op,"\n");
      }

    }

void ComSeq::OutClusterWithNameCombact(FILE * op)
{

    std::vector <string> labels;
    std::vector<int>  Central;
    std::vector<double> CentralValue;
    std::vector<double> PatternNormalized;
    std::vector<double> PatternNormalizedTemp;
    std::vector<double> ComponentNormalizedTemp;

    std::vector<std::vector<double> > Orthogonal;

    Orthogonal.resize(0);

    PatternNormalized.resize(cor.n_rows);
    ComponentNormalizedTemp.resize(cor.n_cols);
    PatternNormalizedTemp.resize(cor.n_rows);
    mat CorTemp;

    /*for(int i=0;i<PatternNormalized.size();i++)
    {
           PatternNormalized[i] = 0;
    }


    for(int i=0;i<PatternNormalized.size();i++)
    {
        for(int j=0;j<cor.n_cols;j++)
        {
           PatternNormalized[i] = PatternNormalized[i] + cor(i,j)*cor(i,j)*Correlation_Matrix[1,i];
        }
    }


        for(int i=0;i<PatternNormalized.size();i++)
       {
           PatternNormalized[i] = sqrt(PatternNormalized[i]);
       }

*/


    fprintf(op," Total Size: %d\n",(int)clusterResult.size());



    Simple_Nozero.resize(Sample.n_cols);


      for(int j=0;j<Sample.n_cols;j++)
     {
         for(int i=0;i<Sample.n_rows;i++)
         {
             if(Sample(i,j)!=0)
              Simple_Nozero[j].push_back(i);
         }
     }


    std::vector <int> OriginalCommonSet;
    std::vector <int> CoVector;
    int filtered_Total_size =0;





    for(int i =0;i<clusterResult.size();i++)
    {

            if( clusterResult[i].PatternIDs.size()==0)
            {
                continue;
            }


             OriginalCommonSet=Simple_Nozero[clusterResult[i].PatternIDs[0]];






              for(int j=1;j< clusterResult[i].PatternIDs.size();j++)
              {
                  CoVector = Common_Set(Simple_Nozero[clusterResult[i].PatternIDs[j]],OriginalCommonSet);
                  OriginalCommonSet = CoVector;
              }

              if(OriginalCommonSet.size()<(Sample.n_rows*NoZeroRatio*0.5))
              {
                    continue;
              }

           fprintf(op,"%d\n",(int)clusterResult[i].PatternIDs.size());
           filtered_Total_size++;

             for(int j=0;j< clusterResult[i].PatternIDs.size();j++)
             {
                 fprintf(op," %s",SampleName[clusterResult[i].PatternIDs[j]].c_str());
                 PatternNormalized[clusterResult[i].PatternIDs[j]] = 0;
                 for(int j1=0;j1< clusterResult[i].IDArray.size();j1++)
                  {
                    int ComponentID = clusterResult[i].IDArray[j1];
                     if(ComponentID>=Rank)
                     {
                        ComponentID = ComponentID - Rank;
                     }
                    PatternNormalized[clusterResult[i].PatternIDs[j]] = PatternNormalized[clusterResult[i].PatternIDs[j]] + cor(clusterResult[i].PatternIDs[j],ComponentID)*cor(clusterResult[i].PatternIDs[j],ComponentID)*Correlation_Matrix[1,clusterResult[i].PatternIDs[j]];
                  }

                   PatternNormalized[clusterResult[i].PatternIDs[j]] = sqrt(PatternNormalized[clusterResult[i].PatternIDs[j]]);

             }

                fprintf(op,"\n");

                fprintf(op,"%d  \n",(int)clusterResult[i].IDArray.size());

             double CorrSum =0;
             double CorrSumAdd =0;
             double InterCor=0;

             for(int j=0;j< clusterResult[i].IDArray.size();j++)
             {
                    int ComponentID = clusterResult[i].IDArray[j];

                    if(ComponentID>=Rank)
                     {
                        ComponentID = ComponentID - Rank;
                     }
                        fprintf(op," %d",ComponentID);


                   double Corr =0;

                   for(int i1=0;i1< clusterResult[i].PatternIDs.size();i1++)
                   {
                       int SampleID = clusterResult[i].PatternIDs[i1];

                       fprintf(op," %lf",cor(SampleID,ComponentID)*sqrt(Correlation_Matrix[1,SampleID]));

                       Corr = cor(SampleID,ComponentID)*sqrt(Correlation_Matrix[1,SampleID])/PatternNormalized[SampleID] + Corr;
                       InterCor = InterCor + cor(SampleID,ComponentID)*cor(SampleID,ComponentID);
                    }



                     CorrSum = Corr*Corr+CorrSum;

                     fprintf(op,"\n");
               }



               if( Orthogonal.size()<clusterResult[i].IDArray.size())
                {
                    Orthogonal.resize(clusterResult[i].IDArray.size());
                }
               for(int j=0;j<clusterResult[i].IDArray.size();j++)
               {
                 if( Orthogonal[j].size()< clusterResult[i].PatternIDs.size())
                 {
                     Orthogonal[j].resize(clusterResult[i].PatternIDs.size());
                 }
               }


               double normalize =0.0;


               CorTemp.resize(clusterResult[i].IDArray.size(),clusterResult[i].PatternIDs.size());


               for(int j=0;j< clusterResult[i].IDArray.size();j++)
                {
                    int ComponentID = clusterResult[i].IDArray[j];

                    if(ComponentID>=Rank)
                    {
                        ComponentID = ComponentID - Rank;
                    }
                    double OrthogonalComponent =0;
                    for(int i1=0;i1< clusterResult[i].PatternIDs.size();i1++)
                    {
                        int SampleID = clusterResult[i].PatternIDs[i1];
                        CorTemp(j,i1) = cor(SampleID,ComponentID)*sqrt(Correlation_Matrix[1,SampleID]);
                    }
               }

               svd(U,s,V,CorTemp);

               double AC = 0;
               double PC = 0;

               for(int i=0;i<s.n_elem;i++)
               {
                 AC = AC+s[i]*s[i];
               }
               int Rank_ = 0;
               for(int i=0;i<s.n_elem;i++)
               {
                  PC = PC+s[i]*s[i];
                  Rank_ ++;

                  if(PC>ScannedTolerance*ScannedTolerance*AC)
                   {
                      break;
                   }
                }

               for(int j=0;j< clusterResult[i].IDArray.size();j++)
               {

                   for(int i1=0;i1< clusterResult[i].PatternIDs.size();i1++)
                   {
                       if(j<Rank_ )
                       Orthogonal[j][i1]=V(i1,j);
                       else
                       Orthogonal[j][i1]=0.0;
                   }
              }


               for(int j=0;j<cor.n_cols;j++)
              {
                   double Corr =0;
                   CorrSumAdd=0;

                    if(std::find(clusterResult[i].IDArray.begin(),clusterResult[i].IDArray.end(),j)!=clusterResult[i].IDArray.end()||std::find(clusterResult[i].IDArray.begin(),clusterResult[i].IDArray.end(),(j+Rank))!=clusterResult[i].IDArray.end())
                   {
                     continue;
                   }



                  for(int i1=0;i1<clusterResult[i].PatternIDs.size();i1++)
                  {
                      PatternNormalizedTemp[clusterResult[i].PatternIDs[i1]] = PatternNormalized[clusterResult[i].PatternIDs[i1]]*PatternNormalized[clusterResult[i].PatternIDs[i1]] + cor(clusterResult[i].PatternIDs[i1],j)*cor(clusterResult[i].PatternIDs[i1],j)*Correlation_Matrix[1,clusterResult[i].PatternIDs[i1]];
                      PatternNormalizedTemp[clusterResult[i].PatternIDs[i1]] = sqrt(PatternNormalizedTemp[clusterResult[i].PatternIDs[i1]]);
                  }


                   for(int j=0;j< clusterResult[i].IDArray.size();j++)
                   {
                        int ComponentID = clusterResult[i].IDArray[j];
                        Corr =0;


                        if(ComponentID>=Rank)
                        {
                            ComponentID = ComponentID - Rank;
                        }

                        for(int i1=0;i1< clusterResult[i].PatternIDs.size();i1++)
                        {
                            int SampleID = clusterResult[i].PatternIDs[i1];
                            Corr = cor(SampleID,ComponentID)*sqrt(Correlation_Matrix[1,SampleID])/PatternNormalizedTemp[SampleID] + Corr;
                        }
                        CorrSumAdd = Corr*Corr+CorrSumAdd;
                 }


                 Corr =0;

                 for(int i1=0;i1< clusterResult[i].PatternIDs.size();i1++)
                 {
                   int SampleID = clusterResult[i].PatternIDs[i1];
                   Corr = cor(SampleID,j)*sqrt(Correlation_Matrix[1,SampleID])/PatternNormalizedTemp[SampleID] + Corr;
                 }

                CorrSumAdd = Corr*Corr+CorrSumAdd;




                double CommmonComponent = 0;

                for(int j1=0;j1<clusterResult[i].IDArray.size();j1++)
                {
                    double OrthogonalComponent=0;

                    for(int i1=0;i1<clusterResult[i].PatternIDs.size();i1++)
                    {
                        int SampleID = clusterResult[i].PatternIDs[i1];

                       OrthogonalComponent =OrthogonalComponent+cor(SampleID,j)*sqrt(Correlation_Matrix[1,SampleID])*Orthogonal[j1][i1];

                    }
                    CommmonComponent = CommmonComponent + OrthogonalComponent*OrthogonalComponent;
                }

                CommmonComponent = sqrt(CommmonComponent);

                double AllComponent=0;

                for(int i1=0;i1<clusterResult[i].PatternIDs.size();i1++)
                {
                    int SampleID = clusterResult[i].PatternIDs[i1];
                    AllComponent = AllComponent+ cor(SampleID,j)*Correlation_Matrix[1,SampleID]*cor(SampleID,j);
                }

                AllComponent = sqrt(AllComponent);





                if((CorrSum<=(CorrSumAdd+0.00001)))
                {
                   fprintf(op," %d",j);

                   double ComponentRatio = (CommmonComponent)/(AllComponent+0.000000001);
                   double InterCor_=0;

                   for(int i1=0;i1< clusterResult[i].PatternIDs.size();i1++)
                   {
                       int SampleID = clusterResult[i].PatternIDs[i1];
                       InterCor_ = InterCor_ + cor(SampleID,j)*cor(SampleID,j);
                       fprintf(op," %lf",cor(SampleID,j)*sqrt(Correlation_Matrix[1,SampleID]));
                   }


                   fprintf(op," %lf",ComponentRatio);
                   InterCor = InterCor + InterCor_*ComponentRatio;

                   fprintf(op,"\n");

                }

             }



             if(clusterResult[i].PatternIDs.size()>0)
             {
                fprintf(op,"%s","CorrelationAverage: ");
                fprintf(op,"%d ",(int)clusterResult[i].PatternIDs.size());
		if( clusterResult[i].PatternIDs.size() > 1.0)
                    fprintf(op,"%lf ",(CorrSum-clusterResult[i].PatternIDs.size())/((clusterResult[i].PatternIDs.size()-1)*(clusterResult[i].PatternIDs.size())));
                else
		    fprintf(op,"%lf ",1.0);
			 
               /* double InterCor1=0;
                for(int j=0;j<cor.n_cols;j++)
                {
                    for(int i1=0;i1< clusterResult[i].PatternIDs.size();i1++)
                   {
                       int SampleID = clusterResult[i].PatternIDs[i1];
                       InterCor1 = InterCor1 + cor(SampleID,j)*cor(SampleID,j);
                   }
                } */
                fprintf(op,"%lf\n",( InterCor/(clusterResult[i].PatternIDs.size()*1.0)));
              }


            }
               fprintf(op,"Filtered Total Size:%d",filtered_Total_size);
               fprintf(op,"\n");
      }













