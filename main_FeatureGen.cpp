#include  <iostream>
#include  "stdio.h"
#include  <cmath>
#include  "getCommonSeq/ComSeq.h"
using namespace std;

int main( int argc , char *argv[] )
{

     ClusterToFeature Test;
     
     if(argv[1]==NULL)
     {
          printf("FeatureGen ToFeature\n");
          printf("FeatureGen Decipher\n");
          printf("FeatureGen Statistics\n");
          printf("FeatureGen CorStatistics\n");
          return 0;
     }
     
   

     if( strcmp(argv[1],"ToFeature")==0)
     {


        FILE * ip =  fopen(argv[2],"r");
        FILE * op =  fopen(argv[3],"w+");
          
        if(ip==NULL||op==NULL||argc<5)
        {
            printf("FeatureGen ToFeature   InputFile OutputFIle AverageCorrelation AverageComponentRatio\n");
            return 0;
        }
        Test.tolerance = 0.01;
        Test.CorThreshold = atof(argv[4]);
        Test.CompRThreshold = atof(argv[5]);
        Test.ClusterRead(ip);
        Test.OutPutFeature(op);
     }

     if( strcmp(argv[1],"Decipher")==0)
     {


        FILE * ip =  fopen(argv[2],"r");
        FILE * ip1 = fopen(argv[3],"r");
        FILE * op =  fopen(argv[4],"w+");
          
         if(ip==NULL||op==NULL||argc<4)
         {
            printf("FeatureGen Decipher FeatureFile   InputFile OutputFile\n");
            return 0;
         }
          
        Test.ReadFeature(ip);
        Test.SampleRead(ip1);
        Test.SampleDecipher(op);

     }



     if( strcmp(argv[1],"Statistics")==0)
     {
        char buffer[1024];
        std::vector < std::vector < unsigned int > > Index;
        std::vector < std::vector < double > > Value;


        FILE * ip =  fopen(argv[2],"r");
        FILE * ip1 = fopen(argv[3],"r");
        FILE * op = fopen(argv[4],"w+");

         if(ip==NULL||op==NULL||argc<4)
         {
            printf("FeatureGen Statistics InputFile InputFileGroupTag OutputFile\n");
            return 0;
              
         }

        std::vector < double > CorSum;
        std::vector < double > SquareSum;

        Test.ReadSpareM(Index,Value,ip);

        std::vector < unsigned int > Flag;
        Flag.resize(Index.size());

        for(int i=0;i<Index.size();i++ )
        {
            fscanf(ip1,"\t");
            fscanf(ip1,"%d",&Flag[i]);
            fscanf(ip1,"\n");
        }

        Test.StatisticSpareMGroup(Index,Value,Flag,CorSum,SquareSum);

        for(int i=0;i<CorSum.size();i++)
        {
            fprintf(op,"%lf", CorSum[i]);
            fprintf(op,"\t");
            fprintf(op,"%lf", SquareSum[i]);
            fprintf(op,"\n");
        }
     }

      if( strcmp(argv[1],"CorStatistics")==0)
     {
        char buffer[1024];
        std::vector < std::vector < unsigned int > > Index;
        std::vector < std::vector < double > > Value;


        FILE * ip =  fopen(argv[2],"r");
        FILE * ip1 = fopen(argv[3],"r");
        FILE * op = fopen(argv[4],"w+");

        
         if(ip==NULL||op==NULL||argc<4)
         {
            printf("FeatureGen CorStatistics InputFile InputFileNumericalTag OutputFile\n");
            return 0;
         }  
         
           

        std::vector < double > CorSum;

        Test.ReadSpareM(Index,Value,ip);

        std::vector < double > Flag;
        Flag.resize(Index.size());

        for(int i=0;i<Index.size();i++ )
        {
            fscanf(ip1,"\t");
            fscanf(ip1,"%lf",&Flag[i]);
            fscanf(ip1,"\n");
        }

        Test.StatisticSpareCor(Index,Value,Flag,CorSum);

        for(int i=0;i<CorSum.size();i++)
        {
            fprintf(op,"%lf",CorSum[i]);
            fprintf(op,"\n");
        }
     }


    if( strcmp(argv[1],"NetworkFilter")==0)
    {
       
        FILE * ip =  fopen(argv[2],"r");
        FILE * ip1 = fopen(argv[3],"r");
        FILE * op = fopen(argv[4],"w+");


       Test.ReadNetWork(ip1); 
       Test.ReadFeature(ip); 
       Test.NetWorkFilteFeature();
       Test.OutPutFeatureS(op);       

 
    }
     
   
    if( strcmp(argv[1],"Help")==0)
	{
	    printf("FeatureGen ToFeature   InputFile OutputFile AverageCorrelation AverageComponentRatio\n"); 
            printf("	'ToFeature' is the keyword.\n");
	    printf("	InputFile is a file which is generated from 'Cluster' commond\n");
            printf("	OutputFile is a output feature record file\n");
	    printf("	AverageCorrelation: The threshold of Average correlation\n"); 
	    printf("	AverageComponentRatio: The threshold of Average Component Ratio\n");

            printf("FeatureGen Decipher FeatureFile   InputFile OutputFile\n");
	    printf("	'Decipher' is the keyword\n");
            printf("	FeatureFile is the file which is generated from 'FeatureGen ToFeature' commond\n"); 
            printf("	InputFile is a gene expression matrix\n"); 
	    printf("	OutputFile is the sparse matrix format file\n");
	}




   return 0;
}




