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




   return 0;
}




