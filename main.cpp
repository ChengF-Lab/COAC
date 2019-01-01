#include  <iostream>
#include  "stdio.h"
#include  <cmath>
#include  "getCommonSeq/ComSeq.h"
using namespace std;

int main( int argc , char *argv[] )
{

   ComSeq test;
   char buffer[1024];

   FILE * ip =  fopen(argv[1],"r");
   FILE * op = fopen(argv[2],"w+");
   FILE * filterIp;  
   FILE * configip=NULL;



 //  printf("%ld",std::numeric_limits<size_t>::max());

   if( argc>2 )  
     filterIp = fopen(argv[3],"r");


  

    for(int i =0;i<argc-1;i++)
    {
       if( strcmp(argv[i],"-config")==0)
       {
           configip = fopen(argv[i+1],"r");
       }
    }

    if( configip == NULL)
    {
        configip = fopen("config","r");
    }



    if( strcmp(argv[1],"Help")==0)
	{
	    printf("./Cluster TestMatrix   OutFile   [ FilterFile ] \n"); 
       printf("   The config file is necessary that exist in the same directory with data file. if not,the path of config file must be specified.\n");
       printf("./Cluster TestMatrix  OutFile  [ FilterFile ]  -config   configPATH/config\n"); 
       printf("TestMatrix is a input matrix \n"); 
	}




   if((ip == NULL)||(op== NULL))
   {
      printf("usage: ./cluster MatrixFile  OuputFile  [ FilterFile ]   [ -config  ConfigurationFile ]");
      printf("\n");
      return 0;
   }

    if( filterIp != NULL)
    {
       test.ReadFilteFile(filterIp);
    }
   //test.matrixReadWithName(ip);
   test.LimitedRank = 1024;
   test.tolerance = 0.001;
   test.range = 0.12;
   test.rankSkip = true;
   test.Centering = false;
// test.matrixReadWithName(ip);
//   M_d.close();


   if(configip!= NULL)
   {


       while (!feof(configip))
       {
               fscanf(configip,"%[\b|\t]*",buffer);
               fscanf(configip,"%s",buffer);

               if(strcmp(buffer,"Tolerance")==0)
               {
                   fscanf(configip,"%[\b|\t]*",buffer);
                   fscanf(configip,"%s",buffer);

                   if(strcmp(buffer,"=")==0)
                   {
                      fscanf(configip,"%lf",&test.tolerance);
                   }
                   fscanf(ip,"%[\b|\t|\n]",buffer);
                   continue;
               }

               if(strcmp(buffer,"ZeroTolerance")==0)
               {
                   fscanf(configip,"%[\b|\t]*",buffer);
                   fscanf(configip,"%s",buffer);

                   if(strcmp(buffer,"=")==0)
                   {
                      fscanf(configip,"%lf",&test.NoZeroRatio);
                   }
                   fscanf(configip,"%[\b|\t|\n]",buffer);
                   continue;
               }
              if(strcmp(buffer,"CutOff")==0)
              {
                   fscanf(configip,"%[\b|\t]*",buffer);
                   fscanf(configip,"%s",buffer);

                   if(strcmp(buffer,"=")==0)
                   {
                      fscanf(configip,"%lf",&test.Log_Pvalue);
                      test.Log_Pvalue = log(test.Log_Pvalue);

                   }
                   fscanf(configip,"%[\b|\t|\n]",buffer);
                   continue;
               }


                 if(strcmp(buffer,"MaxDepth")==0)
               {
                   fscanf(configip,"%[\b|\t]*",buffer);
                   fscanf(configip,"%s",buffer);

                   if(strcmp(buffer,"=")==0)
                   {
                      fscanf(configip,"%d",&test.MaxDepth);
                   }
                   fscanf(configip,"%[\b|\t|\n]",buffer);
                   continue;
               }


                  if(strcmp(buffer,"MinSupport")==0)
               {
                   fscanf(configip,"%[\b|\t]*",buffer);
                   fscanf(configip,"%s",buffer);

                   if(strcmp(buffer,"=")==0)
                   {
                      fscanf(configip,"%d",&test.minsupport);
                   }
                   fscanf(configip,"%[\b|\t|\n]",buffer);
                   continue;
               }

                    
                 if(strcmp(buffer,"CutThresholdRatio")==0)
               {
                   fscanf(configip,"%[\b|\t]*",buffer);
                   fscanf(configip,"%s",buffer);
                   test.CutThresHold = 0.0;

                   if(strcmp(buffer,"=")==0)
                   {
                      fscanf(configip,"%lf",&test.CutThresHold);
                   }
                   fscanf(configip,"%[\b|\t|\n]",buffer);
                   continue;
               }











       }

   }

   test.matrixReadWithName(ip);
   printf("Matrix have been read.\n");
   printf("n_rows: %d\n",(int)test.Sample.n_rows);
   printf("n_cols: %d\n",(int)test.Sample.n_cols);
   test.filter();
   printf("Matrix have been filtered.\n");
   printf("n_rows: %d\n",(int)test.Sample.n_rows);
   printf("n_cols: %d\n",(int)test.Sample.n_cols);
   printf("n_cols: %d\n",(int)test.SampleName.size());
   test.AttributeCluster();
   test.OutClusterWithNameCombact(op);
 

   return 0;

}




