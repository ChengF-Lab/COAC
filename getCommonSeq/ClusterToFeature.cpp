#include "ComSeq.h"
#include "stdio.h"
#include "math.h"
#include <fstream>
#include <queue>
#include <boost/math/distributions/fisher_f.hpp>
using namespace std;
using namespace arma;
using namespace boost::math;


char buffer[1024];


 void ClusterToFeature::ReadSpareM( std::vector < std::vector < unsigned int > > &Index,std::vector < std::vector < double > > &Value,FILE *Op)
 {
     Index.resize(1);
     Value.resize(1);
     char buffer[1024];
     int I =0;
     while((!feof(Op)))
     {
        int index_=-1;
        double Value_=0.0;

        fscanf(Op,"%[\b|\t]*",buffer);
        fscanf(Op,"%[  ]",buffer);

        if(fscanf(Op,"%[\n]",buffer))
        {
            if(Index[I].size()!=0)
            {
                Index.resize(Index.size()+1);
                Value.resize(Value.size()+1);
                I=I+1;
            }

           fscanf(Op,"%[\b|\t|\n]*",buffer);
            continue;
        }

        fscanf(Op,"%d",&index_);

        fscanf(Op,"%[:]",buffer);
        if(fscanf(Op,"%lf",&Value_))
        {
          Index[I].push_back(index_);
          Value[I].push_back(Value_);
        }
     }


      while(Index.size()!=0&&Index[Index.size()-1].size()==0)
      {
          Index.resize(Index.size()-1);
          Value.resize(Value.size()-1);
      }

 }

void ClusterToFeature::StatisticSpareMGroup( std::vector < std::vector < unsigned int > > &Index,std::vector < std::vector < double > > &Value, std::vector<unsigned int> & GroupFlag,std::vector < double > &SumValue,std::vector < double > &SquareValue)
{
    int Max =0;

    for(int i=0;i<GroupFlag.size();i++)
    {
       if(GroupFlag[i]>Max)
       {
           Max = GroupFlag[i];
       }
    }

    int MaxIndex=0;

    for(int j=0;j<Index.size();j++)
      for(int i=0;i<Index[j].size();i++)
        {
            if(Index[j][i]>MaxIndex)
            {
               MaxIndex=Index[j][i];
            }
        }
    SquareValue.resize(MaxIndex*(Max+1));
    SumValue.resize(MaxIndex*(Max+1));

    for(int i=0;i<MaxIndex*Max;i++)
    {
         SquareValue[i]=0;
         SumValue[i]=0;
    }


    for(int j=0;j<Index.size();j++)
     for(int i=0;i<Index[j].size();i++)
       {

          SquareValue[ (Index[j][i]-1)*(Max+1)+GroupFlag[j] ] = SquareValue[Index[j][i]-1]+Value[j][i]*Value[j][i];
          SumValue[ (Index[j][i]-1)*(Max+1)+GroupFlag[j] ] = SumValue[Index[j][i]-1] + Value[j][i];
       }


}

 void ClusterToFeature::StatisticSpareM( std::vector < std::vector < unsigned int > > &Index,std::vector < std::vector < double > > &Value, std::vector < double > &SumValue,std::vector < double > &SquareValue)
 {
     int MaxIndex=0;
     for(int j=0;j<Index.size();j++)
      for(int i=0;i<Index[j].size();i++)
        {
            if(Index[j][i]>MaxIndex)
            {
               MaxIndex=Index[j][i];
            }
        }
     SquareValue.resize(MaxIndex);
     SumValue.resize(MaxIndex);

     for(int i=0;i<MaxIndex;i++)
     {
         SquareValue[i]=0;
         SumValue[i]=0;
     }



      for(int j=0;j<Index.size();j++)
       for(int i=0;i<Index[j].size();i++)
        {

            SquareValue[Index[j][i]-1] = SquareValue[Index[j][i]-1]+Value[j][i]*Value[j][i];
            SumValue[Index[j][i]-1] = SumValue[Index[j][i]-1] + Value[j][i];
        }




 }

void ClusterToFeature::StatisticSpareCor( std::vector < std::vector < unsigned int > > &Index,std::vector < std::vector < double > > &Value, std::vector< double > & CorFlag,std::vector < double > & CorValue)
{

     int MaxIndex=0;
     double CorFlagN=0;
     for(int j=0;j<Index.size();j++)
      for(int i=0;i<Index[j].size();i++)
        {
            if(Index[j][i]>MaxIndex)
            {
               MaxIndex=Index[j][i];
            }
        }
    std::vector < double > ValueN;
    CorValue.resize(MaxIndex);
    ValueN.resize(MaxIndex);
      for(int i=0;i<MaxIndex;i++)
     {
        CorValue[i]=0;
        ValueN[i] =0;
     }

      for(int j=0;j<Index.size();j++)
       for(int i=0;i<Index[j].size();i++)
        {
            CorValue[Index[j][i]-1] = CorValue[Index[j][i]-1]+Value[j][i]*CorFlag[j];
            ValueN[Index[j][i]-1] = ValueN[Index[j][i]-1] + Value[j][i]*Value[j][i];
        }

      for(int i=0;i<CorFlag.size();i++)
      {
        CorFlagN = CorFlagN + CorFlag[i]*CorFlag[i];
      }

     for(int i=0;i<MaxIndex;i++)
     {
        CorValue[i]=CorValue[i]/(sqrt(ValueN[i])+0.000000000001);
        CorValue[i]=CorValue[i]/(sqrt(CorFlagN)+0.00000000001);
     }
}

void  ClusterToFeature::OutPutFeatureS(FILE *Op)
{
		char buffer[1024];
		int Size;
	    int Array_Size;

		 Size = FeatureNameArray.size();
		 for(int i=0;i<Size;i++)
		 {

				fprintf(Op,"%d\b\t",(int)FeatureNameArray[i].size());
				for( int j=0;j<FeatureNameArray[i].size();j++)
				{
				   fprintf(Op,"%s,",FeatureNameArray[i][j].c_str());
                }
				   fprintf(Op,"\b\t");
				   fprintf(Op,"%d\b\t",(int)FeatureM[i].n_rows);
				   for(int i1=0;i1<FeatureM[i].n_rows;i1++)
				   {
							fprintf(Op,"(");
							for(int j=0;j<FeatureM[i].n_cols;j++)
						  {
							 fprintf(Op,"%lf,",FeatureM[i](i1,j));
						  }
							fprintf(Op,")");
				   }
							fprintf(Op,"\n");
		 }
}





void ClusterToFeature::OutPutFeature()
{
  unsigned int Size;
  double AC;
  double PC =0;
  int Rank;

     Size = FeatureNameArray.size();
     FeatureM.resize(Size);

     for(int i=0;i<Size;i++)
     {


            svd_econ(U,s,V,FeatureComponentArray[i]);

            AC = 0;
            PC = 0;

            for(int i=0;i<s.n_elem;i++)
            {
                 AC = AC+s[i]*s[i];
            }
             Rank=0;
             for(int i=0;i<s.n_elem;i++)
            {
                PC = PC+s[i]*s[i];
                Rank ++;

                if(PC>((1.0-tolerance)*(1.0-tolerance)*AC))
                {
                   break;
                }
             }

             FeatureM[i].resize(Rank,FeatureComponentArray[i].n_cols);

             for(int i1=0;i1<Rank;i1++)
             {
                 for(int j=0;j<FeatureComponentArray[i].n_cols;j++)
                 {
                    FeatureM[i](i1,j)=V(j,i1)*(s[i1]/s[0]);
                 }
             }
     }
}

void ClusterToFeature::SampleDecipher(FILE *Op)
{
     int Size = FeatureNameArray.size();
     std::vector< std::vector<int> >featureIndex;
     std::vector< string >::iterator it;
     featureIndex.resize(Size);


     for(int i=0;i<Size;i++)
     {
          for(int j=0;j<FeatureNameArray[i].size();j++)
          {
              it = std::find(SampleName.begin(),SampleName.end(),FeatureNameArray[i][j].c_str());

              if(it != SampleName.end())
                featureIndex[i].push_back(it-SampleName.begin());
          }

           if(featureIndex[i].size()!=FeatureNameArray[i].size())
          {
               featureIndex[i].resize(0);
          }

     }

     for(int i=0;i<Sample.n_rows;i++)
     {
         for(int j=0;j<featureIndex.size();j++)
         {

                double value2 =0.0;
                double value3 =0.0;

                if(featureIndex[j].size()==0)
                 continue;



                for(int j2=0;j2<featureIndex[j].size();j2++)
                {
                    value3 = value3+Sample(i,featureIndex[j][j2])*Sample(i,featureIndex[j][j2]);
                }



                for(int j1=0;j1<FeatureM[j].n_rows;j1++)
                {
                      double value1 =0;
                      double  coexpressP =1.0;
                      double  coexpressN =1.0;
                      unsigned int P_NUM =0;

                      for(int j2=0;j2<featureIndex[j].size();j2++)
                     {
                        value1 = value1 + Sample(i,featureIndex[j][j2])*FeatureM[j](j1,j2);
                        if(FeatureM[j](j1,j2)>0)
                        {
                            coexpressP = coexpressP*Sample(i,featureIndex[j][j2]);
                            P_NUM++;
                        }else
                         {
                            coexpressN = coexpressN*Sample(i,featureIndex[j][j2]);
                         }
                     }

                      coexpressN = pow(coexpressN*coexpressN,(1.0/((featureIndex[j].size()-P_NUM)*1.0)));
                      coexpressP = pow(coexpressP*coexpressP,(1.0/(P_NUM*1.0)));

                      if(P_NUM==0)
                      {
                         coexpressP =0;
                      }
                      if(featureIndex[j].size()==P_NUM)
                      {
                         coexpressN =0;
                      }


                      if( (value3!=0)&(((coexpressN+coexpressP)/value3)>0.0001))
                        value2 = value1*value1+value2;
                }





                 if(value3!=0&&value2!=0)
                  {
                     fprintf(Op,"%d",j+1);
                     fprintf(Op,":%lf",sqrt(value2/value3));
                     fprintf(Op,"  ");

                   }


           }
                fprintf(Op,"\n");
     }



}






void ClusterToFeature::OutPutFeature(FILE *Op)
{
  unsigned int Size;
  double AC;
  double PC =0;
  int Rank;

  Size = FeatureNameArray.size();

     for(int i=0;i<Size;i++)
     {


         fprintf(Op,"%d\b\t",(int)FeatureNameArray[i].size());

         for( int j=0;j<FeatureNameArray[i].size();j++)
         {
             fprintf(Op,"%s,",FeatureNameArray[i][j].c_str());
         }

          fprintf(Op,"\b\t");

            svd_econ(U,s,V,FeatureComponentArray[i]);

            AC = 0;
            PC = 0;

            for(int i=0;i<s.n_elem;i++)
            {
                 AC = AC+s[i]*s[i];
            }
             Rank=0;
             for(int i=0;i<s.n_elem;i++)
            {
                PC = PC+s[i]*s[i];
                Rank ++;

                if(PC>((1.0-tolerance)*(1.0-tolerance)*AC))
                {
                   break;
                }
             }

           fprintf(Op,"%d\b\t",Rank);

             for(int i1=0;i1<Rank;i1++)
             {

                 fprintf(Op,"(");
                 for(int j=0;j<FeatureComponentArray[i].n_cols;j++)
                 {
                    fprintf(Op,"%lf,",V(j,i1)*((s[i1]/s[0])));
                 }

                 fprintf(Op,")");
             }

                fprintf(Op,"\n");
     }

}


void   ClusterToFeature::NetWorkFilteFeature ()
{
    std::vector<int> FilterIndex;
    std::queue<int>  Queue;
    std::vector<int> NetworkNameID;
    std::vector < std::vector <std::string> > FeatureNameArrayNew;
    std::vector < mat> FeatureMNew;
    std::vector<int> key;


	key.resize(2);
    if( NetworkDictionary.size()==0 )
    {
       return;
    }






    for(int i=0;i<FeatureNameArray.size();i++)
    {
         int j;
         bool FILTER = false;
		 int FeatureNameArraySize = 0;

         if( NetworkNameID.size() < FeatureNameArray[i].size())
         {
             NetworkNameID.resize(FeatureNameArray[i].size()+1);
         }

         for(j=0;j<FeatureNameArray[i].size();j++)
         {

            if ( NetworkName.find( FeatureNameArray[i][j].c_str()) ==  NetworkName.end() )
            {
                break;

            }else
            {
                NetworkNameID[j] = NetworkName[FeatureNameArray[i][j].c_str()];
            }
         }
		 
         FeatureNameArraySize = FeatureNameArray[i].size();


         if( j!=FeatureNameArray[i].size())
         {
             continue;

         }else{


            Queue.push(NetworkNameID[0]);
            NetworkNameID[0] =-1;

            while(!Queue.empty())
            {

                int RecordedIndexA;
                int RecordedIndexB;

                RecordedIndexA = Queue.front();
				Queue.pop();
                for(int i=0;i<FeatureNameArraySize;i++)
                {
                    if ( NetworkNameID[i]!=-1)
                    {
                       RecordedIndexB = NetworkNameID[i];



                      if(RecordedIndexA>RecordedIndexB)
                        {
                            key[0] = RecordedIndexA;
                            key[1] = RecordedIndexB;

                        }else
                        {
                            key[0] = RecordedIndexB;
                            key[1] = RecordedIndexA;
                        }



                       if( NetworkDictionary.find(key) != NetworkDictionary.end())
                       {
                           Queue.push(RecordedIndexB);
                           NetworkNameID[i]=-1;
                       }
                    }
                }
            }

            for(int i1=0;i1<FeatureNameArray[i].size();i1++)
            {
                if( NetworkNameID[i1]!=-1)
                {
					FILTER =true;
                    break;
                }
            }

            if( !FILTER )
            FilterIndex.push_back(i);
         }
    }
    for(int i=0;i<FilterIndex.size();i++)
    {
          FeatureNameArrayNew.push_back(FeatureNameArray[FilterIndex[i]]);
          FeatureNameArray[FilterIndex[i]].resize(0);
    }

    for(int i=0;i<FilterIndex.size();i++)
    {
         FeatureMNew.push_back(FeatureM[FilterIndex[i]]);
    }
      FeatureNameArray.resize(0);
      FeatureNameArray = FeatureNameArrayNew;
      FeatureM.resize(0);
      FeatureM = FeatureMNew;

}



void  ClusterToFeature::ReadNetWork(FILE *ip)
{
    std::string  IDPath;
    char buffer[1024];
    int Size;
    int Array_Size;
    int RecordedIndexA;
    int RecordedIndexB;
    std::vector <int> key;

	key.resize(2);



    while ((!feof(ip)))
    {
        fscanf(ip,"%[\b|\t]*",buffer);
        fscanf(ip,"%s",buffer);
        IDPath.assign(buffer);
        fscanf(ip,"%[\b\t]*",buffer);



        if( NetworkName.find(IDPath)== NetworkName.end())
       {

            NetworkName.insert((std::make_pair(IDPath,NetworkName.size()-1)));

       }
       {
            RecordedIndexA = NetworkName[IDPath];
       }

       fscanf(ip,"%s",buffer);
       IDPath.assign(buffer);
       fscanf(ip,"%[\b|\t]*",buffer);



       if( NetworkName.find(IDPath)== NetworkName.end())
       {

            NetworkName.insert((std::make_pair(IDPath,NetworkName.size()-1)));

       }
       {
            RecordedIndexB = NetworkName[IDPath];
       }

       if(RecordedIndexA>RecordedIndexB)
       {
          key[0] = RecordedIndexA;
          key[1] = RecordedIndexB;

       }else
       {
          key[0] = RecordedIndexB;
          key[1] = RecordedIndexA;
       }

    //  fscanf(ip,"%s",buffer);
        IDPath.assign("1");
        fscanf(ip,"%[\b|\t|\n]*",buffer);


        if( NetworkDictionary.find(key) == NetworkDictionary.end())
        {
           NetworkDictionary.insert(std::make_pair(key,IDPath));
        }
    }

}



void ClusterToFeature::ReadFeature(FILE *ip)
{
   char buffer[1024];
   int Size;
   int Array_Size;


   while((!feof(ip)))
  {
     fscanf(ip,"%[\b|\t|\n]*",buffer);
     fscanf(ip,"%d",&Size);
     fscanf(ip,"%[\b|\t|\n]*",buffer);

     std::vector <std::string> FeatureName;

     FeatureName.resize(Size);



     for(int i=0;i<Size;i++)
      {
           fscanf(ip,"%[\b|\t|\n]*",buffer);
           fscanf(ip,"%[^,]",buffer);
           FeatureName[i].assign(buffer);
           fscanf(ip,",",buffer);
      }

        FeatureNameArray.push_back(FeatureName);

       fscanf(ip,"%[\b|\t|\n]*",buffer);
       fscanf(ip,"%d",&Array_Size);
       fscanf(ip,"%[\b|\t|\n]*",buffer);

       mat FeatureComponent;
       FeatureComponent.resize(Array_Size,Size);

       for(int j=0;j<Array_Size;j++)
       {
            double Value;

            fscanf(ip,"%[\b|\t|\n]*",buffer);
            fscanf(ip,"%[(]",buffer);
            for(int i=0;i<Size;i++)
            {
                fscanf(ip,"%[\b|\t|\n]*",buffer);
                fscanf(ip,"%lf,",&Value);
                FeatureComponent(j,i)=Value;
            }

           fscanf(ip,"%[\b|\t|\n]*",buffer);
           fscanf(ip,"%[)]",buffer);
      }
         FeatureM.push_back(FeatureComponent);
}
}

void ClusterToFeature::ClusterRead( FILE *ip)
{


  int Size;
  int FeatureNameSize;
  int FeatureSize;
  int ComponentSize;

  std::vector <int> IDPath;



  while((!feof(ip)))
  {

      fscanf(ip,"%[\b|\t|\n]*",buffer);

      if(fscanf(ip,"Filtered Total Size:%d",&Size))
      {
            printf("%d\n",Size);
            break;
      }
      {
          fscanf(ip,"%[^\n]\n",buffer);
      }
  }

   rewind(ip);

   fscanf(ip,"%[^\n]\n",buffer);

      for(int i=0;i<Size;i++)
      {

        IDPath.resize(0);
        fscanf(ip,"%[\b|\t|\n]*",buffer);
        fscanf(ip,"%d",&FeatureSize);
        std::vector <std::string> FeatureName;
        FeatureName.resize(FeatureSize);
        mat FeatureComponent;

            for(int j=0;j<FeatureSize;j++)
            {
               fscanf(ip,"%s",buffer);
               FeatureName[j].assign(buffer);
            }

        fscanf(ip,"%[\b|\t|\n]*",buffer);
        fscanf(ip,"%d",&ComponentSize);
        FeatureComponent.resize(ComponentSize,FeatureSize);
        for(int j=0;j<ComponentSize;j++)
        {
            int CompID;
            fscanf(ip,"%[\b|\t|\n]*",buffer);
            fscanf(ip,"%d",&CompID);
            if(CompID<0)
            {
              CompID = -CompID;
            }
            IDPath.push_back(CompID);
             for(int j2=0;j2<FeatureSize;j2++)
             {
               double value;

                fscanf(ip,"%[\b|\t|\n]*",buffer);
                fscanf(ip,"%lf",&value);
                FeatureComponent(j,j2)=value;
             }



        }
            double CorAvarage = 0;
            double ComponentRatioAvarage =0;
            int AttributeNum;

            fscanf(ip,"%[\b|\t|\n]*",buffer);

            int AdditionalSize =0;

          while(!feof(ip))
          {
            if(fscanf(ip,"CorrelationAverage: %d",&AttributeNum))
            {
                 fscanf(ip,"%[\b|\t]*",buffer);
                 fscanf(ip,"%lf",&CorAvarage);
                 fscanf(ip,"%[\b|\t]*",buffer);
                 fscanf(ip,"%lf",&ComponentRatioAvarage);
                 break;
            }

            AdditionalSize++;

           // FeatureComponent.resize(ComponentSize+AdditionalSize,FeatureSize);
            int CompID;
            fscanf(ip,"%[\b|\t|\n]*",buffer);
            fscanf(ip,"%d",&CompID);
            //IDPath.push_back(CompID);

            for(int j2=0;j2<FeatureSize;j2++)
            {
               double value;

                fscanf(ip,"%[\b|\t|\n]*",buffer);
                fscanf(ip,"%lf",&value);
               // FeatureComponent(ComponentSize+AdditionalSize-1,j2)=value;
             }


         }

         if((CorAvarage<= CorThreshold)||(ComponentRatioAvarage <= CompRThreshold))
         {
            continue;
          }


       if(ComponentTable.find(IDPath)==ComponentTable.end())
       {
            FeatureNameArray.push_back(FeatureName);
            FeatureComponentArray.push_back(FeatureComponent);

            ComponentTable.insert((std::make_pair(IDPath,FeatureNameArray.size()-1)));

       }else
       {
             int RecordedIndex = ComponentTable[IDPath];

             FeatureNameArray[RecordedIndex].insert(FeatureNameArray[RecordedIndex].end(),FeatureName.begin(),FeatureName.end());

             FeatureComponentArray[RecordedIndex].resize(FeatureComponentArray[RecordedIndex].n_rows,FeatureComponentArray[RecordedIndex].n_cols + FeatureComponent.n_cols);

             for(int j=0;j<FeatureComponent.n_cols;j++)
             {
                  for(int j2=0;j2<FeatureComponentArray[RecordedIndex].n_rows;j2++)
                   FeatureComponentArray[RecordedIndex](j2,(j+FeatureComponentArray[RecordedIndex].n_cols-FeatureComponent.n_cols)) = FeatureComponent(j2,j);
             }

      }

  }

}



void ClusterToFeature::ReadFilteFile(FILE *ip)
{
    char buffer[1024];

    while((!feof(ip)))
    {
         fscanf(ip,"%[\b|\t|\n]*",buffer);
         fscanf(ip,"%s",buffer);
         FilterSet.insert(buffer);
    }

}


void ClusterToFeature::FilteFeature()
{
    std::vector<int> FilterIndex;
    std::vector < std::vector <std::string> > FeatureNameArrayNew;
    std::vector < mat> FeatureMNew;

    if( FilterSet.size()==0)
        return;


    FeatureNameArray.size();
    FeatureM.size();

       for(int i=0;i<FeatureNameArray.size();i++)
      {
         bool UNFILTER = false;
         for(int j=0;j<FeatureNameArray[i].size();j++)
         {
            if( FilterSet.find(FeatureNameArray[i][j].c_str()) != FilterSet.end() )
            {
                UNFILTER = true;
                break;
            }
         }
        if(UNFILTER)
         {
            FilterIndex.push_back(i);
         }
      }
      for(int i=0;i<FilterIndex.size();i++)
      {
          FeatureNameArrayNew.push_back(FeatureNameArray[FilterIndex[i]]);
          FeatureNameArray[FilterIndex[i]].resize(0);
      }

      for(int i=0;i<FilterIndex.size();i++)
      {
          FeatureMNew.push_back(FeatureM[FilterIndex[i]]);
      }
      FeatureNameArray.resize(0);
      FeatureNameArray = FeatureNameArrayNew;
      FeatureM.resize(0);
      FeatureM = FeatureMNew;

}




void ClusterToFeature::SampleRead(FILE *ip)
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

}












