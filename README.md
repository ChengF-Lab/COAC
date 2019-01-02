# compile

     make

When 'make' command is fininshed,'FeatureGen' and 'Cluster' files are generated.    

# config file 
  The parameters for procedure are recorded in config file. 
  
The content of config file:
    
    Tolerance = 0.000001
    CutOff = 0.20
    ZeroTolerance = 0.03
    MaxDepth =  1024
    MinSupport = 3

Tolerance:

   ![first equation](http://latex.codecogs.com/gif.latex?%5Csigma%20_%7Bi%7D) is an eigen value of matrix.  
   
   If Tolerance is specified, threshold = 1-Tolerance; The procedure will find the maxmum r which satisfy  
   ![first equation](http://latex.codecogs.com/gif.latex?%5Cfrac%7B%5Csigma%20_%7B1%7D%5E%7B2%7D&plus;%5Csigma%20_%7B2%7D%5E%7B2%7D&plus;%5Csigma%20_%7B3%7D%5E%7B2%7D%20...%20&plus;%5Csigma%20_%7Br%7D%5E%7B2%7D%7D%7B%5Csigma%20_%7B1%7D%5E%7B2%7D&plus;%5Csigma%20_%7B2%7D%5E%7B2%7D&plus;%5Csigma%20_%7B3%7D%5E%7B2%7D%20...%20&plus;%5Csigma%20_%7Bn%7D%5E%7B2%7D%7D%5Cleq%20threshold).Then the first r principal components are selected. 
 
CutOff: 
   
   There a F-distribution for all attribute values in a component.The P-value of each selected attribute caculated by F-distribution must be smaller than CutOff. So the number of selected attribute in each component is depended on CutOff. Smaller CutOff can shrink the number of final result.If the number of result exponentially expand too large,the CutOff should be tuned smaller.
   
MaxDepth:

   The maximum number of principal components of Attribute Clustering. 

MinSupport:

   The minmum size of Attribute Clustering.The bigger MinSupport will introduce smaller the number of final result.
   
   
# Commands

    ./Cluster TestMatrix   OutFile   [ FilterFile ] 
    
The config file is necessary that exist in the same directory with data file. if not,the path of config file must be specified.

    ./Cluster TestMatrix  OutFile  [ FilterFile ]  -config   configPATH/config 
    
'FilterFile' is optional, the column name which not be contained in this file will be ignored,the format of the 'FilterFile' is described as follow: 

    column_name1
    column_name2
    ...
    column_namen
         
'TestMatrix' is the input file, the format of the input file is described as follow: 
       
       column_name1  column_name2  column_name3  column_name4 ...
        0.12         0.11           0.34          0.98        ...
        ... 
  
 'OutFile' is output file name. 
 
        grep CorrelationAverage OutFile > log 
  
 After above command is executed,the content of log file is described as following: 
 
        CorrelationAverage: 3 0.999124 0.148882
        CorrelationAverage: 3 0.989092 0.142556
        ......
  The first column numeric indicate the number of elements in each clustering.The second column numeric indicate Average correlation of this clustering and the third column numeric indicate Average Component Ratio of this clustering.  
  The  distribution of Average correlation and Average Component Ratio can be obtained from log file,The thresholds of Average correlation and Average Component Ratio for filtration are evaluated from this distribution. 
 
      
     FeatureGen ToFeature   InputFile OutputFile AverageCorrelation AverageComponentRatio
 
 'ToFeature' is the keyword. InputFile is the filename of Cluster commond generated file.  OutputFile is output file name for features record.AverageCorrelation is a numeric indicate the threshold of Average correlation and AverageComponentRatio is a  numeric indicate the threshold of Average Component Ratio. this commond will generate the feature file. The format of feature file is described as follows: 
            
      3       MIR1244-1,MIR1244-3,MIR1244-2,  1       (0.577350,0.577350,0.577350,)
 
 The first numeric column  indicate the number of elements in each clustering.The next columns indicate the names of elements. the second numeric column indicate the number of components of this clustering and the next are the vectors for each component.
 
 
    FeatureGen Decipher FeatureFile   InputFile OutputFile
   
  'Decipher' is the keyword,FeatureFile is the filename of 'FeatureGen ToFeature' commond generated file.InputFile is a matrix file with the same format with the input file of 'Cluster' command, OutputFile is the filename of a sparse matrix which is a mapping of outputfile. This mapping is implemented using features which recorded in FeatureFile.The OutputFile with sparse matrix format can be loaded into R environment directly for further analysis.     
  
    FeatureGen NetworkFilter  FeatureFile  NetworkFile OutputFile

  'NetworkFilter' is the keyword,NetworkFile is two column file which indicate such as protein-protein interaction or protein-gene interaction. FeatureFile file is the same as input feature file of "FeatureGen Decipher" command. OutputFile is the feature file after fitering. The filtering rule is definited that if elements in a clustering is connect to each another in the network which  is provided in the network file, then this clustering is kept. 


  # Optional Commands
  
    FeatureGen Statistics InputFile InputFileGroupTag OutputFile
    
  'Statistics' is the keyword.This command is used to caculate the mean and the variance for each group. Also this can be caculated with R scripts, as the number of features increasing,it is time consuming.  InputFile is a sparse matrix format file, InputFileGroupTag is a file to record the group tags. The output are two columnes, For a feature, one is the sum of feature value in each group,and another is the sum of square of feature value in each group.    
   
    InputFileGroupTag:
    
     s1     group1 
     s2     group1
     s3     group2
       ....

    OutputFile:
                   sum of feature value       sum of square of feature value
       group1            232.3                           534.21
    F1 group2            182.24                          832.11
       group3            122.11                          77.9
       .....
    
       group1           142.1                            88.76
    F2 group2           265.11                           12.22
       group3           761.1                            102.11
       ...
       
       ....... 
    
  
     FeatureGen CorStatistics InputFile InputFileNumericalTag OutputFile
     
   'CorStatistics' is the keyword, InputFileNumericalTag record a numerical vector indicate the value of each record, And OutputFile   record the correlation between the numerical vector and each features. 
   
    OutputFile 
    F1  0.2
    F2  0.6
    F3  0.2
    ....
    
  
  
  
  

