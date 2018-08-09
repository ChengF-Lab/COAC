# sub-network
#config file 
  The parameters for procedure are recorded in config file. 
  
The content of config file:
    
    Tolerance = 0.000001
    CutOff = 0.20
    ZeroTolerance = 0.03
    MaxDepth =  1024
    MinSupport = 3
    CutThresholdRatio = 0.0

Tolerance is defined as follow:

   ![first equation](http://latex.codecogs.com/gif.latex?%5Csigma%20_%7Bi%7D) is an eigen value of matrix.  
   
   If Tolerance is specified, threshold = 1-Tolerance; The procedure will find the maxmum r which satisfy  
   ![first equation](http://latex.codecogs.com/gif.latex?%5Cfrac%7B%5Csigma%20_%7B1%7D%5E%7B2%7D&plus;%5Csigma%20_%7B2%7D%5E%7B2%7D&plus;%5Csigma%20_%7B3%7D%5E%7B2%7D%20...%20&plus;%5Csigma%20_%7Br%7D%5E%7B2%7D%7D%7B%5Csigma%20_%7B1%7D%5E%7B2%7D&plus;%5Csigma%20_%7B2%7D%5E%7B2%7D&plus;%5Csigma%20_%7B3%7D%5E%7B2%7D%20...%20&plus;%5Csigma%20_%7Bn%7D%5E%7B2%7D%7D%5Cleq%20threshold).Then the first r principal components are selected. 
 
CutOff is defined as follow: 
   
   There a F-distribution for all attribute values in a component.The P-value of each selected attribute caculated by F-distribution must be smaller than CutOff. So the number of selected attribute in each component is depended on CutOff. Small CutOff can shrink the number of final result.If the number of result exponentially expand too large,the CutOff should be tuned smaller.
   
   
 
  
# Comands
