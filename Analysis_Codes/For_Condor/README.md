1. cmsrel CMSSW_10_5_0

2. cd CMSSW_10_5_0/src

3. Put the codes here

4. Compile the code:
   ./Makefile Anal_Sig_Wp

5. Create a log file containing the list of root files you want to run the code on. An example file is given here. 

6. Run the executable; 

   ./Anal_Sig_Wp.exe

7. It will ask for the name of the mass point, put 1000 for this example file. 

   If you want to use the full name of the log file, change L1359-1393 of Anal_Sig_Wp.C according to https://github.com/Sumantifr/Wprimetotb/blob/master/Analysis_Codes/Anal_Sig_Wp.C#L1361-L1390

   To use it in condor add the path before the file name in L1383 of Anal_Sig_Wp.C

8. Enjoy the output histograms

9. In case of an error, delete lines dealing with LHAPDF

10. Different conditions are applied depending on the year. To use the correct year, just uncomment corresponding line in L11-16
