temp=LH_Wp_Job_submit.sh
echo "#!/bin/sh "  | cat >>$temp

njobs=0
mass=0
ntotjobs=$1
while [ ${njobs} -lt ${ntotjobs} ] 
do
njobs=`expr $njobs + 1`
printf "njobs = $njobs \n"
fil1=LH_Wp_Job${njobs}
mass=`expr $njobs \* 100 + 1000`
mass=`expr $mass - 100`
echo $mass
echo "Universe = vanilla
getenv = TRUE
Executable = Anal_Sig_Wp.exe
Arguments = ${mass}
Log = ${fil1}.log
Output = ${fil1}.out
Error = ${fil1}.error
notification = never
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
#request_memory = 2000
#+MaxRuntime = 21600
Queue"| cat >>${fil1}.sh

chmod 744 ${fil1}.sh

echo "condor_submit ${fil1}.sh "  | cat >>$temp

done 

chmod 755 $temp
