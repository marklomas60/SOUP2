#!/bin/bash
#$ -l h_rt=48:00:00

#How many qsubs
nprocesses=500
dir=/data/sm1epk/SDGVM_runs
#Input file for run.Several values should be set to ARGUMENT
inputfile=$dir/test_two_crops_par.dat
#Output final directory
outputdir=$dir/test_two_crops_par
#File with sites at resolution
sites=/home/sm1epk/SDGVM/SDGVM_misc/SDGVM/data/sites/global_30min.dat

#Make temporary output folder
tmpdir=$dir/tempoutput
rm -fr $tmpdir
mkdir -p $tmpdir

#Count sites and adjust number of processes
count=`wc -l < $sites`
if [ $nprocesses -gt $count ]
then
  nprocesses=$count
fi

#How many sites per process
chunk=$[($count-1)/$nprocesses+1]
echo sites=$count processes=$nprocesses chunk=$chunk

#Create batch file and folder for each chunk
#Creates the strings needed for the prompt
for i in $(seq 1 $nprocesses)
do
  outd=$tmpdir/r$i
  batchfile=$tmpdir/batchh-$i
  mkdir $outd
  let start=($i-1)*$chunk+1
  let end=$i*$chunk
  if [ $end -gt $count ]
  then
    end=$count
  fi
  echo "#!/bin/bash">> $batchfile
  echo "#$ -l h_rt=08:00:00">> $batchfile
  echo "module load compilers/gcc/5.2">> $batchfile
  echo "./sdgvm.exe $inputfile $outd $sites $start $end">> $batchfile
done

#Launch the jobs
for i in $(seq 1 $nprocesses)
do
  prunning=`Qstat | grep 'batchh-'|grep 'sm1epk'| wc -l`
  while [ $prunning -ge $nprocesses ]
  do
    sleep 1
    prunning=`Qstat | grep 'batchh-'|grep 'sm1epk'| wc -l`
  done
  echo submitting job $i
  qsub $tmpdir/batchh-$i 1> $tmpdir/output-$i 2> $tmpdir/error-$i
  sleep 1
done 

#Wait for final job to finish
sleep 30
prunning=`Qstat | grep 'batchh-'|grep 'sm1epk'| wc -l`
while [ $prunning -gt 0 ]
do
  sleep 30
  prunning=`Qstat | grep 'batchh-'|grep 'sm1epk'| wc -l`
done  
  
#Create output folder
mkdir -p $outputdir    
rm -f $outputdir/*

#Removes first line from site_info.dat file
#It just keeps the header from the 1st one
for i in $(seq 2 $nprocesses)
do
  sed 1d $tmpdir/r$i/site_info.dat > $tmpdir/r$i/temp
  cp $tmpdir/r$i/temp $tmpdir/r$i/site_info.dat
  rm $tmpdir/r$i/temp
done
 
#Gets the names of the output files
files=`ls $tmpdir/r1`

#Merges output files
for i in $(seq 1 $nprocesses)
do
  for file in $files
  do
    cat $tmpdir/r$i/$file >> $outputdir/$file
  done
done
  




