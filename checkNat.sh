#!/bin/bash

# written by Daniyar Bakir

# this file executes the code to run the program under several configurations
# write correct path to the project, and everyting will be OK

# add databases here as an array, make sure arrays of $database and $input_files have the amount of elements and correspond to each other
#input_files="vandeVijverBreast2002_5003x266.txt"
src_folder=src_experiments_C++

#input_files=("BhattacherjeeLung2001_12600x203.txt"
	     #"ChenLiver2004_10237x157.txt"
	     #"DesmedtBreast2007_22215x175.txt"
             #"Dubois2016_54675x223.txt"
	     #"NatsoulisRats2005_8491x181.txt"
	     #"Rosenwald2002_5013x203.txt"
	     #"SuCancer2001_12533x174.txt"
	     #"ValkLeukemia2004_22215x273.txt"
	     #"vandeVijverBreast2002_5003x266.txt"
	     #"WangBreast2005_22215x276.txt"
	     #"YeohLeukemia2002_5077x248.txt"
	     #"ZhanMyeloma2006_54613x234.txt")

# short names for databases
#database=("bhatta"
	  #"chen"
	  #"desmedt"
          #"dubois"
	  #"natsoulis"
	  #"rosenwald"
	  #"su"
	  #"valk"
	  #"vijver"
	  #"wang"
	  #"yeoh"
	  #"zhan")


input_files="Dubois2016_54675x223.txt"
database="dubois"

data_num=${#database[@]}

# use that just in case
#MYPATH=/path/to/project/folder/
#cd $MYPATH

#feature_size=(5 20 50 100)
#sample_size=(30 40 50 60 70 80 90 100)
#feature_size=20
#sample_size=50
feature_size=20
sample_size=50

feature_num=${#feature_size[@]}
sample_num=${#sample_size[@]}

myNP=1
myREP=10

rm ./${src_folder}/*.o
rm ./optgamma
cd $src_folder
#make mpi # to run with MPI settings
#make time # to record time during run this is not compatible with MPI should be studied more
make roc
#make # simple run
cd ..

for i_data in `seq 0 $((${data_num} - 1))`; do
    
    file_cur=${input_files[${i_data}]}
    data_cur=${database[${i_data}]}
    
    for i_sam in `seq 0 $((${sample_num} - 1))`; do
	
	sample_cur=${sample_size[${i_sam}]}
	echo sample_cur is ${sample_cur}
	
	for i in `seq 0 $((${feature_num} - 1))`; do 
	    
	    feature2=${feature_size[${i}]}
	    
	    #output_filename=${data_cur}/${data_cur}_s${sample_cur}_p${feature2}_$(date +%b_%d)
	    output_filename=bbb
	    #mylog_filename=${output_filename}.log
	    mylog_filename=aaa
	    #echo log_filename is $mylog_filename
	    
	    echo "output_filename is ${output_filename}" > ./logs/${mylog_filename}
	    echo "time is $(date +%b_%d_%Y_%H:%M:%S)" >> ./logs/${mylog_filename}
	    echo "tr=${sample_cur}" >> ./logs/${mylog_filename}
	    echo "d=${feature2}" >> ./logs/${mylog_filename}
	    echo "mr=${myREP}" >> ./logs/${mylog_filename}
	    
	    echo mpirun -np ${myNP} ./optgamma -i data/${file_cur} -o ./results/${output_filename} -tr ${sample_cur} -d ${feature2} -mr ${myREP}
	    
	    #mpirun -np ${myNP} ./optgamma -i data/${file_cur} -o ./results/${output_filename} -tr ${sample_cur} -d ${feature2} -mr ${myREP}
	    ./optgamma -i data/${file_cur} -o ./results/${output_filename} -tr ${sample_cur} -d ${feature2} -mr ${myREP}
	    
	    
	done # for i_feature
	
    done # for i_sam
    
done # for i_data
echo $database
