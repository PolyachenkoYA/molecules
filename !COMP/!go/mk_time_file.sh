#!/bin/bash
if [[ $# != 0 && $# != 3 && $# != 4 ]]; then
  echo "format:"  
  echo "script_name (start = 0) (step = 0.1) (end = 100) (file_name = T.txt)"
  exit 1;
fi

start=$1
step=$2
end=$3
Fname=$4
if [[ $# == 0 ]]; then
  Fname="T.txt"
  start=0
  step=0.1
  end=100  
elif [[ $# == 3 ]]; then
  Fname="T.txt"  
fi
touch $Fname
for i in $(seq $start $step $end)
do
   echo $i >> $Fname
done

exit 0