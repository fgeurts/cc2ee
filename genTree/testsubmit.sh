#!/bin/bash

#totalJobs=10000 first submitted 10000 jobs, 500000 events/per job
#totalJobs=10000 second submitted 10000 jobs,2500000 events/per job

for num in {1..101..10}
do	echo "submit for the ${num}th job !!! "
	echo "job id will be: " ${num} "to" $((${num}+10))
	./submit.sh ${num} $((${num}+10))
	echo "-------------------go to sleep for 1 hour------------------------------"
	sleep 1s
done

echo "all jobs submitted"



#./submit.sh 1 1000
#echo " 1000 jobs submitted, now go to sleep for 4 hours"
#sleep 2h
#
#./submit.sh 1001 2000 
#echo " second 1000 jobs submitted, now go to sleep for 4 hours"
#sleep 2h
#
#./submit.sh 2001 3000 
#echo " third 1000 jobs submitted, now go to sleep for 4 hours"
#sleep 2h
