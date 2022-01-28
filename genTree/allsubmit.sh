#!/bin/bash

#totalJobs=10000 first submitted 10000 jobs, 500000 events/per job
#totalJobs=10000 second submitted 10000 jobs,2500000 events/per job
# since there are so many jobs no finished in second round, need to submit the third round
# 900000 events/per job

#for num in {1..10001..1000}
#for num in {10001..50001..1000}
#for num in {50001..100001..2000}


nEvtStep=1000 # step of num should be always same as nEvtStep!!!!!

for num in {1..4001..1000}
#for num in {4001..10001..1000}
#for num in {11001..20001..1000}
#for num in {21001..40001..1000}
#for num in {41001..60001..1000}
do echo "submit for the ${num}th job !!! "
	iStart=${num}
	iStop=$((${num}+${nEvtStep}))
	echo "job id will be: " ${num} "to" ${iStop}
	./submit.sh ${num} ${iStop}
	echo "-------------------go to sleep for 3 hour------------------------------"
	sleep 3h
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
