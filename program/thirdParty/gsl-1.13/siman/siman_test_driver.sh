#! /bin/sh

# assume good result from tests; increment it if any test fails
EXIT_STATUS=0

for seed in "" 12345 ;
do 
./siman_test > siman_test.out 2>&1
SECOND_LAST_ENERGY=`tail -2 siman_test.out1 | head -1 | awk '{print $4}'`
LAST_ENERGY=`tail -1 siman_test.out1 | awk '{print $4}'`
# echo " " $SECOND_LAST_ENERGY $LAST_ENERGY
if [ $SECOND_LAST_ENERGY = $LAST_ENERGY ];
then
    echo -n "PASS: "
else
    echo -n "FAIL: "
    EXIT_STATUS=`expr $EXIT_STATUS + 1`
fi
echo "simulated annealing test (travelling salesman problem) seed=${seed:-default}"
done

exit $EXIT_STATUS
