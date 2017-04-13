#!/bin/sh
for i in "adv" "diff"; do
	for j in "quad" "tri"; do
		for k in "20" "40" "60" "80"; do
			for m in "1" "2" "3"; do
				mpirun -n 2 -host localhost ../../../bin/conv2d_st -run "${i}_${j}_${k}_${m}.inc" > "${i}_${j}_${k}_${m}.log"
				echo "finish ${i}_${j}_${k}_${m} test"
			done
		done
	done
done