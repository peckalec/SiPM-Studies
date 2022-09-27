export DIR="LongbarTiming-part2"
mkdir /home/arratialab2/fcal/$DIR/
export subDIR="LongBar_newpositions"
mkdir /home/arratialab2/fcal/$DIR/$anotherDIR

for i in {00000007..00000007}
	do
		for j in $(seq 2 1000 2) # loop over events if >10,000 events
		do
			
			new=$(echo $i-$j | sed 's/^0*//')
			ddump -i -e $j -n 1000 /home/arratialab2/rcdaq-$i-0000.evt > /home/arratialab2/fcal/$DIR/$subDIR/$new.txt			
			echo $new
		done
	done
