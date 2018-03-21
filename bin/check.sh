#!/bin/sh



hello=1

mandatory(){
	while [ $# -gt 0 ]; 
	do
		x=$1
		if [ ! -z ${!x} ]; then
			echo $x
		fi
		shift
	done
}

mandatory hello world
