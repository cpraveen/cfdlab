#!/bin/bash
# Developed by : Pritam Giri
# Date : 8.10.2016
# TIFR-CAM
i=0;
flag=1;
while [ "$flag" = "1" ];
do
	if ls $(printf "sol-%0.3d*" $i) >/dev/null 2>&1;
	then
		cat $(printf "sol-%0.3d*" $i) > $(printf "sol-%0.3d.plt" $i)
		rm $(printf "sol-%0.3d-*" $i)
		i=$((i + 1))
	else
		flag=0
	fi
done
