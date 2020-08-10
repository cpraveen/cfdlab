#!/bin/bash
# Developed by : Pritam Giri
# Date : 8.10.2016
# TIFR-CAM
# Modified by Biswarup Biswas, IIT Delhi, creates single file for each timestep
# To use, this script:
#   Use DM_STENCIL_STAR
#   When saving solution of each partition, use 
#        iend = ibeg+nlocx
#        jend = jbeg+nlocy

# Compute total nx,ny using first time step
inputfile=$(printf "sol-%0.3d*" 0)
for filename in $inputfile              # loop over subdomains
do
	# awk will look at 3rd line and find the 6th element which is the
   # value of nxloc
	nxloc=$(awk -F"[=,]" 'NR==3{print $6}' $filename)
	# awk will look at 3rd line and find the 8th element which is the 
   # value of nyloc
	nyloc=$(awk -F"[=,]" 'NR==3{print $8}' $filename)
	# we are adding all nx and ny defined in the subdomains
	nx=$((nx + nxloc));	ny=$((ny + nyloc))
done
# global nx and ny. We have actually added nx and ny two times, so we 
# need to divide by 2.
nx=$((nx/2)); ny=$((ny/2));
		
i=0;
flag=1;
# flag=1 means we have files to merge
while [ "$flag" = "1" ];
do
	# if some files found in the directory of the pattern "sol-%0.3d*", i
	# i varies over timesteps
	if ls $(printf "sol-%0.3d*" $i) >/dev/null 2>&1;
	then
		# Now we are at a specific time say, t=tstar. We will merge all the 
      # subdomains corresponding to t=tstar into one domain.
		# this collect all the files saved at t=tstar
		inputfile=$(printf "sol-%0.3d*" $i)
		# we will merge in the file name outputfile
		outputfile=$(printf "sol-%0.3d.plt" $i)
		for filename in $inputfile              # loop over subdomains
		do
			# Starting from 4th line copy everything and paste to patch. As the 
         # first 3 lines are header
			awk 'NR>=4' $filename >> patch
		done
		
		# awk will look at the 3rd line and get the defined time tstar.
		soltime=$(awk -F"[=,]" 'NR==3{print $4}' $filename)
		# writing the output file
		# upto 2nd line copy everything and paste to outputfile.
		awk 'NR<3' $filename > $outputfile
		# Now paste the 3rd line with new global nx and ny
		echo "ZONE STRANDID=1, SOLUTIONTIME=$soltime, I=$nx, J=$ny, DATAPACKING=POINT" >> $outputfile
		# sorting is important for plotting in visit
		# finally from 4th line paste the numeric data after sorting over Y.
		sort -g -k2 -k1 patch >> $outputfile
		# end writing
		# remove all subdomains at time tstar and also the patch file 
		rm $(printf "sol-%0.3d-*" $i) patch
		i=$((i + 1)) # go for next timestep
	else
		flag=0
	fi
done
