#!/bin/bash
#
# copygraphics.sh texfile directory
#
fnames=`grep includegraphics $1 | sed 's/^[^{]*{\([^{}]*\)}.*/\1/'`
for f in $fnames
do
   echo "Copying $f"
   cp $f.* $2/
done
