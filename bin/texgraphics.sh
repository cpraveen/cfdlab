#!/bin/bash
#
# copygraphics.sh texfile directory
#
fnames=`grep includegraphics $1 | sed 's/^[^{]*{\([^{}]*\)}.*/\1/'`
files=
for f in $fnames
do
   echo ">>> Copying $f"
   files+=" "$f.*
done
echo ">>> All files"
echo $files
echo ">>> Creating archive"
tar zcvf graphics.tgz $files
echo ">>> See graphics.tgz"
