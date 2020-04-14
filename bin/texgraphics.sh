#!/bin/bash
#
# Usage
#    texgraphics.sh <tex file>
#    texgraphics *.tex
#    texgraphics.sh paper1.tex paper2.tex ...
#
fnames=`grep includegraphics $@ | sed 's/^[^{]*{\([^{}]*\)}.*/\1/'`
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
