#!/bin/bash
#
# Usage
#    texgraphics.sh <tex file>
#    texgraphics *.tex
#    texgraphics.sh paper1.tex paper2.tex ...
#
# Ignore comments, but will not work if comment % is at end of line
fnames=`grep includegraphics $@ | grep -v % | sed 's/^[^{]*{\([^{}]*\)}.*/\1/'`
files=
for f in $fnames
do
   # Check file actually exists
   if [ -f $f.ps ] || [ -f $f.eps ] || [ -f $f.pdf ]; then
      echo ">>> Copying $f"
      files+=" "$f.*
   else
      echo ">>> No file called: $f"
   fi
done
echo ">>> All files"
echo $files
echo ">>> Creating archive"
tar zcvf graphics.tgz $files
echo ">>> See graphics.tgz"
