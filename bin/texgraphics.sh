#!/bin/bash
#
# Usage
#    texgraphics.sh <tex file>
#    texgraphics *.tex
#    texgraphics.sh paper1.tex paper2.tex ...
#
# Do not use tab/indenting in tex files.
# Do not put graphics filename extension, use it like this
#
#    \includegraphics[width=0.5\textwidth]{foo}
#
# Ignores commented includegraphics line, but will not work if comment % is at 
# end of line, e.g., this line will be ignored also
#
#    \includegraphics[width=0.5\textwidth]{foo}  % this is a comment
#
fnames=`grep includegraphics $@ | grep -v % | sed 's/^[^{]*{\([^{}]*\)}.*/\1/'`
files=
for f in $fnames
do
   # Check file actually exists
   if [ -f $f.ps ] || [ -f $f.eps ] || [ -f $f.pdf ] || [ -f $f.png ]; then
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
