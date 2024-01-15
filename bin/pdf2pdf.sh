# Converts
#   pdf --> ps --> eps --> pdf 
# and then runs pdfbb, which must be in your path.
# This ensures all fonts/symbols are embedded in the pdf.
#
# NOTE THAT IF YOU HAVE A PS OR EPS FILE OF SAME NAME, IT WILL BE DELETED.

# Find all pdf in this and all subdirectories
PDF=`find . -name "*.pdf"`

for f in $PDF
do
   BASENAME="${f%.*}"
   echo "===> Processing " $BASENAME
   rm -f $BASENAME.ps $BASENAME.eps
   pdf2ps $BASENAME.pdf $BASENAME.ps
   ps2eps -O $BASENAME.ps 
   ps2pdf $BASENAME.eps $BASENAME.pdf
   pdfbb $BASENAME.pdf
   rm -f $BASENAME.ps $BASENAME.eps
done
