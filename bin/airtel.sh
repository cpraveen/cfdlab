FILES=`ls -tU 14713537_*.pdf`

COMMAND=""

for f in $FILES
do
   COMMAND="$COMMAND $f 1"
done

pdfjam $COMMAND -o out.pdf
