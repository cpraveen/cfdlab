gfortran -o tidy tidy.for

# Tidy files
FILES=`ls -1 *.f`

for f in $FILES
do
   echo $f
   ./tidy $f
done

read -p "Press enter to rename files: "

# Rename them
FILES=`ls -1 *.tid`

for f in $FILES
do
   echo "Rename $f"
   base="${f%.*}"
   mv $f $base.f
   rm -f $base.lis
done
