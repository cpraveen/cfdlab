#!/bin/sh
# Convert ps or eps files to pdf
# 06 August, 2009

PN=`basename "$0"`

usage () {
   echo >&2 "usage: $PN <eps file>

   Convert ps or eps files to pdf
   example: $PN foo.eps
            $PN foo1.eps foo2.eps
            To convert all ps files
            $PN \"*.ps\"

   WARNING: THIS PROGRAM OVERWRITES THE PDF FILE !!!"

    exit 1
}

if [ $# -lt 1 ]
then
	usage
fi

for file in "$@"
do
   if [ -f $file ]
   then
      echo "Converting $file"
      ps2pdf $file
   else
      echo "File $file does not exist !"
   fi
done
