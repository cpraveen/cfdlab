# convert pdf to odp
# Needs
#   gs
#   soffice (from libreoffice)

FILENAME=$1  # Include .pdf extension
RAND=$RANDOM # Needed to avoid overwriting 
NEW_FILE=/tmp/${RAND}_${FILENAME}
gs -o $NEW_FILE -dNoOutputFonts -sDEVICE=pdfwrite \
   -dPDFSETTINGS=/prepress -dEmbedAllFonts=true \
   $FILENAME
#gs -o $NEW_FILE -dBATCH -dNoOutputFonts -dPDFSETTINGS=/prepress -dEmbedAllFonts=true -sDEVICE=pdfwrite $FILENAME
#gs -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=$NEW_FILE -c ".setpdfwrite </NeverEmbed []>> setdistillerparams" -f $FILENAME
soffice --infilter=impress_pdf_import --convert-to odp --outdir /tmp $NEW_FILE
ODP_NAME_WITH_RAND="${NEW_FILE%.*}"
ODP_NAME_WITHOUT_RAND="${FILENAME%.*}"
mv $ODP_NAME_WITH_RAND.odp $ODP_NAME_WITHOUT_RAND.odp
echo "Converted file is ${ODP_NAME_WITHOUT_RAND}.odp"
