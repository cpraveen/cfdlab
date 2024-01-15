# convert pdf to pptx
# Needs
#   gs
#   soffice (from libreoffice)

FILENAME=$1  # Include .pdf extension
RAND=$RANDOM # Needed to avoid overwriting 
NEW_FILE=/tmp/${RAND}_${FILENAME}
gs -o $NEW_FILE -dNoOutputFonts -sDEVICE=pdfwrite $FILENAME
#gs -o $NEW_FILE -dBATCH -dNoOutputFonts -dPDFSETTINGS=/prepress -dEmbedAllFonts=true -sDEVICE=pdfwrite $FILENAME
#gs -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=$NEW_FILE -c ".setpdfwrite </NeverEmbed []>> setdistillerparams" -f $FILENAME
soffice --infilter=impress_pdf_import --convert-to pptx --outdir /tmp $NEW_FILE
PPTX_NAME_WITH_RAND="${NEW_FILE%.*}"
PPTX_NAME_WITHOUT_RAND="${FILENAME%.*}"
mv $PPTX_NAME_WITH_RAND.pptx $PPTX_NAME_WITHOUT_RAND.pptx
echo "Converted file is ${PPTX_NAME_WITHOUT_RAND}.pptx"
