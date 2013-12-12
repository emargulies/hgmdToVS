#!/bin/sh

### first run the bat file on windows.
### second copy that file back here.

if [ ! -e "Convert HGMD to VariantStudio Custom Annotation.exe" ]; then
    echo "you need to run the bat file on windows first. see ../README"
    exit
fi

zip convert.zip "Convert HGMD to VariantStudio Custom Annotation.exe"
mv convert.zip ../dist/
rm "Convert HGMD to VariantStudio Custom Annotation.exe"
rm hgmd2VScustomAnno.pl
