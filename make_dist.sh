#!/bin/sh

script=hgmd2VScustomAnno.pl
doc='Instructions for hgmd2VSca.tex'
version=`git describe --abbrev=6 HEAD`

### perl script
cat $script | sed s/\$VERSION\\$/$version/ > exe/$script
chmod 775 exe/$script

mkdir -p dist

rsync -a exe/$script dist/

### documentation
rsync -a doc dist/
cd dist/doc
cat "$doc" | sed s/\$VERSION\\$/$version/ > "$doc.tmp"
mv "$doc.tmp" "$doc"
pdflatex "$doc" > pdflatex.log
cd ../
base=`basename "$doc" .tex`
mv "doc/$base.pdf" .
rm -r doc
cd ../

echo "Now run build_exe.bat in exe/ on a Windows box, followed by exe/dist.sh"
