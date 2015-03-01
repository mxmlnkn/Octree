#!/bin/bash
#cat test.txt | sed 's/\[n_, j_, np_\]/\( int n, int j, int np \) {\n/g' | sed 's/\[n, j, np\]/\(n,j,p\)/g' | sed -e :a -e '/[[:space:]]*:=[[:space:]]*/N; s/[[:space:]]*:=[[:space:]]*\n/   return/; ta' | sed -e 's/;/\n}\n/g'

cat test.txt | sed -r -f convert.sed | sed -r '/( int n, int j, int np )/! s/^[[:blank:]]*/    /' > converted.txt

