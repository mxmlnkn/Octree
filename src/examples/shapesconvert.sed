#!/bin/sed -r

# Load all lines into pattern space !
:notlastline
N
# don't really modify line, but give 's' a successful find to make t work
$! s/^//
t notlastline

s/\[n, j, np\]/\(n,j,p\)/g

# Replace := with {
s/[[:blank:]]*:=[[:blank:]]*\n[ \t]*//g
# jump to label a if last substitute command modified pattern space

# delete unnecessary Module command
s/Module\[\{[^{]*\},[[:blank:]]*[\n]*//g

# Replace ; with newline if necessary
s/\;[[:blank:]]*[\n]*/\;\n/g

# replace (...)^2 with pow(...,2)
s/\(([^(]*)\)\^2\)/pow\(\1,2\)/g

# replace var^2 with pow(var,2)
s/([[:alpha:]]{1,})\^2/pow\(\1,2\)/g

s/\[n_, j_, np_\]/\( int n, int j, int np \) {\n/g


# make ...["..."] a valid variable
s/([[:alnum:]]+)\[\"([[:alnum:]]+)\-([[:alnum:]]+)\"\]/\1\2m\3/g
s/([[:alnum:]]+)\[\"([[:alnum:]]+)\+([[:alnum:]]+)\"\]/\1\2p\3/g
s/([[:alnum:]]+)\[\"([[:alnum:]]+)\"\]/\1\2/g

# remove leading whitespace
s/\n^[[:blank:]]*/\n/g
# remove wrongly added 'int' prefix if not at beginning of line do this until no more ints found
#:intfound
#s/\n[[:blank:]]*([^\n]+)int /\n\1/g
#t intfound

s/Ceiling\[\j\/2\]/ceil( (double)j \/ 2.0 )/g

# Remove Mathematica Comments
s/\(\*[[:blank:]]*([^*]*)[[:blank:]]*\*\)/\/\* \1 \*\//g

:noheader

