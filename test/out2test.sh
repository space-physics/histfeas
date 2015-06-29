#!/bin/bash
# condenses transcar output
# Michael Hirsch

if (( $# != 3)); then
echo "data_dir test_dir datetime"
exit 1
fi

datdir=$1
testdir=$2
datetime=$3

nhead=5
#nen=170
NprecipCol=2
Ndatacol=11
Nalt=123 #a priori from conttanh.dat

Nprecip=$(($NprecipCol*$Ndatacol))
ndat=$(($Ndatacol*$Nalt))
size_record=$(($ndat+$Nprecip+$nhead))
nrow=$(($size_record/$nhead*2))
nrow=339
echo $nrow

fp="dir.output/emissions.dat"

ddirs=$(find $datdir -maxdepth 1 -mindepth 1 -type d)
for d in ${ddirs[*]}; do
    odir="$testdir/${d##/*/}"
	echo -e "\n$odir"
    mkdir -p $odir/dir.{input,output}
    cp DATCAR "$odir/dir.input/" #not ln, since Cygwin doesn't work

	txt=$(grep -A$nrow "$datetime" $d/$fp)
	echo "$txt" > "$odir/$fp"
done
