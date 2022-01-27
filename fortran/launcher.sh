#!/bin/bash

InputFile=$1
OutputFile=${InputFile%%.inp}.out
ScratchFile=${InputFile%%.inp}.integrals
QuoteFile=$((1 + $RANDOM % 45))

cat ./writer/main_to_write/header > ${OutputFile}

echo $1 > tmp.1
echo `./HartreeFockRoothaan.exe < tmp.1`

cat ./tmp/out.out >> ${OutputFile}

echo "*******UNIQUE OVERLAP MATRIX VALUES*******" >> ${OutputFile}
cat ./tmp/Overlap.int >> ${OutputFile}

echo "*******UNIQUE KINETIC MATRIX VALUES*******" >> ${OutputFile}
cat ./tmp/Kinetic.int >> ${OutputFile}

echo "******UNIQUE POTENTIAL MATRIX VALUES******" >> ${OutputFile}
cat ./tmp/Potential.int >> ${OutputFile}

echo "*****UNIQUE TWO ELECTRON MATRIX VALUES****" >> ${OutputFile}
cat ./tmp/TwoElectron.int >> ${OutputFile}

cat ./writer/succesfull_quotes/${QuoteFile} >> ${OutputFile}

rm tmp.1 ./tmp/out.out ./tmp/Overlap.int ./tmp/Kinetic.int ./tmp/Potential.int ./tmp/TwoElectron.int
