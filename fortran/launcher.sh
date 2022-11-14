#!/bin/bash

InputFile=$1

if $2; then
  OutputFile=$2
else
  OutputFile=${InputFile%%.inp}.out
fi

ScratchFile=${InputFile%%.inp}.integrals
QuoteFile=$(($RANDOM % 45))

cat ./writer/header > ${OutputFile}

./bin/xRHFR.exe $1

if [[ ! -f "./tmp/summary.out" ]]; then

  cat ./tmp/summary.out >> ${OutputFile}

  if ( $2 && $2 == "integrals" ); then

    echo "*******UNIQUE OVERLAP MATRIX VALUES*******" >> ${OutputFile}
    cat ./tmp/Overlap.int >> ${OutputFile}
  
    echo "*******UNIQUE KINETIC MATRIX VALUES*******" >> ${OutputFile}
    cat ./tmp/Kinetic.int >> ${OutputFile}
    
    echo "******UNIQUE POTENTIAL MATRIX VALUES******" >> ${OutputFile}
    cat ./tmp/Potential.int >> ${OutputFile}
    
    echo "*****UNIQUE TWO ELECTRON MATRIX VALUES****" >> ${OutputFile}
    cat ./tmp/TwoElectron.int >> ${OutputFile}
    
  fi

  cat ./writer/succesfull_quotes/${QuoteFile} >> ${OutputFile}

else

  if 
  cat ./tmp/error.out >> ${OutputFile}

fi

rm -f ./tmp/out.out ./tmp/Overlap.int ./tmp/Kinetic.int ./tmp/Potential.int ./tmp/TwoElectron.int
