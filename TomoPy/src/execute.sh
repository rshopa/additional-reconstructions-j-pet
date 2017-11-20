echo $'\n==================================='
echo $' STARTING RECONSTRUCTION PROCESS'
echo $'===================================\n'
echo $'Current directory: \n'$PWD$'\n'

# input files have got 'PSF-' prefix
CENTERED=$1  # optional parameter (z centered)
# echo $OPTIONS
INPUT_LIST=$(ls -t -U | grep "PSF")
VECTORS=0  # whether vectors of coordinates have been saved

for x in $INPUT_LIST
do
  # validate if vectors are already exported into files
  if [ $VECTORS = 0 ]
  then
    INPUT_STRING=$x
    VECTORS=1   # assign 1 to avoid duplicates further
  else
    INPUT_STRING=$x$' novectors'
  fi

  # validate if z position is non-zero
  if echo $x | grep -q '1875'
  then
    INPUT_STRING=$INPUT_STRING$' biased'
  fi

  # validate if z axis is set as centered (i.e. bash shell.sh centered)
  if [ "$CENTERED" = "centered" ]
  then
    INPUT_STRING=$INPUT_STRING$' centered'
  fi

  echo $'Process started for '$x
  echo $'---------------------------------------------------------------\n'
  python ../get_slices.py $INPUT_STRING
  # echo "python ../get_slices.py "$INPUT_STRING
  echo $'\nProcess ended for '$x
  echo $'===============================================================\n'
done
