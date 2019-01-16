#!/bin/bash
# script to enclose all source files in a namespace
#
# NB: this does not always work, i.e. it still failed for PvalueVectors.h 
# because the last include was inside an if statement. However, in most cases 
# it works, so the exceptions can always be fixed by hand...

# NB2: for the files containing a callable main function, we should manually
# remove the namespace scope and then add the namespace to the relevant objects
cd ../src/

for f in MaRaCluster.cpp MaRaCluster.h; do
  if (grep "namspace maracluster" $f);then 
    echo "No need to insert namespace in $f"
  else
    line=$(($(grep -n "#include" $f | cut -f1 -d: | tail -n1) + 1))
    echo $line
    sed "${line}i\\\\nnamespace maracluster {" $f > $f.new
    
    if [[ "$f" == *.cpp ]]
    then
      echo -e '\n} /* namespace maracluster */' >> $f.new
      #mv $f.new $f
    else
      line=$(grep -n "#endif" $f.new | cut -f1 -d: | tail -n1)
      echo $line
      #sed "${line}i} /* namespace maracluster */\n" $f.new > $f
      #rm $f.new
    fi    
    echo "Namespace inserted in $f"
  fi 
done
