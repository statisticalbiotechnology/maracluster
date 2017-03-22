#!/bin/sh
# script to copy the headers to all the source files and header files

cd ../src/

for f in *.cpp *.h; do
  if (grep Copyright $f);then 
    echo "No need to copy the License Header to $f"
  else
    cat license.txt $f > $f.new
    mv $f.new $f
    echo "License Header copied to $f"
  fi 
done
