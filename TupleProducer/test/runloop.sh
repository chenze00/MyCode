#!/bin/sh

FILEPATH="/eos/cms/store/group/phys_tau/TauML/prod_2018_v2/full_tuples/DYJetsToLL_M-50"

for FILE in $FILEPATH/*.root
do
    NAME=${FILE##*/}
    NAME=${NAME%.*}
    echo $FILE $NAME >> jobs/items.txt
#    break
done
