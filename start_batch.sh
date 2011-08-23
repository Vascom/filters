#!/bin/bash

PROJECT_NAME="leon3mp"
EN64="--64bit"
READ_WRITE_SETTINGS="--read_settings_files=on --write_settings_files=off"

echo "Start mapping"
quartus_map $EN64 --parallel=on $READ_WRITE_SETTINGS $PROJECT_NAME -c $PROJECT_NAME &>longlog.map.log

ERROR_PRESENT=`cat $PROJECT_NAME.map.rpt | grep -c Error:`
if [ $ERROR_PRESENT == 0 ]
then
    echo "Start merging"
    quartus_cdb $EN64 $READ_WRITE_SETTINGS $PROJECT_NAME -c $PROJECT_NAME --merge=on 2>&1 1>>longlog.map.log
else
    echo "Finished with errors in mapping"
    exit
fi

ERROR_PRESENT=`cat $PROJECT_NAME.merge.rpt | grep -c Error:`
if [ $ERROR_PRESENT == 0 ]
then
    echo "Start fitting"
    quartus_fit $EN64 --parallel=4 $READ_WRITE_SETTINGS $PROJECT_NAME -c $PROJECT_NAME &>longlog.fit.log
else
    echo "Finished with errors in merging"
    exit
fi

ERROR_PRESENT=`cat $PROJECT_NAME.fit.rpt | grep -c Error:`
if [ $ERROR_PRESENT == 0 ]
then
    echo "Start STA and assembling"
    quartus_sta $EN64 --parallel=4 $PROJECT_NAME -c $PROJECT_NAME 1>/dev/null
    quartus_asm $EN64 $READ_WRITE_SETTINGS $PROJECT_NAME -c $PROJECT_NAME 1>/dev/null &
    cat $PROJECT_NAME.sta.rpt | head -n 2
    cat $PROJECT_NAME.sta.rpt | head -n 159 | tail -n 9

    echo "All finished"
else
    echo "Finished with errors in fitting"
fi

