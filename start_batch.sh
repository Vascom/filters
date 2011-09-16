#!/bin/bash

PROJECT_NAME="leon3mp"
EN64="--64bit"
READ_WRITE_SETTINGS="--read_settings_files=on --write_settings_files=off"

function get_start_time() {
    START_TIME=`date | cut -d " " -f4 | cut -d ":" -f1`
    START_TIME_H=`expr $START_TIME \* 60`
    START_TIME_M=`date | cut -d " " -f4 | cut -d ":" -f2`
    START_TIME_FULL=`expr $START_TIME_H + $START_TIME_M`
}

function get_stop_time() {
    STOP_TIME=`date | cut -d " " -f4 | cut -d ":" -f1`
    STOP_TIME_H=`expr $STOP_TIME \* 60`
    STOP_TIME_M=`date | cut -d " " -f4 | cut -d ":" -f2`
    STOP_TIME_FULL=`expr $STOP_TIME_H + $STOP_TIME_M`
}

function print_time() {
    PROCESSING_TIME=`expr $STOP_TIME_FULL - $START_TIME_FULL`
    PROCESSING_TIME_H=`expr $PROCESSING_TIME / 60`
    PROCESSING_TIME_M=`expr $PROCESSING_TIME % 60`
    echo "processing time: $PROCESSING_TIME_H:$PROCESSING_TIME_M"
}
#======================================================================================================
date
echo -n "Start mapping "

get_start_time
quartus_map $EN64 --parallel=on $READ_WRITE_SETTINGS $PROJECT_NAME -c $PROJECT_NAME &>longlog.map.log
get_stop_time
print_time

ERROR_PRESENT=`grep -c "Error:" "$PROJECT_NAME.map.rpt"`
if [ $ERROR_PRESENT == 0 ]
then
    echo "Start merging"
    quartus_cdb $EN64 $READ_WRITE_SETTINGS $PROJECT_NAME -c $PROJECT_NAME --merge=on 2>&1 1>>longlog.map.log
else
    echo "Finished with errors in mapping"
    exit
fi
#======================================================================================================
ERROR_PRESENT=`grep -c "Error:" "$PROJECT_NAME.merge.rpt"`
if [ $ERROR_PRESENT == 0 ]
then
    echo -n "Start fitting "

    get_start_time
    quartus_fit $EN64 --parallel=4 $READ_WRITE_SETTINGS $PROJECT_NAME -c $PROJECT_NAME &>longlog.fit.log
    get_stop_time
    print_time
else
    echo "Finished with errors in merging"
    exit
fi
#======================================================================================================
ERROR_PRESENT=`grep -c "Error:" "$PROJECT_NAME.fit.rpt"`
if [ $ERROR_PRESENT == 0 ]
then
    echo "Start STA and assembling"
    quartus_sta $EN64 --parallel=4 $PROJECT_NAME -c $PROJECT_NAME 1>/dev/null
    quartus_asm $EN64 $READ_WRITE_SETTINGS $PROJECT_NAME -c $PROJECT_NAME 1>/dev/null &
    head -n 2 $PROJECT_NAME.sta.rpt
    head -n 159 $PROJECT_NAME.sta.rpt | tail -n 9

    echo "All finished"
else
    echo "Finished with errors in fitting"
fi
