#!/bin/bash

PROJECT_NAME="leon3mp"
EN64="--64bit"
READ_WRITE_SETTINGS="--read_settings_files=on --write_settings_files=off"

function get_start_time() {
    if [ -e "$PROJECT_NAME.$1.rpt" ]
    then
        mv "$PROJECT_NAME.$1.rpt" "$PROJECT_NAME.$1.old.rpt"
    fi
    PTIME=0

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

    PROCESSING_TIME=`expr $STOP_TIME_FULL - $START_TIME_FULL`
    PROCESSING_TIME_H=`expr $PROCESSING_TIME / 60`
    PROCESSING_TIME_M=`expr $PROCESSING_TIME % 60`
    echo "processing time: $PROCESSING_TIME_H:$PROCESSING_TIME_M"
}

function timer() {
    if [ -e "$PROJECT_NAME.$1.rpt" ]
    then
        mv "$PROJECT_NAME.$1.rpt" "$PROJECT_NAME.$1.old.rpt"
    fi
    PTIME=0
    while [ ! -e "$PROJECT_NAME.$1.rpt" ]
    do
        echo -ne "\r$PTIME min "
        PTIME=`expr $PTIME + 1`
        sleep 1m
    done
}
#======================================================================================================
date
echo "Start mapping "

# get_start_time map
quartus_map $EN64 --parallel=on $READ_WRITE_SETTINGS $PROJECT_NAME -c $PROJECT_NAME &>longlog.map.log &
timer map
# get_stop_time

ERROR_PRESENT=`grep -c "Error:" "$PROJECT_NAME.map.rpt"`
if [ $ERROR_PRESENT == 0 ]
then
    echo -e "\nStart merging"
    quartus_cdb $EN64 $READ_WRITE_SETTINGS $PROJECT_NAME -c $PROJECT_NAME --merge=on 2>&1 1>>longlog.map.log
else
    echo -e "\nFinished with errors in mapping"
    exit
fi
#======================================================================================================
ERROR_PRESENT=`grep -c "Error:" "$PROJECT_NAME.merge.rpt"`
if [ $ERROR_PRESENT == 0 ]
then
    echo "Start fitting "

#     get_start_time fit
    quartus_fit $EN64 --parallel=4 $READ_WRITE_SETTINGS $PROJECT_NAME -c $PROJECT_NAME &>longlog.fit.log &
    timer fit
#     get_stop_time
else
    echo "Finished with errors in merging"
    exit
fi
#======================================================================================================
ERROR_PRESENT=`grep -c "Error:" "$PROJECT_NAME.fit.rpt"`
if [ $ERROR_PRESENT == 0 ]
then
    echo -e "\nStart STA and assembling"
    quartus_sta $EN64 --parallel=4 $PROJECT_NAME -c $PROJECT_NAME 1>/dev/null
    quartus_asm $EN64 $READ_WRITE_SETTINGS $PROJECT_NAME -c $PROJECT_NAME 1>/dev/null &
    head -n 2 $PROJECT_NAME.sta.rpt
    head -n 161 $PROJECT_NAME.sta.rpt | tail -n 11 | cut -d ";" -f1,2,3,4

    echo "All finished"
else
    echo -e "\nFinished with errors in fitting"
fi
