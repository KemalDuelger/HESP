#!/bin/bash

EXECUTABLE="./stream-omp-host"

value=65536
i=1
while [ "$i" -le  13 ]
do
    echo "Starte mit $value"
    $EXECUTABLE $value
    value=$(($value*2));
    i=$((i+1));
done