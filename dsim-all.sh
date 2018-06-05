#!/bin/sh
./dsim -gb 0
./dsim -gb 1
./dsim -gb 2
./dsim -gb 3
./dsim -gb 4
wait
wait
wait
wait
wait
cat rle4.txt rle6.txt rle8.txt rle12.txt rle20.txt > rle-all.txt
