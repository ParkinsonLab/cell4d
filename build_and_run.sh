#!/bin/bash
mkdir obj
set -e
if [ "$1" == "clean" ]; then
	make clean
fi
make
./CellApp
