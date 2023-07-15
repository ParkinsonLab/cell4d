#!/bin/bash
cd "${0%/*}"

ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
brew install wxmac --force-bottle

echo "Now please install libSBML from https://sourceforge.net/projects/sbml/files/libsbml/5.15.0/stable/Mac%20OS%20X/"
echo "Press any key to open the webpage"
read -n 1 -s
open "https://sourceforge.net/projects/sbml/files/libsbml/5.15.0/stable/Mac%20OS%20X/"

echo "Press any key to continue after installing libSBML"
read -n 1 -s

echo "`pwd`/simulation &" > ~/Desktop/Cell_Simulation
chmod +x ~/Desktop/Cell_Simulation
./build_and_run.sh clean
