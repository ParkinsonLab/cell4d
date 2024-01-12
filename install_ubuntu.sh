#!/bin/bash
cd "${0%/*}"

sudo apt-get -q -y -o Dpkg::Use-Pty=0 install g++ freeglut3-dev libsbml5-dev libwxgtk3.0-gtk3-dev

