#!/bin/bash -x

# This script is executed on the virtual machine during the Installation phase (need to be ran as root!).
# It is used to record a predefined VM-image of the appliance.
# Otherwise executed first during a cloud deployement in IFB-Biosphere

# Setting up the tutorial environment

APP_SRC="/home/ubuntu/opt/src"
BIN="/home/ubuntu/opt/bin"
mkdir -p "${APP_SRC}"
mkdir -p "${BIN}"

cd "${APP_SRC}" || exit

git clone https://github.com/metaDMG-dev/metaDMG-cpp.git

cd metaDMG-cpp || exit

make clean && make -j 4

mv metaDMG-cpp "${BIN}"

rm -rf "${APP_SRC}/metaDMG-cpp"
