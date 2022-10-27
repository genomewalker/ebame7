#!/bin/bash -x

# This script is executed on the virtual machine during the Installation phase (need to be ran as root!).
# It is used to record a predefined VM-image of the appliance.
# Otherwise executed first during a cloud deployement in IFB-Biosphere

source /etc/profile

apt-get update
apt-get install -y nginx

# Setting up the tutorial environment

mamba env create -f https://raw.githubusercontent.com/genomewalker/ebame7/main/ebame7-aDNA-tutorial/environment.yaml

conda activate ebame7-aDNA-tutorial

APP_SRC="${IFB_MAIN}"/src
BIN="${IFB_MAIN}"/bin
mkdir -p "${APP_SRC}"

cd "${APP_SRC}" || exit

git clone https://github.com/metaDMG-dev/metaDMG-cpp.git

cd metaDMG-cpp || exit

make clean && make -j 4

mv metaDMG-cpp "${BIN}"
