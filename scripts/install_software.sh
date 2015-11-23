#!/usr/bin/env bash

mkdir -p src

pushd src

# Screed
echo "Installing screed"
git clone https://github.com/ged-lab/screed.git
pushd screed
git checkout protocols-v0.8.3
python setup.py install
popd



# Khmer
echo "Installing the khmer protocol"
git clone https://github.com/ged-lab/khmer.git
cd khmer
git checkout protocols-v0.8.3
make
make install
popd



# Trimmomatic v0.33
cd /root
curl -O http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.33.zip
unzip Trimmomatic-0.32.zip
cd Trimmomatic-0.30/
cp trimmomatic-0.30.jar
