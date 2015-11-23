#!/usr/bin/env bash

mkdir -p src
mkdir -p data/adapters

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
unzip Trimmomatic-0.33.zip
rm Trimmomatic-0.33.zip
cp Trimmomatic-0.33/trimmomatic-0.33.jar
cp Trimmomatic-0.33/adapters/*.fa ../data/adapters/


# libgtexutils
curl -O http://hannonlab.cshl.edu/fastx_toolkit/libgtextutils-0.6.1.tar.bz2
tar xvjf libgtextutils-0.6.1.tar.bz2
rm libgtextutils-0.6.1.tar.bz2
pushd libgtextutils-0.6.1/
./configure --prefix=$(pwd)/../.. && make
popd


# fastx-toolkit
curl -O http://hannonlab.cshl.edu/fastx_toolkit/fastx_toolkit-0.0.13.2.tar.bz2
tar xjvf fastx_toolkit-0.0.13.2.tar.bz2
rm fastx_toolkit-0.0.13.2.tar.bz2
pushd fastx_toolkit-0.0.13.2/
./configure --prefix=$(pwd)/../.. && make && make install


