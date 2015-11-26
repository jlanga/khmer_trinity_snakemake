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
pushd khmer
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



# Trinity
curl -L -O https://github.com/trinityrnaseq/trinityrnaseq/archive/v2.1.1.tar.gz
tar xzvf v2.1.1.tar.gz
rm v2.1.1.tar.gz
pushd trinityrnaseq-2.1.1
make -j 24
make -j 24 plugins
make -j 24 test_all
ln -s ../src/trinityrnaseq-2.1.1/Trinity ../../bin/Trinity
popd



# Bowtie
curl -O -L http://sourceforge.net/projects/bowtie-bio/files/bowtie/0.12.7/bowtie-0.12.7-linux-x86_64.zip
unzip bowtie-0.12.7-linux-x86_64.zip
rm bowtie-0.12.7-linux-x86_64.zip
pushd bowtie-0.12.7
cp bowtie bowtie-build bowtie-inspect ../../bin
popd



# samtools
curl -L http://sourceforge.net/projects/samtools/files/latest/download?source=files >samtools.tar.bz2
tar xvjf samtools.tar.bz2
rm samtools.tar.bz2
mv samtools-* samtools-latest
pushd samtools-latest/
make -j 24
cp samtools ../../bin
popd
