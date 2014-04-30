#!/bin/sh

wget  http://fastjet.fr/repo/fastjet-3.0.6.tar.gz
tar -xzf fastjet-3.0.6.tar.gz
mkdir FASTJET

export CXXFLAGS="-m64"
export FFLAGS="-m64"
export CPPFLAGS="-m64"
export CFLAGS="-m64"

export DIR=$PWD/FASTJET
# Compile and install Fastjet.
cd fastjet-3.0.6
./configure --prefix=$DIR
make 
make check
make install

# Come back to the original directory, and clean up.
cd ../
\rm -r fastjet-3.0.6
\rm fastjet-3.0.6.tar.gz

export version=`$PWD/FASTJET/bin/fastjet-config --version`
echo "*******************************************************************"
echo "Fastjet version installed in :$PWD/FASTJET : $version"
echo "*******************************************************************"

#source setup_root_64.csh
