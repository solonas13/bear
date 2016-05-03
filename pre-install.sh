#! /bin/sh

unzip libFLASM.zip
cd libFLASM
unzip seqan.zip
make clean
make
cd ..

tar -xvf sdsl-lite.tar.gz
cd sdsl-lite
./install.sh "$(pwd)"/libsdsl
mv libsdsl/ ..
cd ..

tar -xvf libdatrie_0.2.8.orig.tar.gz
cd libdatrie-0.2.8
./configure --prefix="$(pwd)"/libdatrie
make
make install
mv libdatrie ..
cd ..

tar -xvf ahocorasick.tar.gz
cd ahocorasick
make clean
make
cd ..
