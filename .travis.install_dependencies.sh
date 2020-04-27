#!/usr/bin/env bash
set -ex
wget https://github.com/khowe/quicktree/archive/v2.5.tar.gz
tar xf v2.5.tar.gz
cd quicktree-2.5/
make
export PATH=$PWD:$PATH
