#!/usr/bin/env bash
set -vexu

install_root=$1

apt-get install -y software-properties-common
apt-add-repository universe
apt-get update

apt-get install -y \
  build-essential \
  python3 \
  python3-pip \
  python3-setuptools \
  wget

pip3 install tox

if [ ! -d $install_root ]; then
  mkdir $install_root
fi

#________________________ quicktree _________________________#
cd $install_root
wget https://github.com/khowe/quicktree/archive/v2.5.tar.gz
tar xf v2.5.tar.gz
cd quicktree-2.5/
make
cd ..
cp -s quicktree-2.5/quicktree .
