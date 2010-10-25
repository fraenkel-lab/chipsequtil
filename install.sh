#!/bin/bash

# this script installs chipsequtils into /usr/local on the cluster nodes
# since the cluster nodes do not have root write access to the network
# volumes from nodes other than node 9, sudo ... doesn't work with
# setuptools because it writes egg_info to the source directory and I
# can't figure out how to get it to write to a local directory.
#
# this script copies the entire chipsequtil source tree to /tmp on the
# local machine, runs sudo ./setup.py install --prefix=/usr/local, and,
# on success, deletes the temporary source directory
#
# it _must_ be run from the source directory

TMPDIR="/tmp/chipsequtil_tmp_$(date +%F)"
if [ ! -d $TMPDIR ]; then
    echo "temporary source dir $TMPDIR does not exist, creating"
    mkdir $TMPDIR
fi

cd ../
echo "copying source tree to $TMPDIR"
cp -vr -t $TMPDIR chipsequtil/{setup.*,ez_setup.py,src,scripts,setuptools*}
cd $TMPDIR
echo "cd'ed to $PWD, installing"
sudo ./setup.py install --prefix=/usr/local
if [ $? -eq 0 ]; then
        echo "install successful, removing $TMPDIR"
        cd
        sudo rm -r $TMPDIR
fi
