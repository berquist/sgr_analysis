#!/usr/bin/env python3

from __future__ import print_function
from glob import glob
import os
from mbe.utils import pad_left_zeros

def rename_snapshot_files():
    filenames = glob('*_snap.*')
    print(len(filenames))

    maxlen = 0

    for oldfilename in filenames:
        newlen = len(oldfilename.split('_')[0])
        if newlen > maxlen:
            maxlen = newlen

    newfilenames = []

    for oldfilename in filenames:
        ext = os.path.splitext(oldfilename)[1][1:]
        splitname = oldfilename.split('_')
        filenumlen = len(splitname[0])
        if filenumlen < maxlen:
            splitname[0] = pad_left_zeros(splitname[0], maxlen)
        newfilename = 'drop_{}.{}'.format(splitname[0], ext)
        os.rename(oldfilename, newfilename)
        print(oldfilename + ' -> ' + newfilename)
        newfilenames.append(newfilename)

    return newfilenames

if __name__ == '__main__':
    rename_snapshot_files()
