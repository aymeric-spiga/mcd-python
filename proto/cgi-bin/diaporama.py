#! /usr/bin/env python

from os import listdir
from random import choice

ext2conttype = {"jpg": "image/jpeg",
                "jpeg": "image/jpeg",
                "png": "image/png",
                "gif": "image/gif"}

def content_type(filename):
    return ext2conttype[filename[filename.rfind(".")+1:].lower()]

def isimage(filename):
    """true if the filename's extension is in the content-type lookup"""
    filename = filename.lower()
    return filename[filename.rfind(".")+1:] in ext2conttype

def random_file(dir):
    """returns the filename of a randomly chosen image in dir"""
    images = [f for f in listdir(dir) if isimage(f)]
    return choice(images)

if __name__ == "__main__":
    dir = "/home/aymeric/Images/wallpapers/"
    r = random_file(dir)
    print "Content-type: %s\n" % (content_type(r))
    print file(dir+r, "rb").read()

