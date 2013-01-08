#!/usr/bin/python

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
    dir = "../img/"
    r = random_file(dir)
    #print "Content-type: %s\n" % (content_type(r))

    print "Content-type:text/html\n"
    print "  "  #Apache needs a space after content-type
    header="""<html><head><title>Mars Climate Database: The Web Interface</title></head><body>"""
    print header

    print "<form name='diapo' action='./diapo.py' method='post'>"
    print "<center>"
    print "<b style='font-size: 125%;'>Mars Climate Database: The Web Interface.</b><br />"
    print "<b style='font-size: 125%;'>Demo mode!</b>"
    print "<input type='submit' value='Click here for another random example' style='font-weight:bold'/><br />"
    print "<a href='../index.html'>back to the main interface.</a><br />"
    print "</center>"
    print "<hr />"

    #dafile = file(dir+r, "rb").read()
    dafile = dir+r
    print "<img src='"+dafile+"'><br />"

    bottom = "</body></html>"
    print bottom
