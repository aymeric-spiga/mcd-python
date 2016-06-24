#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import cgi, cgitb 
#from modules import *

########################################

# for debugging in web browser
cgitb.enable()

## Create instance of FieldStorage 
#form = cgi.FieldStorage() 

##### NOW WRITE THE HTML PAGE TO USER
print "Content-type:text/html;charset=utf-8\n"
print     #Apache needs a space after content-type
header="""<html><head><title>Mars Climate Database: The Web Interface</title></head><body>"""
print header
print "THIS IS A TEST!"
#print "<img src='"+figname+"'><br />"
bottom = "</body></html>"
print bottom

