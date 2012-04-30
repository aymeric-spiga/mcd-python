#!/usr/bin/python3
# -*- coding: UTF-8 -*-
'''
Created on 30 mars 2012
 
@author: Nicolas Vergnes
'''
 
from http.server import HTTPServer, CGIHTTPRequestHandler
 
# Port du serveur http
port = 8080
# Allocation de l'objet de gestion du serveur
httpd = HTTPServer(('', port), CGIHTTPRequestHandler)
 
# On lance le serveur indéfiniement
print('Lancement du serveur http sur le port: ' + str(httpd.server_port))
httpd.serve_forever()
