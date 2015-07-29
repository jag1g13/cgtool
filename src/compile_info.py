#!/usr/bin/python

from __future__ import print_function
import sys
import os
import pwd
import datetime
import socket

try:
    compilefile = sys.argv[1]
    versionfile = sys.argv[2]
except:
    print("USAGE: compile_info.py info_file version_file")
    sys.exit(0)

# open the include file
try:
    idfile = open(compilefile, "w")
    verfile = open(versionfile, "w")
except:
    print("Cannot open info files for writing")
    sys.exit(0)

# now write down the date that this code was compiled, and who has compiled it
user = pwd.getpwuid(os.getuid())[0]
today = datetime.datetime.today()
today = today.replace(microsecond=0)
host = socket.gethostname()

print("\"Compiled on {0} by {1} at {2}\\n\"".format(host, user, today), file=idfile)

try:
    identify = os.popen("hg identify", "r").readline()
    identify = identify.strip("\n ")
    version = os.popen("hg identify --num", "r").readline()
    version = version.strip("\n ")
    path = os.popen("hg paths default", "r").readline()
    path = path.strip("\n ")
except:
    print("\"Mercurial information not available\\n\"", file=idfile)
else:
    if(identify.startswith("abort") or path == "not found!"):
        print("\"Mercurial information not available\\n\"", file=idfile)
    else:
        print("\"Revision: {0} {1}\\n\"".format(version, identify), file=idfile)
        print("\"From: {0}\\n\"".format(path[path.find("@")+1:]), file=idfile)
        print("\".{0}\"".format(version), file=verfile)
