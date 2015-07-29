#!/usr/bin/python

from __future__ import print_function
import sys
import os
import pwd
import datetime
import socket

try:
    filename = sys.argv[1]
except:
    print("USAGE: compile_info.py outputfile")
    sys.exit(0)

# open the include file
try:
    idfile = open(filename, "w")
except:
    print("Cannot open info file %s for writing" % idfile)
    sys.exit(0)

# now write down the date that this code was compiled, and who has compiled it
user = pwd.getpwuid(os.getuid())[0]
today = datetime.datetime.today()
today = today.replace(microsecond=0)
host = socket.gethostname()

print("\"Compiled on %s by %s at %s\\n\"" % (host, user, today), file=idfile)

try:
    identify = os.popen("hg identify", "r").readline()
    identify = identify.strip("\n ")
    path = os.popen("hg paths default", "r").readline()
    path = path.strip("\n ")
except:
    print("\"Mercurial information not available\\n\"", file=idfile)
else:
    if(identify.startswith("abort") or path == "not found!"):
        print("\"Mercurial information not available\\n\"", file=idfile)
    else:
        print("\"Revision: %s\\n\"" % identify, file=idfile)
        print("\"From: %s\\n\"" % path[path.find("@")+1:], file=idfile)
