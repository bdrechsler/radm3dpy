"""
 RADMC-3D Python module
 (c) Attila Juhasz, Leiden, 2011,2012,2013

 Generate source files to set the pythonpath variable
"""
import os

wfile=open("sourceme.sh", 'w')
wfile.write("if [ -z ${PYTHONPATH} ]\n")
wfile.write("    then\n")
wfile.write("         export PYTHONPATH="+os.getcwd()+"\n")
wfile.write("    else\n")
wfile.write("         export PYTHONPATH=${PYTHONPATH}:"+os.getcwd()+"\n")
wfile.write("fi\n")
wfile.close()

wfile=open("sourceme.csh", 'w')
wfile.write("if ! $?PYTHONPATH then\n")
wfile.write("     setenv PYTHONPATH "+os.getcwd()+"\n")
wfile.write("else\n")
wfile.write("     setenv PYTHONPATH ${PYTHONPATH}:"+os.getcwd()+"\n")
wfile.write("endif\n")

wfile.close()


