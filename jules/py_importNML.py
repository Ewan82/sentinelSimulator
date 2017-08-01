#!/home/db903833/dataLand01/enthoughtDistros/epd-7.2-2-rh5-x86_64/bin/python

import re
import sys


def importJulesNML(nml):
    """Parse a JULES nml file and write it in the style used for the
    julesNML class.

    :param nml: JULES NML file name.
    :type nml: str.
    :return: None
    """

    variables = {}
    members = {}
    inNmlMember = False

    print "%s=\"\"\"" % nml.replace(".nml", "_txt")
    nmlText = open(nml).readlines()
    for line in nmlText:

        if re.match(r'&', line):
            # this line defines a namelist member but there can be mulitple
            # namelist members of the same name and we need to deal with that

            i = 1
            memName = line.replace('&', '').replace('\n', '') + "_%d" % i
            while True:
                if memName in members:
                    i += 1
                    memName = line.replace('&', '').replace('\n', '') + "_%d" % i
                else:
                    break
            members[memName] = True
            inNmlMember = True
            print line,
            continue

        if re.match(r'/', line):
            # we have left a namelist member
            inNmlMember = False
            print line
            continue

        if re.search(r"\=", line) != None:
            # this line a variable so we must store its values
            varName = memName + "_" + line.split('=')[0]
            variables[varName] = line.split('=')[1].replace('\n', '')
            print "%s=" % line.split('=')[0],
            print "${%s}," % varName
            continue

        if re.search(r"\=", line) == None and inNmlMember == True:
            # this is probably a continuation line so just add
            # it ont to the previous variable.
            # Don't print though! (it comes out later!)
            variables[varName] = variables[varName] + "" + line.replace('\n', '')
            continue

    print "\"\"\""

    print "%s=julesNML(%s,\"%s\")" % (nml.replace(".nml", "_nml"), nml.replace(".nml", "_txt"), nml)

    # set mappings for default variable values:
    for varName in variables:
        print "%s.mapping[\"%s\"]=\"%s\"" % (nml.replace(".nml", "_nml"), varName, variables[varName])


if __name__ == "__main__":

    import os

    if len(sys.argv) > 1:
        os.chdir(sys.argv[1])

    from glob import glob

    NML_files = glob("*nml")

    # write the julesNML class

    print """
#!/usr/bin/python
\"\"\"This module holds JULES namelist files. 
It has been automatically generated.
\"\"\"

import sys
from string import Template

class julesNML:
  \"\"\"
  This is the base class for storing
  and writing JULES namelist files
  \"\"\"

  def __init__( self, template, filename ):
  
    self.t=Template(template)
    self.filename=filename
    self.mapping={}
       
  def update( self, template ): 
    self.t=Template(template)      
        
  def write( self ):
  
    f=open( self.filename, "w" )
    print >> f , self.t.safe_substitute( self.mapping )
    f.close( )
"""

    # parse the individual NML files

    for nml in NML_files:
        importJulesNML(nml)
        print "\n\n\n"

    # write the AllNML class

    print """
class julesAllNML:

  def __init__( self ):
"""

    for nml in NML_files:
        print "    self.%s=%s" % (nml.replace(".nml", "_nml"), nml.replace(".nml", "_nml"))
        print "    self.%s=%s" % (nml.replace(".nml", "_txt"), nml.replace(".nml", "_txt"))

    print "\n"
    print "  def writeNML( self ):"

    for nml in NML_files:
        print "    self.%s.write()" % (nml.replace(".nml", "_nml"))

    print ""
