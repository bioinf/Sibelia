#***************************************************************************
# Copyright (c) 2012 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
#****************************************************************************

#Text resource compiler for Sibelia. Compiler version 0.666

escape = dict()
escape['\\'] = '\\\\'
escape['"'] = '\\"'
maxBlockSize = 1024
namespace = 'namespace SyntenyFinder'
resourceCpp = open('../resource.cpp', 'w')
resourceHeader = open('../resource.h', 'w')

def wrapLine(line):
	token = [escape[ch] if ch in escape else ch for ch in line]
	return 'std::string("' + ''.join(token) + '")'
	
def compileResource(textFileName, varName):
	block = []
	nowBlock = []
	for line in open(textFileName):
		wrapped = wrapLine(line.rstrip())
		if len(nowBlock) == maxBlockSize:		
			block.append(nowBlock)
			nowBlock = []
		nowBlock.append(wrapped)

	if nowBlock != []:
		block.append(nowBlock)
	blockName = []
	for i in xrange(len(block)):
		nowBlock = block[i]
		nowBlockName = varName + 'BlockNo' + str(i)
		blockName.append(nowBlockName)
		print >> resourceCpp, '\t\tconst std::string ' + nowBlockName + '[] = '
		print >> resourceCpp, '\t\t{'
		for line in nowBlock:
			print >> resourceCpp, '\t\t\t' + line + ','
		print >> resourceCpp, '\t\t};\n'
	
	functionName = 'GlueResourceV' + varName
	print >> resourceCpp, '\t\tstd::string ' + functionName + '()\n\t\t{'
	print >> resourceCpp, '\t\t\tstd::stringstream buf;'
	for i in xrange(len(block)):
		name = blockName[i]
		arg = name + ', ' + name + ' + ' + str(len(block[i])) + ', std::ostream_iterator<std::string>(buf, "\\n")'
		print >> resourceCpp, '\t\t\tstd::copy(' + arg + ');'
	print >> resourceCpp, '\t\t\treturn buf.str();'
	print >> resourceCpp, '\t\t}\n'
	return functionName	


print >> resourceHeader, '//****************************************************************************'
print >> resourceHeader, '//* Copyright (c) 2012 Saint-Petersburg Academic University                   '
print >> resourceHeader, '//* All Rights Reserved                                                       '
print >> resourceHeader, '//* See file LICENSE for details.                                             '
print >> resourceHeader, '//****************************************************************************'
print >> resourceHeader, '#include "common.h"\n'

resList = []
print >> resourceCpp, '#include "resource.h"\n'
print >> resourceCpp, namespace + '\n{'
print >> resourceCpp, '\tnamespace\n\t{'
print >> resourceHeader, namespace + '\n{'

for line in open('resource.cfg'):
	line = line.strip()
	if len(line) > 0 and line[0] != '#':
		line = line.split()
		resList.append((line[1], compileResource(line[0], line[1])))
print >> resourceCpp, '\t}\n'

for res in resList:
	varDef = 'const std::string ' + res[0]
	print >> resourceCpp, '\t' + varDef + ' = ' + res[1] + '();'
	print >> resourceHeader, '\textern ' + varDef + ';'

print >> resourceCpp, '}\n'
print >> resourceHeader, '}\n'