#!/usr/bin/env python2
#encoding: ISO-8859-15

# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.

#This will take a gro and top file, and take the point charges from there.
#THIS IS NOT YET IMPLEMENTED#This list of charges and coordinates will then be trimmed by a set of point charges which then will be concatenated to a charge shift model, like in a QM/MM approach.
#THIS IS NOT YET IMPLEMENTED#The trimming is optional.
#will print output list as x y z charge in Angstrom
#THIS IS NOT YET TRUE#needs input pdb, file containing groups/atoms to trim (format: <number of groups>\n<for each group:><amount of atoms in group>\n<indices of atoms, one per line>
#needs input .gro file, input .top file, PCF output file

__author__ = "jangoetze"
__date__ = "$15-May-2018 17:02:17$" #during a rain storm

def checkformol(molname,inp):
	import re
	with open(inp) as ifile:
		correct=False
		for line in ifile:
			match=re.search(r'^;', line,flags=re.MULTILINE)
			if match:
				continue
			match=re.search(r'^\[\s*moleculetype\s*\]', line,flags=re.MULTILINE)
			if match:
				for line in ifile:
					match=re.search(r'^;', line,flags=re.MULTILINE)
					if match:
						continue
					else:
						matchstring=r"^\s*" + re.escape(molname)
						match=re.search(matchstring, line,flags=re.MULTILINE)
						if match:
							correct=True
						break
			if correct:
				break
	return correct

def getincludelist(inp):
	import re
	import os.path
	gmxpath="/usr/local/sw/gromacs/share/gromacs/top/"
	toplist=[]
	with open(inp) as ifile:
		for line in ifile:
			match=re.search(r'^;', line,flags=re.MULTILINE)
			if match:
				continue
			match=re.search(r'^#include\s+\"(\S+)\"', line,flags=re.MULTILINE)
			if match:
				match2=re.search('ffbonded',match.group(1))
				if match2:
					continue
				match2=re.search('ffnonbonded',match.group(1))
				if match2:
					continue
				match2=re.search('forcefield.itp',match.group(1))
				if match2:
					continue
				match2=re.search('posre.itp',match.group(1))
				if match2:
					continue
				foundname=match.group(1)
				check=os.path.isfile(foundname)
				if not check:
					foundname=gmxpath+foundname
					check=os.path.isfile(foundname)
					if not check:
						print "File " + foundname + " was not found. Maybe update the gmxpath variable in the script? Exiting."
						exit(1)
				toplist.append(foundname)
				toplist.extend(getincludelist(foundname))
	return toplist	

def readcharges(molvecentry,top):
	import re
	cvec=[]
	curr_top=top
	molname=molvecentry[0]
	molcount=molvecentry[1]
	found=checkformol(molname,top)
	#if found:
		#print "Found " + str(molname) + " in " + str(top) + "."
	if not found:
		toplist=getincludelist(top)
		for element in toplist:
			found=checkformol(molname,element)
			if found:
				#print "Found " + str(molname) + " in " + str(element) + "."
				curr_top=element
				break
	if not found:
		print "No charges found for " + str(molname) + ". Exiting."
		exit(1)
	with open(curr_top) as ifile:
		for line in ifile:
			match=re.search(r'^;', line,flags=re.MULTILINE)
			if match:
				continue
			match=re.search(r'^\[\s*moleculetype\s*\]', line,flags=re.MULTILINE)
			if match:
				for line in ifile:
					match=re.search(r'^;', line,flags=re.MULTILINE)
					if match:
						continue
					matchstring=r"^\s*" + re.escape(molname)
					match=re.search(matchstring, line,flags=re.MULTILINE)
					if match:
						found=True
						#print "Found " + str(molname) + " in " +str(curr_top)
						for line in ifile:
							match=re.search(r'^\[\s*atoms\s*\]', line,flags=re.MULTILINE)
							if match:
								break
						break
					else:
						found=False
						break
				if found:
					break
		for line in ifile:
			match=re.search(r'^\[', line,flags=re.MULTILINE)
			if match:
				break
			match=re.search(r'^;', line,flags=re.MULTILINE)
			if match:
				continue
			match=re.search(r'^\s*\d+\s+\S+\s+\d+\s+\S+\s+\S+\s+\d+\s+([-]*\d+[\.]*[\d+]*)*', line,flags=re.MULTILINE)
			if match:
				cvec.append(float(match.group(1)))
#			else:
#				print "Did not find topology data for " + str(molname) + " in " + str(curr_top) + ". Line was:"
#				print line
	finalcvec=[]
	for i in range(0,int(molcount)):
		finalcvec.extend(cvec)
	return finalcvec

def readg96(inp):
        import re
        coords=[]
        with open(inp) as ifile:
                count=0
                for line in ifile:
                        count+=1
                        if count==4:
                            break
                count=1
                for line in ifile:
                        match=re.search(r'^(.{5})\s(.{5})\s(.{5})\s(.{6})\s*([-]*\d+\.*\d*)\s*([-]*\d+\.*\d*)\s*([-]*\d+\.*\d*)', line,flags=re.MULTILINE)
                        if match:
                                coords.append(float(match.group(5))*10.)
                                coords.append(float(match.group(6))*10.)
                                coords.append(float(match.group(7))*10.)
                        else:
                                break
        return coords

def readgeo(inp):
	import re
	coords=[]
	n_a=0
	with open(inp) as ifile:
		for line in ifile:
			break
		for line in ifile:
			match=re.search(r'^\s*(\d+)', line,flags=re.MULTILINE)
			if match:
				n_a=int(match.group(1))
				break
			else:
				print ".gro is corrupt (no number of atoms found, second line). Exiting."
				exit(1)
		count=1
		for line in ifile:
                        match=re.search(r'^(.{5})(.{5})(.{5})(.{5})\s*([-]*\d+\.*\d*)\s*([-]*\d+\.*\d*)\s*([-]*\d+\.*\d*)', line,flags=re.MULTILINE)
			if match:
				coords.append(float(match.group(5))*10.)
				coords.append(float(match.group(6))*10.)
				coords.append(float(match.group(7))*10.)
			else:
				print ".gro is corrupt. Exiting."
				print "Last line:"
				print line
				exit(1)
			count+=1
			if count>n_a:
				break
	return coords

def readmols(top):
	import re
	mollist=[]
	with open(top) as ifile:
		found=False
		for line in ifile:
			match=re.search(r'^\[ molecules \]', line,flags=re.MULTILINE)
			if match:
				found=True
				break
		if not found:
			print "No \"molecules\" entry in " + str(top) + " found. Exiting."
			exit(1)
		for line in ifile:
			match=re.search(r'^;', line,flags=re.MULTILINE)
			if match:
				continue
			else:
				match=re.search(r'^(\S+)\s+(\d+)', line,flags=re.MULTILINE)
				if match:
					mollist.append([match.group(1),match.group(2)])
				else:
					print "Found an incomprehensible line in molecules list. Exiting."
					print "Last line was:"
					print line
					exit(1)
	return mollist
				
def makeout(coords,charges,name):
	ofile=open(name,"w")
	for i in range(0,len(charges)):
		#ofile.write("{:>8} ".format(str(i+1)))
		for j in range(0,3):
			ofile.write("{:>16.8f} ".format(coords[i*3+j]))
		ofile.write("{:>16.8f}\n".format(charges[i]))
	ofile.close()

def generate_pcf_from_top(gro,top,out):
	import re
	#from numpy import reshape, array
	chargevec=[]
	mollist=readmols(top)
	for element in mollist:
		chargevec.extend(readcharges(element,top))
        term=str(str(gro[-3])+str(gro[-2])+str(gro[-1]))
        geo=[]
        if term=="g96":
            geo=readg96(gro)
        else:
	    geo=readgeo(gro)
	#print array(geo).reshape(-1,3)
	#print array(chargevec).reshape(-1,1)
#	testarray=[]
#	for i in range(1971,len(chargevec)):
#		testarray.append(chargevec[i])
#	print array(testarray).reshape(-1,1)
	if len(geo)!=3*len(chargevec):
		print "Not all atoms (" + str(len(geo)/3.) + ") were replaced by charges (" + str(len(chargevec)) + ") in full iteration step! Exiting."
		exit(1)
	makeout(geo,chargevec,out)
		
if __name__ == '__main__':
        import sys
	generate_pcf_from_top(sys.argv[1],sys.argv[2],sys.argv[3])
