#! /usr/bin/env python
"""This script was developed to get the average profile 
of specific genome regions, which could do by sitepro or 
siteproBW in CEAS suite.
Xiaoqin Yang @ Tongji University Feb-7-2013"""


import re

def linenumber(infile):
	count = 0
	for line in open(infile).xreadlines(): 
		count += 1
	return count


def chrprof(wigfile):
	values = {}
	for line in open(wigfile):	
		if re.search("track", line):
			continue
		elif re.search("chrom=chr", line):
			parts = line.strip().split("\t")
			pieces = parts[1].split("=")
			chrname = pieces[1]
			values[chrname] = []
		elif re.search('^\d+', line):
			 values[chrname].append(int(line.strip()))
		else:
			continue
	return values	 


def average(bedfile, wigfile, span):

	pileup = chrprof(wigfile)	
	
	with open(bedfile, 'r') as f:
		first_line = f.readline()

		parts = first_line.strip().split("\t")
		chrname = parts[0]
		startp = int(parts[1])
		endp = int(parts[2])
		midp = float(startp + endp)/2
		if float(startp + endp)%2 == 0:
			profile = [0] * (int(span) * 2 )
		else:
			profile = [0] * (int(span) * 2 + 1)

	
	for line in open(bedfile):
		parts = line.strip().split("\t")
		chrname = parts[0]
		startp = int(parts[1])
		endp = int(parts[2])
		midp = float(startp + endp)/2
		if float(startp + endp)%2 == 0:
			start = midp - (span - 1)
			end = midp + (span - 1)				
		else:
			start = midp - 0.5 - (span -1)
			end = midp + 0.5 + (span -1) 
		
		


	linenum = linenumber(bedfile)
	i = 0




	while i < linenum:
		
		for location in profile:
			profile[i] += location[i]
		i += 1			

	return profile


if __name__ == '__main__':
	import sys
	from optparse import OptionParser

	usage = "usage: python average_profile.py -c chr22 -i input.bed -o output.R"
	parser = OptionParser(usage)
	parser.add_option("-b",dest="bedfile",help="input Bed file for peak region")
	parser.add_option("-w",dest="wigfile",help="input wig file for reads")
	parser.add_option("-s",dest="span",help="Span from the center of each BED region in both directions")
	parser.add_option("-o",dest="output", help="output R script to plot the average profile")

	(options, args) = parser.parse_args()    
	if not options.spename or not options.inbed or not options.outwig:
        	parser.print_help()
        	sys.exit()      
	
	o = open(options.output, 'w')
	
	averpro(options.bedfile, options.wigfile, options.span)
	
	