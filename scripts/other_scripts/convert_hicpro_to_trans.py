import argparse, re, sys
from argparse import RawDescriptionHelpFormatter

#def getOptions():
description="""This script converts hic-pro valid pairs to trans reads. Input: SAM and a min percent"""

parser = argparse.ArgumentParser(description=description, formatter_class=RawDescriptionHelpFormatter)
parser.add_argument('-i','--input', dest="input", action='store', required=True,  help='Input file [REQUIRED]', metavar="VALIDPAIRS")
parser.add_argument('-o','--output', dest="output", action='store', required=True, help='Output file [REQUIRED]', metavar="OUTTRANS")
args = parser.parse_args()

inhandle = args.input
outhandle = args.output
with open(inhandle) as infile, open(outhandle,'w') as outfile:
	for line in infile:
		fields = line.strip().split('\t')
		if fields[1] != fields[4]:
			outfile.write('{0}\tNA\t{1}\t{2}\tNA\tNA\t{3}\tNA\t{4}\t{5}\tNA\tNA\n'.format(fields[0],fields[1],fields[2],fields[0],fields[4],fields[5]))
