#!/usr/bin/env python

#e.g. python unique_fn_parts.py -s $INPUT_FQ_SUFFIX -d input_files -o input_filenames_auto.txt

import os, glob, getopt, sys

options, remainder = getopt.getopt(sys.argv[1:], 's:d:o:', ['suffix=', 
                                                         'dir=',
                                                         'outfile=',
                                                         ])
for opt, arg in options:
    if opt in ('-s', '--suffix'):
        search_suffix = arg
    elif opt in ('-d', '--dir'):
        indir = arg
    elif opt in ('-o','--outfile'):
        outfn = arg

TCPdir=os.getcwd()
os.chdir(indir)
allfiles=glob.glob('*'+search_suffix)
print str(len(allfiles))+' files found: in \''+indir+ '\' with suffix \'' + search_suffix + '\':'

com_pref = len(os.path.commonprefix(allfiles))
com_suff = len(os.path.commonprefix([x[::-1] for x in allfiles]))

os.chdir(TCPdir)
with open(outfn, 'w') as fout:
	print >> fout, '\t'.join(['file_name','file_uniq','fraction', 'sample'])
	for f in allfiles:
		uniqpart = f[com_pref:len(f)-com_suff]
		print >> fout, '\t'.join([f, uniqpart])
