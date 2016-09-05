#!/usr/bin/env python
import sys
fasta=sys.stdin.read()
fasta_list = fasta.split('>')
fasta_uniq = dict()
for entry in fasta_list:
    if entry.rstrip('\n') == '':
        continue
    output = '>' + entry.rstrip('\n')
    if output[-3:] == "CCA":
        output = output[0:-1]
    else:
        output += 'CC'
    splut = output.split("\n",1)
    splut[1] = splut[1].upper()
    fasta_uniq[splut[1]]=splut[0]  # key=seq; value=name
for seq in fasta_uniq:
    print fasta_uniq[seq]+"\n"+seq
    