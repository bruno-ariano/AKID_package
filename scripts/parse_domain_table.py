import argparse
from itertools import groupby
import re

usage = """%(prog)s reads .domtblout file and returns non-overlapped (or
specified percent of overlapping) domains below given e-value.
"""
p = argparse.ArgumentParser(description=usage)
p.add_argument("-f", "--f", "-file", "--file", dest="file",
                  help=".domtblout file")
p.add_argument("-e", "--e", "--evalue", "-evalue", type=float, dest="evalue",
                   help="E-value cut-off", default=1e-05)

args = p.parse_args()

fh = open(args.file)
oh = open(args.file + '.out', 'w')

grouper = lambda x: re.split('\s+', x)[3] if not x.startswith('#') else x[0]
# Iterate through proteins.
for k, g in groupby(fh, grouper):
    if k.startswith('#'): continue  # Skip header and footer of a domtblout file.
    l = []
    for line in g:
        sl = re.split('\s+', line)
        #print sl
        pid = sl[3]
        dname = sl[0]
        did = sl[1]
        hmmstart=int(sl[15])
        hmmend=int(sl[16])
        dstart = int(sl[17])
        dend = int(sl[18])
        evalue = float(sl[12])
        if evalue <= args.evalue:
             l.append([dname,hmmstart,hmmend, did, dstart, dend, evalue])
    
    # Filter domains by a given E-value and overlapping cut-off.
    filtered=[]
    # Save domains to a file.
    #filtered.sort(key=lambda x: x[2])
    for i in range(0, len(l)):
        if i:
            if l[i][1]>=l[i-1][2]-200 and l[i][4]-l[i-1][5]<=200:   
                #val1=max(l[i-1][5],l[i][5])
                l[i-1][5]=l[i][5]
                continue
        filtered.append(l[i])

    # Save domains to a file.
    filtered.sort(key=lambda x: x[2])
    for d in filtered:
        oh.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n'.format(pid, *d))

fh.close()
oh.close()
