#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pysam
import sys
from math import log

minlen = 300

samfile = pysam.AlignmentFile(sys.stdin)
reftotlen = sum(samfile.lengths)
rdstotlen = 0

reads = []
fragments = []
for read in samfile.fetch(until_eof=True):
    if read.is_unmapped:
        rdstotlen += read.query_length
    else:
        if not read.is_secondary:
            rdstotlen += read.infer_read_length()
        if read.reference_length>=minlen:
             reads.append(read)
             fragments.append((read.reference_start, read.reference_end))
samfile.close()

contig_overlaps = []
for i, frag1 in enumerate(fragments):
    for frag2 in fragments[i+1:]:
        s,e = max(frag1[0], frag2[0]),min(frag1[1], frag2[1])
        if s<e:
            contig_overlaps.append((s,e))
contig_overlaps.sort()

redundant = []
ls,le = -1,-1
for cs,ce in contig_overlaps:
    if cs<=le:
        le = max(le, ce)
    else:
        if le>=0:
            redundant.append((ls,le))
        ls,le = cs,ce
if le>=0:
    redundant.append((ls,le))

almtotlen = 0
almmmcount = 0
almcount = 0
for read in reads:
    s,e = read.reference_start, read.reference_end
    alms = s
    alpairs = read.get_aligned_pairs(True, True)
    for rs, re in redundant + [(e,e)]:
        if alms<rs:
            alme = min(e, rs)
            if alme-alms >= minlen:
                almtotlen += alme-alms
                almcount+=1
                i = alms-s
                while alpairs[i][1] is None or alpairs[i][1]<alms:
                    i += 1
                while i<len(alpairs) and (alpairs[i][1] is None or alpairs[i][1]<alme):
                    if alpairs[i][1] is None or alpairs[i][2] is None or alpairs[i][2].islower():
                        almmmcount += 1
                    i += 1
        if re<e:
            alms = max(alms, re)
        else:
            break

almtotlen = float(almtotlen)
refcoverage = almtotlen/reftotlen if almtotlen else 0
rdscoverage = almtotlen/rdstotlen if almtotlen else 0
ident_score = max(0.5, 1-10*almmmcount/almtotlen) if almtotlen else 0
count_score = 1/log(4+almcount, 5) if almtotlen else 0

score = refcoverage*rdscoverage*ident_score*count_score

print("Pokrycie referencji:", refcoverage)
print("Pokrycie odczytów:"  , rdscoverage)
#print("Błędy uliniowień:"   , almmmcount/almtotlen if almtotlen else 0)
print("Ocena identyczności:", ident_score)
#print("Liczba uliniowień:"  , almcount)
print("Ocena rozdrobnienia:", count_score)
print("Łączna ocena:"       , score)
