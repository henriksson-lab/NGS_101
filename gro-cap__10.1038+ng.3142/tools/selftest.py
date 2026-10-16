#!/usr/bin/env python3
"""Compare local declarations with original source text when available."""
from pathlib import Path
import os
import re
import sys
import xml.etree.ElementTree as ET

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / 'lib')]
import chem
from checks import Check, Source, have

c = Check()
c.section('GRO-cap source transcription')
bundle = Path(os.environ.get('CHEM_DATA', HERE.parents[1] / '_data')) / 'sources' / HERE.parent.name

vendor = Source(bundle / 'illumina_TruSeq-SmallRNA.htm.txt',
                'https://support-docs.illumina.com/SHARE/AdapterSequences/Content/TruSeq-SmallRNA.htm')
if have(c, vendor, label='current vendor adapter references'):
    text = vendor.path.read_text()
    for name, declared in (('RA3', chem.RA3), ('RA5', chem.RA5)):
        match = re.search(r'Adapter \(' + name + r'\)\s+5′\s+([ACGTU]+)', text)
        c(name + ' vendor entry located', match is not None)
        if match:
            c(name + ' current reference bases', declared, match[1].replace('U', 'T'))

kruesi = Source(bundle / 'elife00808.pdf.txt', 'https://doi.org/10.7554/eLife.00808')
if have(c, kruesi, label='original predecessor mapping method'):
    text = kruesi.path.read_text().split('Mapping of GRO-seq and GRO-cap reads', 1)[1]
    target = re.search(r'cutadapt -a ([ACGT]+)', text)
    c('original joint mapping adapter target located', target is not None)
    if target:
        c('RA3 reference agrees with original read-through trimming target', chem.RA3, target[1])

core = Source(bundle / '594009_9144e1ea58ad4573871ec024827290f3.pdf.txt',
              'https://www.cecerelab.com/_files/ugd/594009_9144e1ea58ad4573871ec024827290f3.pdf')
if have(c, core, label='Core 2008 original inherited completion method'):
    text = core.path.read_text().split('Reverse transcription and amplification and PAGE purification of NRO-RNA', 1)[1]
    rt = re.search(r'at (\d+)°C for (\d+)min, followed\s+by (\d+)°C for (\d+) min', text)
    c('RT programme source entry located', rt is not None)
    if rt:
        a, b, d, e = map(int, rt.groups())
        c('inherited RT programme', chem.RT_PROGRAM, ((a, b), (d, e)))
    cycles = re.search(r'amplified (\d+) cycles', text)
    c('PCR cycle source entry located', cycles is not None)
    if cycles:
        c('inherited PCR cycles', chem.PCR_CYCLES, int(cycles[1]))
    threshold = re.search(r'cDNAs greater than (\d+)\s+nucleotides', text)
    c('PAGE selection source entry located', threshold is not None)
    if threshold:
        c('inherited PAGE lower threshold', chem.PAGE_MIN_NT, int(threshold[1]))

experiments=[Source(bundle / f'SRX{number}.xml',
                    f'https://www.ebi.ac.uk/ena/browser/api/xml/SRX{number}?format=xml')
             for number in range(682016,682020)]
if have(c,*experiments,label='deposited Core 2014 acquisition metadata'):
    originals=[ET.parse(source.path).getroot().find('EXPERIMENT') for source in experiments]
    c('experiment originals retain their requested accession identities',
      [x.get('accession') for x in originals], [s.path.stem for s in experiments])
    layouts={child.tag for x in originals for child in x.find('.//LIBRARY_LAYOUT')}
    c('single-end declaration follows every deposited Core GRO-cap experiment',
      (layouts,chem.ACQUIRED_ROLES), ({'SINGLE'},('Read 1',)))

c.report()
