"""Source-linked index oligos, not a whole Ståhl/Visium library recipe.

Engblom final SR indexing Methods explicitly delegates to Ståhl 2016 ref.6.
Ståhl published supplement PDF pp8–9 gives three PCR oligos. The six Xs
are an unassigned selected index, not a random or degenerate reagent pool.
Concentrations/enzymes from 2016 are not silently assigned to the 2023 reaction.
"""
import illumina
from chemdraw import Scene, Segment, feature

SOURCE_URL = 'https://doi.org/10.1126/science.aaf2403'
SOURCE = 'Ståhl et al.2016 published supplement, PDF pp8–9, library preparation and PCR oligo list; Engblom2023 final SR indexing Methods reference6'
INPE1 = illumina.TRUSEQ_P5_FULL
INPE2 = illumina.TRUSEQ_READ2
INDEX_FOOT = illumina.TRUSEQ_READ2[:13]
INDEX_NT = 6
SAMPLE = feature('spatialvdj.sample','sample_index','unknown',
    note='Six-base selected index in the cited Ståhl PCR Index order; actual sample assignment unassigned; no i5 index added')


def index_segments():
    return [Segment('P7',illumina.P7,'p7'),
            Segment('sample index in primer order','X'*INDEX_NT,'idx',placeholder=True,
                    length_bp=INDEX_NT,feature=SAMPLE),
            Segment('3′ Read2 binding foot',INDEX_FOOT,'r2')]


def primer_scene():
    sc=Scene()
    sc.strand('InPE1.0',[Segment('InPE1.0',INPE1,'r1')])
    sc.strand('InPE2.0',[Segment('InPE2.0',INPE2,'r2')])
    sc.strand('PCR_Index',index_segments(),label='PCR Index; XXXXXX = chosen six-base index, not a mixed oligo')
    return sc
