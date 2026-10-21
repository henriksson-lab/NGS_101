"""Synthetic molecular examples for paper figures, not protocol/source facts."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'lib'))

from chemdraw import Construct, Scene, Segment, complement_segments, feature, revcomp
from topology import HairpinTopology, PadlockCaptureTopology


def annotation_example():
    parts = [
        Segment('handle', 'ACGTCGATGCAT', 'p5'),
        Segment('UMI', 'NNNN', 'umi', placeholder=True,
                feature=feature('paper_umi', 'umi', 'random')),
        Segment('insert', 'GATCAGTACG'),
        Segment('tail', 'TTA', 'r2', inferred=True),
    ]
    construct = Construct(parts, name='synthetic annotation example')
    scene = Scene()
    scene.strand('top', parts, mod5='P')
    scene.anneal('bottom', complement_segments(parts), to='top', pair=("handle'", 'handle'))
    scene.mark('top', 'UMI', 'UMI span')
    scene.mark('bottom', "tail'", 'inferred tail')
    return construct, scene


def padlock_example():
    ext_site = 'ACGTCAGT'
    lig_site = 'TGCATCGA'
    return PadlockCaptureTopology(
        'synthetic padlock', 'TT' + ext_site + 'GGAT' + lig_site + 'AA',
        revcomp(lig_site), 'GACTAACGTAGCTAGTCA', revcomp(ext_site))


def hairpin_example():
    stem = 'ACGTCAGT'
    loop = 'GATTACAG'
    prefix, suffix = 'TTA', 'AAC'
    return HairpinTopology('synthetic hairpin',
                           prefix + stem + loop + revcomp(stem) + suffix,
                           len(prefix), len(prefix) + len(stem) + len(loop), len(stem))
