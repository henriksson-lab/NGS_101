#!/usr/bin/env python3
"""Sample public, processed Spatial VDJ target-sequence TSV rows without a full ZIP.

Raw EGA FASTQs are controlled. This accesses only targetSequences fields in a
public Zenodo TSV export following the reported pRESTO/UMI-processing workflow. It cannot certify an
original indexing-primer order. Samples/manifests stay outside the checkout.
"""
from pathlib import Path
import argparse
import csv
import io
from collections import Counter
import hashlib
import json
import re
import struct
import sys
import urllib.request
import zlib
from datetime import datetime, timezone

HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE),str(HERE.parents[1]/'lib')]
from storage import data_root

URL='https://zenodo.org/records/7961605/files/SpatialVDJ_forZenodo_resubmission.zip?download=1'
PREFIX='SpatialVDJ_forZenodo/data/tonsil/2_SR-SpatialVDJ/'
MAX_DIRECTORY=1024*1024
MAX_OUTPUT=4*1024*1024


def part(start, end=None):
    if (end is None and start>=0) or (end is not None and (start<0 or end<start)):
        raise ValueError('Expected a negative suffix or an ordered non-negative byte interval')
    spec=str(start) if end is None else f'{start}-{end}'
    request=urllib.request.Request(URL, headers={'Range':'bytes='+spec,'User-Agent':'NGS101 bounded public alignment sample'})
    with urllib.request.urlopen(request, timeout=25) as response:
        match=re.fullmatch(r'bytes (\d+)-(\d+)/(\d+)',response.headers.get('Content-Range',''))
        if response.status!=206 or match is None:
            raise RuntimeError('Range not honored; refusing a full-archive download')
        actual_start,actual_end,total=map(int,match.groups())
        if total<=0 or not 0<=actual_start<=actual_end<total:
            raise RuntimeError('Invalid returned byte interval or total')
        expected=(max(0,total+start),total-1) if end is None else (start,end)
        if (actual_start,actual_end)!=expected:
            raise RuntimeError('Wrong byte interval returned')
        limit=actual_end-actual_start+1
        content=response.read(limit+1)
        if len(content)!=limit:raise RuntimeError('Wrong bounded response length')
        return content


def complete_lines(content, max_reads):
    # A range/deflate prefix may stop inside a UTF-8 character in the final row.
    # Drop that entire incomplete row as bytes before strict decoding; malformed
    # UTF-8 inside retained complete rows must still fail.
    end=content.rfind(b'\n')
    if end<0:raise RuntimeError('No complete lines in bounded prefix')
    lines=content[:end+1].decode('utf-8',errors='strict').splitlines()
    if len(lines)<2:raise RuntimeError('No complete alignment rows in bounded prefix')
    return lines[:max_reads+1]


def check_existing_sample(target, source, new_sample):
    manifest_path=target.with_suffix('.manifest.json')
    if not target.exists() and not manifest_path.exists():return
    if not target.is_file() or not manifest_path.is_file():
        raise RuntimeError('Existing sample/manifest pair is incomplete')
    prior=json.loads(manifest_path.read_text())
    if not isinstance(prior,dict):raise RuntimeError('Existing sample manifest is not an identity record')
    identity=('source_url','member','declared_crc32','compressed_bytes',
              'uncompressed_bytes','compression_method','local_header_offset')
    if any(prior.get(key)!=source.get(key) for key in identity):
        raise RuntimeError('Existing sample belongs to another source identity or member revision')
    with target.open('rb') as handle:existing=handle.read(MAX_OUTPUT+1)
    expected=prior.get('sample_sha256','')
    if (len(existing)>MAX_OUTPUT or not isinstance(expected,str)
            or re.fullmatch(r'[0-9a-f]{64}',expected)is None
            or hashlib.sha256(existing).hexdigest()!=expected):
        raise RuntimeError('Existing sample does not match its manifest hash')
    if hashlib.sha256(new_sample).hexdigest()!=expected:
        raise RuntimeError('New sample differs from the previously manifested sample')


def member_metadata(name):
    trailer=part(-131072)
    index=trailer.rfind(b'PK\x06\x06')
    if index<0:raise RuntimeError('Expected ZIP64 end record')
    values=struct.unpack_from('<4sQ2H2I4Q',trailer,index)
    size,offset=values[-2:]
    if size>MAX_DIRECTORY:raise RuntimeError('Directory byte limit exceeded')
    directory=part(offset,offset+size-1)
    index=0
    while index+46<=len(directory) and directory[index:index+4]==b'PK\x01\x02':
        n,e,c=struct.unpack_from('<3H',directory,index+28)
        filename=directory[index+46:index+46+n].decode('utf-8')
        if filename==name:
            compressed,uncompressed=struct.unpack_from('<2I',directory,index+20)
            local=struct.unpack_from('<I',directory,index+42)[0]
            method=struct.unpack_from('<H',directory,index+10)[0]
            crc=struct.unpack_from('<I',directory,index+16)[0]
            # Selected metadata has ordinary sizes/offset below ZIP64 sentinel;
            # refuse a changed member layout rather than guessing extra fields.
            if 0xffffffff in (compressed,uncompressed,local):raise RuntimeError('Selected member layout changed')
            return {'compressed_bytes':compressed,'uncompressed_bytes':uncompressed,
                    'local_header_offset':local,'compression_method':method,'declared_crc32':f'{crc:08x}'}
        index+=46+n+e+c
    raise RuntimeError('Exact source member not found')


def summarize(sample, sample_sha256):
    rows=list(csv.DictReader(io.StringIO(sample.decode()),delimiter='\t'))
    if 'targetSequences' not in rows[0]:raise RuntimeError('Expected source targetSequences field missing')
    sequences=[sequence.upper()for row in rows for sequence in row['targetSequences'].split(',')if sequence]
    from illumina import P5,P7,TRUSEQ_READ1,TRUSEQ_READ2
    from chemdraw import revcomp
    import oligos_spatial_vdj as oligos
    anchors={'P5':P5,'P7':P7,'TruSeq_Read1':TRUSEQ_READ1,'TruSeq_Read2':TRUSEQ_READ2,'source_PartRead1':oligos.PARTREAD1}
    summary={'sample_sha256':sample_sha256,'alignment_rows':len(rows),
             'target_sequences':len(sequences),'length_histogram':dict(Counter(map(len,sequences))),
             'exact_anchor_matches':{name:{'as_written':sum(anchor in seq for seq in sequences),
                 'reverse_complement':sum(revcomp(anchor) in seq for seq in sequences)}for name,anchor in anchors.items()},
             'partial_anchor_matches':{str(width):{name:{orientation:sum(any(candidate[k:k+width]in seq for k in range(len(candidate)-width+1))for seq in sequences)
                 for orientation,candidate in (('as_written',anchor),('reverse_complement',revcomp(anchor)))}for name,anchor in anchors.items()}for width in (12,16,20)},
             'interpretation':'First rows only, not a random or exhaustive sample. A partial technical-anchor match is only a candidate; matching final copied sequence cannot uniquely reconstruct original primer synthesis order or modifications.'}
    return summary


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--chain',choices=('TRA','TRB'),default='TRB')
    parser.add_argument('--max-reads',type=int,default=100,help='Maximum processed alignment rows; not original sequencer read count (1–100)')
    parser.add_argument('--max-compressed-bytes',type=int,default=131072,help='Strict compressed member-prefix cap, maximum131072')
    args=parser.parse_args()
    if not 1<=args.max_reads<=100:parser.error('--max-reads must be 1–100')
    if not 1024<=args.max_compressed_bytes<=131072:parser.error('compressed cap must be1024–131072')
    name=PREFIX+args.chain+'/Tonsil_merged.matrix.tsv'
    meta=member_metadata(name);offset=meta['local_header_offset']
    header=part(offset,offset+29)
    if header[:4]!=b'PK\x03\x04':raise RuntimeError('Wrong source local header')
    n,e=struct.unpack_from('<2H',header,26)
    if n+e>4096:raise RuntimeError('Oversized member header')
    filename=part(offset+30,offset+29+n+e)[:n].decode('utf-8')
    if filename!=name:raise RuntimeError('Source member identity mismatch')
    start=offset+30+n+e
    length=min(args.max_compressed_bytes,meta['compressed_bytes']-1)
    compressed=part(start,start+length-1)
    if meta['compression_method']!=8:raise RuntimeError('Unexpected source compression method')
    decoder=zlib.decompressobj(-15)
    content=decoder.decompress(compressed,MAX_OUTPUT)
    selected=complete_lines(content,args.max_reads)
    sample=('\n'.join(selected)+'\n').encode()
    folder=data_root()/'read_samples'/'spatial-vdj__10.1126+science.adf8486'/'index-tail-probe-2026-11-16'
    if folder.resolve().is_relative_to(HERE.parents[1]):
        raise RuntimeError('Read samples must stay outside the repository checkout')
    folder.mkdir(parents=True,exist_ok=True)
    target=folder/(args.chain+'_first_'+str(args.max_reads)+'_processed_alignments.tsv')
    manifest={'source_url':URL,'record':'10.5281/zenodo.7961605','member':name,**meta,
        'scope':'Public TSV targetSequences export rows after reported pRESTO/UMI processing; not inspected VDJCA saved-input contents, raw EGA FASTQs or original-primer evidence.',
        'requested_max_reads':args.max_reads,'sampled_alignment_rows':len(selected)-1,
        'member_prefix_start':start,'member_prefix_bytes':length,'max_decompressed_bytes':MAX_OUTPUT,
        'full_member_crc_verified':False,'partial_crc_limit':'Declared whole-member CRC retained from exact central-directory entry; only a bounded prefix was read.',
        'sample_sha256':hashlib.sha256(sample).hexdigest(),'sample_file':target.name,
        'retrieved_at':datetime.now(timezone.utc).isoformat()}
    summary=summarize(sample,manifest['sample_sha256'])
    check_existing_sample(target,manifest,sample)
    target.write_bytes(sample)
    target.with_suffix('.manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    target.with_suffix('.summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps({'file':str(target),'rows':len(selected)-1,'columns':selected[0].split('\t'),'prefix_bytes':length}))
    print(json.dumps(summary))


if __name__=='__main__':main()
