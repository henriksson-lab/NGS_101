#!/usr/bin/env python3
"""Offline failure-mode checks for the bounded public TSV sampler; no network."""
from pathlib import Path
import hashlib
import io
import json
import sys
from tempfile import TemporaryDirectory
from unittest.mock import patch

HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE),str(HERE.parents[1]/'lib')]
from checks import Check
import sample_public_alignments as sampler


class Response:
    def __init__(self, content_range, body, status=206):
        self.status=status;self.headers={'Content-Range':content_range}
        self.body=io.BytesIO(body);self.read_limits=[]
    def __enter__(self):return self
    def __exit__(self,*args):self.body.close()
    def read(self,limit):self.read_limits.append(limit);return self.body.read(limit)


def rejected(fn,error=RuntimeError):
    try:fn()
    except error:return True
    return False


def run(c):
    def request(response,start,end=None):
        with patch.object(sampler.urllib.request,'urlopen',return_value=response):
            return sampler.part(start,end)
    body=b'x'*100
    response=Response('bytes 900-999/1000',body)
    c('Correct archive suffix is accepted',request(response,-100),body)
    c('Suffix body read is bounded',response.read_limits,[101])
    c('Equal-length wrong suffix is rejected',rejected(lambda:request(Response('bytes 800-899/1000',body),-100)),True)
    c('Suffix inconsistent with reported total is rejected',rejected(lambda:request(Response('bytes 900-999/1100',body),-100)),True)
    c('Unknown archive total cannot validate a suffix',rejected(lambda:request(Response('bytes 900-999/*',body),-100)),True)
    c('Requested suffix larger than archive uses complete reported interval',request(Response('bytes 0-9/10',b'abcdefghij'),-100),b'abcdefghij')
    c('Forward range checks exact offsets too',rejected(lambda:request(Response('bytes 20-29/100',b'0123456789'),10,19)),True)
    response=Response('bytes 0-9/10',b'abcdefghij',status=200)
    c('Non-partial response is refused before body reading',rejected(lambda:request(response,-10)),True)
    c('Refused full response body is untouched',response.read_limits,[])
    c('Truncated bounded body is rejected',rejected(lambda:request(Response('bytes 900-999/1000',b'x'*99),-100)),True)

    header=b'targetSequences\tdescription\n'
    complete='ACGT\tvalid α\n'.encode()
    prefix=header+complete+b'ACGT\tpartial \xce'
    c('Split trailing UTF8 character does not discard valid rows',sampler.complete_lines(prefix,100),['targetSequences\tdescription','ACGT\tvalid α'])
    c('Invalid UTF8 in a complete row is still rejected',rejected(lambda:sampler.complete_lines(header+b'ACGT\t\xce\n',100),UnicodeDecodeError),True)
    c('Header plus incomplete row is not a read sample',rejected(lambda:sampler.complete_lines(header+b'ACGT',100)),True)
    c('Complete prefix is capped at requested rows',len(sampler.complete_lines(header+b'ACGT\tone\nACGT\ttwo\n',1)),2)

    with TemporaryDirectory() as tmp:
        target=Path(tmp)/'sample.tsv';manifest=target.with_suffix('.manifest.json')
        sample=header+b'ACGT\tfixture\n'
        source={'source_url':'https://example.org/immutable.zip','member':'processed.tsv',
                'declared_crc32':'0123abcd','compressed_bytes':50,'uncompressed_bytes':100,
                'compression_method':8,'local_header_offset':500}
        prior={**source,'sample_sha256':hashlib.sha256(sample).hexdigest()}
        def reset():target.write_bytes(sample);manifest.write_text(json.dumps(prior))
        reset()
        c('Unchanged source and verified sample may be reused',sampler.check_existing_sample(target,source,sample),None)
        for field,new in [('declared_crc32','deadbeef'),('compressed_bytes',51),('uncompressed_bytes',101)]:
            changed={**source,field:new}
            c('Changed member '+field+' is refused at the same URL/name',rejected(lambda:sampler.check_existing_sample(target,changed,sample)),True)
        c('A changed source URL is still rejected',rejected(lambda:sampler.check_existing_sample(target,{**source,'source_url':'https://other.org/data.zip'},sample)),True)
        target.write_bytes(sample+b'corrupted')
        c('Corrupted prior sample is detected against its manifest hash',rejected(lambda:sampler.check_existing_sample(target,source,sample)),True)
        c('Rejected corrupted sample is not overwritten',target.read_bytes(),sample+b'corrupted')
        reset()
        c('A changed incoming sample cannot silently replace manifested bytes',rejected(lambda:sampler.check_existing_sample(target,source,sample+b'changed')),True)
        c('Prior sample is preserved after incoming hash mismatch',target.read_bytes(),sample)
        manifest.unlink()
        c('Existing sample without identity manifest is refused',rejected(lambda:sampler.check_existing_sample(target,source,sample)),True)


if __name__=='__main__':
    checks=Check();run(checks);checks.report()
