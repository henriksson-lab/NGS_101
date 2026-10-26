#!/usr/bin/env python3
"""Build the unlisted source-request page from the outstanding-document ledger."""
from collections import OrderedDict
from dataclasses import dataclass
import html
from pathlib import Path
import re
import sys
from urllib.parse import urlsplit

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / 'lib'))
from page import head

OUT = ROOT / 'source-downloads.html'
URL = re.compile(r'https?://[^\s<>]+')


@dataclass(frozen=True)
class Request:
    number: int
    title: str
    directory: str
    links: tuple[tuple[str, str], ...]
    details: str


def parse_requests(text):
    starts = list(re.finditer(r'(?m)^(\d+)\.\s+(.+)$', text))
    requests = []
    for i, start in enumerate(starts):
        block = text[start.end():starts[i+1].start() if i+1 < len(starts) else len(text)]
        destination = re.search(r'(?m)^\s+Directory:\s*(\S+)\s*$', block)
        if destination is None:
            raise ValueError(f'request {start[1]} needs a Directory field')
        directory = destination[1]
        if not re.fullmatch(r'[A-Za-z0-9_+().-]+', directory) or directory in {'.', '..'}:
            raise ValueError(f'unsafe source directory: {directory}')
        lines = [line.strip() for line in block[:destination.start()].strip().splitlines()]
        title = start[2].rstrip(':')
        # Preserve a wrapped title/revision before the first link label or URL.
        while lines and lines[0] and not URL.search(lines[0]) and not lines[0].endswith(':'):
            if re.match(r'(Status|Needed|DOI|PDF|Published original):', lines[0]):
                break
            title += ' ' + lines.pop(0)
        links = OrderedDict()
        previous = ''
        for line in lines:
            for match in URL.finditer(line):
                url = match[0]
                prefix = line[:match.start()].strip().rstrip(':')
                label = prefix or (previous.rstrip(':') if previous.endswith(':') else '')
                if not label:
                    parsed = urlsplit(url)
                    label = 'DOI' if parsed.netloc.lower() == 'doi.org' else 'PDF' if parsed.path.lower().endswith('.pdf') else parsed.netloc
                links.setdefault(url, label)
            previous = line
        if not links:
            raise ValueError(f'request {start[1]} needs at least one source URL')
        requests.append(Request(int(start[1]), title, directory,
                                tuple(links.items()), '\n'.join(lines)))
    if len({r.number for r in requests}) != len(requests):
        raise ValueError('request numbers must be unique')
    return tuple(requests)


def render_page(text):
    requests = parse_requests(text)
    groups = OrderedDict()
    for request in requests:
        groups.setdefault(request.directory, []).append(request)
    parts = ['<!doctype html><html lang="en"><head>',
             '<meta name="robots" content="noindex,nofollow">',
             head('Source downloads · NGS 101', home='index.html', include_header=False),
             '<style>.source-links{display:flex;flex-wrap:wrap;gap:8px 16px;list-style:none;padding:0}'
             '.source-details{white-space:pre-line;overflow-wrap:anywhere}'
             '.source-request{border-top:1px solid var(--rule);padding:16px 0}'
             '.source-directory{color:var(--ink-muted);overflow-wrap:anywhere}</style>',
             '</head><body><main class="wrap"><h1>Source downloads</h1>',
             '<p>Outstanding originals and supplements. Use the specified published version or manual revision.</p>',
             f'<p>{len(requests)} document requests · {len(groups)} protocols</p>']
    for directory, group in groups.items():
        for request in group:
            links = ''.join(f'<li><a href="{html.escape(url, quote=True)}">{html.escape(label)}</a></li>'
                            for url, label in request.links)
            parts.append(f'<section class="source-request" id="request-{request.number}">'
                         f'<h2>{html.escape(request.title)}</h2>'
                         f'<p class="source-directory">Source folder: <code>{html.escape(directory)}</code></p>'
                         f'<ul class="source-links">{links}</ul>'
                         '<details><summary>Request details</summary>'
                         f'<p class="source-details">{html.escape(request.details)}</p></details></section>')
    if not requests:
        parts.append('<p>No outstanding document requests.</p>')
    return '\n'.join([*parts, '</main></body></html>']) + '\n'


def main():
    OUT.write_text(render_page((ROOT / 'geturl.txt').read_text(encoding='utf-8')), encoding='utf-8')
    print(OUT.name)


if __name__ == '__main__':
    main()
