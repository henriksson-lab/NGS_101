#!/usr/bin/env python3
"""Opt-in Chrome restriction-tab checks; uses standard-library CDP and a local server.

Run after building: python3 tools/benchmark_restriction_browser.py --site _site
No packages are installed. A temporary profile, deliberate fetch failures and a
GitHub-Pages-style URL prefix exercise workers, filters, exports and cancellation.
"""
import argparse
import base64
import functools
import http.server
import json
import os
from pathlib import Path
import socket
import struct
import subprocess
import tempfile
import threading
import time
import urllib.request

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--site',type=Path,default=Path(__file__).resolve().parents[1]/'_site')
parser.add_argument('--browser',default='google-chrome')
parser.add_argument('--screenshot',type=Path,help='Save the initial duplex view as a PNG for visual inspection')
parser.add_argument('--report',type=Path,default=Path('/tmp/restriction-browser-report.json'))
args=parser.parse_args()
if not (args.site/'restriction.html').is_file():parser.error('Build the restriction page/site first.')

class WS:
 def __init__(self,url):
  from urllib.parse import urlparse
  u=urlparse(url);self.s=socket.create_connection((u.hostname,u.port),timeout=5);self.buf=b''
  key=base64.b64encode(os.urandom(16)).decode()
  self.s.sendall(('GET '+u.path+' HTTP/1.1\r\nHost: '+u.netloc+'\r\nUpgrade: websocket\r\nConnection: Upgrade\r\nSec-WebSocket-Key: '+key+'\r\nSec-WebSocket-Version: 13\r\n\r\n').encode())
  data=b''
  while b'\r\n\r\n' not in data:data+=self.s.recv(4096)
  header,self.buf=data.split(b'\r\n\r\n',1)
  assert b' 101 ' in header,header
 def read(self,n):
  while len(self.buf)<n:
   chunk=self.s.recv(max(n-len(self.buf),4096))
   if not chunk:raise ConnectionError('Browser debugger disconnected')
   self.buf+=chunk
  result,self.buf=self.buf[:n],self.buf[n:];return result
 def recv(self):
  a,b=self.read(2);length=b&127
  if length==126:length=struct.unpack('!H',self.read(2))[0]
  elif length==127:length=struct.unpack('!Q',self.read(8))[0]
  mask=self.read(4) if b&128 else None
  data=self.read(length)
  if mask:data=bytes(v^mask[i%4] for i,v in enumerate(data))
  return json.loads(data)
 def send(self,value):
  data=json.dumps(value).encode();length=len(data);header=bytes([129,128|length]) if length<126 else bytes([129,254])+struct.pack('!H',length)
  mask=os.urandom(4);self.s.sendall(header+mask+bytes(v^mask[i%4] for i,v in enumerate(data)))
 def call(self,method,params=None):
  self.ident=getattr(self,'ident',0)+1;ident=self.ident;self.send({'id':ident,'method':method,'params':params or {}})
  while True:
   reply=self.recv()
   if reply.get('id')==ident:return reply
class Handler(http.server.SimpleHTTPRequestHandler):
 def log_message(self,*args):pass
 def do_GET(self):
  if not self.path.startswith('/NGS_101/'):
   self.send_error(404);return
  self.path=self.path[len('/NGS_101'):]
  if self.path=='/assets/restriction_worker.js?forced-import-error=1':
   self.send_response(200);self.send_header('Content-Type','text/javascript');self.end_headers();self.wfile.write(b'importScripts("forced-missing-helper.js");');return
  if self.path=='/assets/sequence_annotation.js' and getattr(self.server,'fail_worker_once',False):
   self.server.fail_worker_once=False;self.send_error(503);return
  if self.path=='/assets/restriction_catalogue.json' and self.server.fail_once:
   self.server.fail_once=False;self.send_error(503);return
  super().do_GET()
server=http.server.ThreadingHTTPServer(('127.0.0.1',0),functools.partial(Handler,directory=str(args.site.resolve())));server.fail_once=True
threading.Thread(target=server.serve_forever,daemon=True).start()
checks=[]
with tempfile.TemporaryDirectory(prefix='chem-restriction-browser-') as profile:
 browser=subprocess.Popen([args.browser,'--headless','--no-sandbox','--disable-gpu','--disable-dev-shm-usage','--user-data-dir='+profile,'--remote-debugging-port=0','--window-size=1440,1000','about:blank'],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
 try:
  portfile=Path(profile)/'DevToolsActivePort';deadline=time.monotonic()+30
  while not portfile.exists():
   if browser.poll() is not None:raise RuntimeError('Chrome exited')
   if time.monotonic()>deadline:raise TimeoutError('No debugger port')
   time.sleep(.05)
  base='http://127.0.0.1:'+portfile.read_text().splitlines()[0]
  tab=json.load(urllib.request.urlopen(urllib.request.Request(base+'/json/new?about:blank',method='PUT')));ws=WS(tab['webSocketDebuggerUrl']);ws.call('Runtime.enable');ws.call('Page.enable');ws.call('Network.enable');ws.call('Network.setCacheDisabled',{'cacheDisabled':True})
  ws.call('Page.addScriptToEvaluateOnNewDocument',{'source':'''window.__errors=[];window.addEventListener('error',e=>window.__errors.push(e.message));window.addEventListener('unhandledrejection',e=>window.__errors.push(String(e.reason)));window.__blobs=[];const create=URL.createObjectURL.bind(URL);URL.createObjectURL=b=>{window.__blobs.push(b);return create(b)};const click=HTMLAnchorElement.prototype.click;HTMLAnchorElement.prototype.click=function(){if(this.download){window.__downloadName=this.download;return}return click.call(this)};window.__workerCount=0;const W=window.Worker;window.Worker=class extends W{constructor(...args){if(window.__failNextWorkerStart){window.__failNextWorkerStart=false;args[0]+='?forced-import-error=1'}super(...args);window.__workerCount++}};window.__beats=0;setInterval(()=>window.__beats++,10);'''})
  def evaluate(js,await_promise=False):
   r=ws.call('Runtime.evaluate',{'expression':js,'returnByValue':True,'awaitPromise':await_promise})
   if 'exceptionDetails' in r.get('result',{}):raise RuntimeError(r)
   return r.get('result',{}).get('result',{}).get('value')
  def wait(js):
   deadline=time.monotonic()+40
   while not evaluate(js):
    if time.monotonic()>deadline:raise TimeoutError(js+'; '+str(evaluate("document.querySelector('#restriction-status')?.textContent")))
    time.sleep(.04)
  def check(label,js):
   assert evaluate(js),label+' '+str(evaluate('window.__errors'));checks.append(label)
  def setv(name,value,event='change'):
   evaluate(f"(()=>{{const e=document.getElementById('restriction-{name}');e.value={json.dumps(value)};e.dispatchEvent(new Event('{event}',{{bubbles:true}}));return true}})()")
  def click(name):evaluate("document.getElementById('restriction-"+name+"').click();true")
  def submit():evaluate("document.getElementById('restriction-form').requestSubmit();true")
  def ready():wait("!document.getElementById('restriction-results').hidden && document.getElementById('restriction-cancel').disabled")
  ws.call('Page.navigate',{'url':'http://127.0.0.1:'+str(server.server_port)+'/NGS_101/restriction.html'});wait("document.readyState==='complete' && !!window.RestrictionAnalysis")
  check('restriction tab and all navigation work under GitHub Pages prefix',"[...document.querySelectorAll('nav a')].every(a=>new URL(a.href).pathname.startsWith('/NGS_101/')) && !!document.querySelector('nav a[href=\"restriction.html\"]')")
  click('example');submit();wait("document.getElementById('restriction-status').textContent.includes('Could not load') && !document.getElementById('restriction-analyze').disabled")
  checks.append('catalogue loading failure is recoverable');submit();ready()
  check('example scans in worker and displays certain unique cutters',"window.__workerCount>0 && document.querySelectorAll('#restriction-rows > tr').length>0 && [...document.querySelectorAll('#restriction-rows > tr')].every(r=>r.children[2].textContent==='1')")
  check('table is collapsed and interactive duplex is the main view',"!document.getElementById('restriction-table').open && !!document.querySelector('#restriction-duplex svg') && document.querySelectorAll('#restriction-duplex [data-strand=top]').length===80")
  if args.screenshot:
   evaluate("document.querySelector('.restriction-viewer').scrollIntoView();true")
   shot=ws.call('Page.captureScreenshot',{'format':'png'})['result']['data'];args.screenshot.write_bytes(base64.b64decode(shot))
  check('licensed snapshot attribution and complete license are visible on demand',"document.getElementById('restriction-data-info').textContent.includes('1088') || document.getElementById('restriction-data-info').textContent.includes('1,088')")
  check('copyright permission notice is bundled',"document.getElementById('restriction-license').textContent.includes('Permission to use, copy, modify, and distribute')")
  evaluate('window.__failNextWorkerStart=true');submit();wait("document.getElementById('restriction-status').textContent.includes('worker could not start') && !document.getElementById('restriction-analyze').disabled");check('worker import failure is handled without an uncaught error',"window.__errors.length===0");submit();ready();checks.append('worker startup recovers after retry')
  click('all');setv('search','HpaII,MspI','input');check('name search combines with all counts',"[...document.querySelectorAll('#restriction-rows .restriction-name')].map(e=>e.textContent).join(',')==='HpaII,MspI'")
  setv('cpg','unaffected');check('CpG-insensitive filter selects MspI independently of shared motif',"document.querySelectorAll('#restriction-rows > tr').length===1 && document.querySelector('#restriction-rows .restriction-name').textContent==='MspI'")
  setv('cpg','sensitive');check('CpG-sensitive filter selects HpaII',"document.querySelectorAll('#restriction-rows > tr').length===1 && document.querySelector('#restriction-rows .restriction-name').textContent==='HpaII'")
  setv('cpg','any');setv('context','overlap_cpg');check('actual motif context includes both isoschizomers',"document.querySelectorAll('#restriction-rows > tr').length===2")
  setv('context','avoid_cpg');check('motif avoidance is separate from methylation resistance',"document.querySelectorAll('#restriction-rows > tr').length===0")
  click('reset-filters');setv('search','EcoRI','input');evaluate("document.querySelector('#restriction-rows .restriction-name').click();true")
  check('EcoRI recognition span and both strand boundaries appear on duplex',"!!document.querySelector('#restriction-duplex [data-enzyme=EcoRI][data-site=\"10\"]') && !!document.querySelector('#restriction-duplex [data-cut-strand=top][data-cut-boundary=\"10\"]') && !!document.querySelector('#restriction-duplex [data-cut-strand=bottom][data-cut-boundary=\"14\"]')")
  check('bottom bases are aligned plain complements',"[...document.querySelectorAll('#restriction-duplex [data-strand=top]')].map(e=>e.textContent).join('').split('').map(c=>({A:'T',C:'G',G:'C',T:'A'})[c]).join('')===[...document.querySelectorAll('#restriction-duplex [data-strand=bottom]')].map(e=>e.textContent).join('')")
  evaluate("document.querySelector('#restriction-duplex [data-enzyme=EcoRI]').dispatchEvent(new MouseEvent('click',{bubbles:true}));true")
  check('clicking a span selects enzyme and explicit cut coordinates',"document.getElementById('restriction-view-enzyme').value==='EcoRI' && document.getElementById('restriction-selection').textContent.includes('10 / 14')")
  setv('view-width','40');check('zoom changes the base viewport',"document.querySelectorAll('#restriction-duplex [data-strand=top]').length===40")
  click('pan-next');check('pan moves the sequence viewport',"Number(document.getElementById('restriction-view-start').value)>1")
  click('export-view');duplex=evaluate("window.__blobs[window.__blobs.length-1].text()",True)
  import xml.etree.ElementTree as ET
  assert ET.fromstring(duplex).tag=='{http://www.w3.org/2000/svg}svg';checks.append('interactive duplex exports standalone SVG')
  evaluate("document.getElementById('restriction-table').open=true;true")
  check('enzyme details expose correctly placed cuts and a compact vector map',"document.querySelector('.restriction-sites').textContent.includes('top/bottom 10 / 14') && !!document.querySelector('svg.restriction-map')")
  evaluate("[...document.querySelectorAll('.restriction-sites button')].find(b=>b.textContent==='Export map SVG').click();true")
  svg=evaluate("window.__blobs[window.__blobs.length-1].text()",True)
  import xml.etree.ElementTree as ET
  assert ET.fromstring(svg).tag=='{http://www.w3.org/2000/svg}svg';checks.append('exported map is a standalone valid SVG')
  click('export');text=evaluate("window.__blobs[window.__blobs.length-1].text()",True);assert '\nEcoRI\tGAATTC\t1\t0\t' in text;checks.append('summary TSV exports filtered counts')
  click('export-sites');text=evaluate("window.__blobs[window.__blobs.length-1].text()",True);assert '\t10\t14\ttrue\t' in text;checks.append('site TSV exports explicit cut coordinates')
  setv('sequence','TTCAGAA','input');check('editing input clears stale results',"document.getElementById('restriction-results').hidden")
  click('all');submit();ready();check('linear mode excludes origin-spanning EcoRI site',"document.querySelector('#restriction-rows tr').children[2].textContent==='0'")
  setv('topology','circular');submit();ready();check('circular mode finds the origin-spanning site',"document.querySelector('#restriction-rows tr').children[2].textContent==='1'")
  check('circular recognition wraps across the displayed origin',"!!document.querySelector('#restriction-duplex [data-enzyme=EcoRI]') && document.querySelectorAll('#restriction-duplex [data-strand=top]').length===7")
  setv('sequence','TTTTTTTTTTTTTTTGAGACCTTTTTTTTTTTTTTTTTTTTTT','input');setv('topology','linear');setv('search','BsaI','input');submit();ready()
  check('reverse type IIS display keeps cuts outside recognition with correct polarity',"!!document.querySelector('#restriction-duplex [data-enzyme=BsaI]') && [...document.querySelectorAll('#restriction-duplex [data-cut-strand=top]')].some(e=>Number(e.getAttribute('data-cut-boundary'))<15)")
  click('reset-filters');setv('sequence','TGAATTCA','input');setv('custom','OnlyMine GAATTC 1 5 blocked unaffected unknown','input');evaluate("document.getElementById('restriction-custom-only').checked=true");submit();ready()
  check('custom enzymes and their explicit sensitivity states are supported',"document.querySelectorAll('#restriction-rows > tr').length===1 && document.querySelector('#restriction-rows .restriction-name').textContent==='OnlyMine' && document.querySelector('#restriction-rows tr').children[5].textContent==='Blocked'")
  setv('custom','EcoRI GAATTC 1 5','input');submit();wait("document.getElementById('restriction-status').textContent.includes('Duplicate enzyme')");checks.append('custom definitions cannot silently replace catalogue enzymes')
  setv('custom','','input');click('reset-filters');evaluate("document.getElementById('restriction-sequence').value='ACGT'.repeat(25000);document.getElementById('restriction-sequence').dispatchEvent(new Event('input',{bubbles:true}));true");previous=evaluate('window.__workerCount');submit();wait('window.__workerCount>'+str(previous));click('cancel');check('active analysis cancels without retaining results',"document.getElementById('restriction-status').textContent==='Cancelled.' && document.getElementById('restriction-results').hidden && !document.getElementById('restriction-analyze').disabled")
  click('example');submit();ready();checks.append('analysis recovers after cancellation')
  setv('sequence','>one\nACGT\n>two\nACGT','input');submit();wait("document.getElementById('restriction-status').textContent.includes('single FASTA')");checks.append('multiple FASTA records fail rather than concatenate')
  check('browser had no uncaught errors',"window.__errors.length===0")
  report=json.dumps({'checks':checks,'passed':len(checks)},indent=2)+'\n'
  args.report.write_text(report);print(report)
 finally:
  browser.terminate();browser.wait(timeout=10);server.shutdown()
