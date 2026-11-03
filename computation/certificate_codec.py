#!/usr/bin/env python3
"""Losslessly encode the frozen graph-search JSONL certificate format.

The encoded file preserves every original field. Pair order and parity are
specified by the header; candidate and automorphism-call count use two uint16s.
Every encode is followed by a full streaming decode and original-SHA256 check.
Original files are never removed. This codec is not used for finding proofs.
"""
from pathlib import Path
import argparse,hashlib,json,re,struct,sys,time

HERE=Path(__file__).resolve().parent
MAGIC=b'K213B04\n'
RECORD=struct.Struct('<HH')
LENGTH=struct.Struct('<I')
PATTERN=re.compile(rb'\{"n":(\d+),"i":(\d+),"j":(\d+),"parity":(-?1),"candidate":(\d+),"proposals":(\d+),"automorphism_calls":(\d+)\}\n')

def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  for block in iter(lambda:f.read(1<<20),b''):h.update(block)
 return h.hexdigest()

def keys(odd):
 k=len(odd)
 for i in range(k):
  yield i+1,i+1,1
  if not odd[i]:yield i+1,i+1,-1
 for i in range(k):
  for j in range(i+1,k):
   yield i+1,j+1,1
   if not odd[i] and not odd[j]:yield i+1,j+1,-1

def json_line(n,i,j,p,c,a):
 return (f'{{"n":{n},"i":{i},"j":{j},"parity":{p},"candidate":{c},"proposals":{c},"automorphism_calls":{a}}}\n').encode('ascii')

def read_header(f):
 if f.read(len(MAGIC))!=MAGIC:raise ValueError('wrong certificate format/version')
 raw=f.read(LENGTH.size)
 if len(raw)!=LENGTH.size:raise ValueError('missing header length')
 length=LENGTH.unpack(raw)[0]
 if length<1 or length>16*1024*1024:raise ValueError('invalid header size')
 raw=f.read(length)
 if len(raw)!=length:raise ValueError('truncated header')
 h=json.loads(raw)
 odd=h['has_odd_generator'];k=h['cover_count']
 if h['format']!='K213B04' or h['pair_order']!='diagonal_then_strict_upper_triangle':raise ValueError('unsupported header convention')
 if not isinstance(k,int) or not 1<=k<=1000000 or len(odd)!=k or any(type(o) is not bool for o in odd):raise ValueError('invalid envelope flags')
 e=sum(not o for o in odd)
 if h['witness_count']!=k*(k+1)//2+e*(e+1)//2:raise ValueError('wrong witness count')
 if not 1<=h['pool_size']<=65535 or not isinstance(h['degree'],int) or h['degree']<1:raise ValueError('invalid dimensions')
 return h

def stream_decode(path,output=None,mode='verify'):
 """Reconstruct JSONL bytes for hashing; optionally write JSONL or queries."""
 path=Path(path);target=Path(output) if output is not None else None
 if target and target.exists():raise FileExistsError(f'refusing to replace {target}')
 temporary=target.with_name(target.name+'.codec_tmp') if target else None
 digest=hashlib.sha256();total=0;count=0;start=time.monotonic()
 out=temporary.open('wb') if temporary else None
 try:
  with path.open('rb') as f:
   h=read_header(f);n=h['degree'];pair_iterator=iter(keys(h['has_odd_generator']))
   expected=h['witness_count']
   if path.stat().st_size!=f.tell()+expected*RECORD.size:raise ValueError('encoded byte count mismatch')
   if out and mode=='queries':out.write(f'{expected}\n'.encode('ascii'))
   while True:
    raw=f.read(1<<20)
    if not raw:break
    if len(raw)%RECORD.size:raise ValueError('truncated body record')
    for c,a in RECORD.iter_unpack(raw):
     if not 1<=a<=c<=h['pool_size']:raise ValueError('invalid encoded counters')
     try:i,j,p=next(pair_iterator)
     except StopIteration as exc:raise ValueError('extra record') from exc
     line=json_line(n,i,j,p,c,a)
     digest.update(line);total+=len(line);count+=1
     if out:
      out.write(line if mode=='decode' else f'{i} {j} {p} {c}\n'.encode('ascii'))
   if count!=expected:raise ValueError('missing record')
   try:next(pair_iterator);raise ValueError('pair iterator not exhausted')
   except StopIteration:pass
   actual=digest.hexdigest()
   if actual!=h['jsonl_sha256'] or total!=h['jsonl_bytes']:raise ValueError('original JSONL reconstruction hash/size mismatch')
  if out:out.close();out=None
  if target:temporary.replace(target)
  return {'degree':n,'records':count,'reconstructed_jsonl_bytes':total,'reconstructed_jsonl_sha256':actual,'matches_original_sha256':True,'seconds':round(time.monotonic()-start,3)}
 finally:
  if out:out.close()
  if temporary and temporary.exists():temporary.unlink()

def encode(n):
 start=time.monotonic();prefix=HERE/f'graph_{n}'
 audit=json.loads((HERE/f'graph_{n}_replay_audit.json').read_text())
 if not (audit.get('nauty_search_pass') and audit.get('traces_full_replay_pass') and audit.get('coverage_exact')):raise ValueError('degree not fully certified')
 ga=json.loads((HERE/f'graph_{n}_generator_audit.json').read_text())
 source=HERE/f'graph_{n}_certificates.jsonl';encoded=HERE/f'graph_{n}_certificates.kc4'
 if encoded.exists():raise FileExistsError(f'refusing to replace {encoded}')
 odd=ga['has_odd_generator'];k=len(odd)
 h={'format':'K213B04','degree':n,'cover_count':k,'has_odd_generator':odd,'witness_count':audit['recorded_witness_count'],'pool_size':audit['shared_permutation_pool_size'],'pair_order':'diagonal_then_strict_upper_triangle','body_layout':'little-endian uint16 candidate, little-endian uint16 automorphism_calls','jsonl_file':source.name,'jsonl_bytes':source.stat().st_size,'jsonl_sha256':audit['certificate_sha256'],'pool_sha256':audit['permutation_pool_sha256'],'generator_audit_sha256':sha(HERE/f'graph_{n}_generator_audit.json')}
 raw=json.dumps(h,separators=(',',':')).encode('utf-8')
 temporary=encoded.with_name(encoded.name+'.codec_tmp')
 digest=hashlib.sha256();count=0;buffer=bytearray();iterator=iter(keys(odd))
 try:
  with source.open('rb') as f,temporary.open('wb') as out:
   out.write(MAGIC);out.write(LENGTH.pack(len(raw)));out.write(raw)
   for line in f:
    match=PATTERN.fullmatch(line)
    if match is None:raise ValueError(f'noncanonical source record {count+1}')
    nn,i,j,p,c,proposals,a=map(int,match.groups())
    if nn!=n or (i,j,p)!=next(iterator) or proposals!=c or not 1<=a<=c<=h['pool_size']<=65535:raise ValueError(f'source record ordering/counter mismatch {count+1}')
    if line!=json_line(n,i,j,p,c,a):raise ValueError('source is not bytewise canonical')
    digest.update(line);buffer+=RECORD.pack(c,a);count+=1
    if len(buffer)>=1<<20:out.write(buffer);buffer.clear()
   if buffer:out.write(buffer)
  if count!=h['witness_count'] or digest.hexdigest()!=h['jsonl_sha256']:raise ValueError('source count/hash mismatch')
  try:next(iterator);raise ValueError('missing source pairs')
  except StopIteration:pass
  temporary.replace(encoded)
  roundtrip=stream_decode(encoded)
  report={'degree':n,'original_file':source.name,'encoded_file':encoded.name,'original_bytes':h['jsonl_bytes'],'encoded_bytes':encoded.stat().st_size,'original_sha256':h['jsonl_sha256'],'encoded_sha256':sha(encoded),'codec_source_sha256':sha(Path(__file__)),'lossless_roundtrip':roundtrip,'all_original_fields_preserved':True,'original_file_retained':True,'total_seconds':round(time.monotonic()-start,3),'scope':'Storage encoding only; the mathematical certificate is unchanged. A decoder reconstructs the exact original JSONL byte stream or the replay queries.'}
  (HERE/f'graph_{n}_certificate_codec_audit.json').write_text(json.dumps(report,indent=2)+'\n')
  print(json.dumps(report),flush=True)
 finally:
  if temporary.exists():temporary.unlink()

if __name__=='__main__':
 parser=argparse.ArgumentParser();subs=parser.add_subparsers(dest='command',required=True)
 p=subs.add_parser('encode');p.add_argument('degrees',nargs='+',type=int)
 for name in ['decode','queries','verify']:
  p=subs.add_parser(name);p.add_argument('encoded_file');p.add_argument('--output')
 args=parser.parse_args()
 if args.command=='encode':
  for degree in args.degrees:encode(degree)
 else:
  if args.command!='verify' and not args.output:parser.error('--output is required for decode/queries')
  if args.command=='verify' and args.output:parser.error('verify does not take --output')
  print(json.dumps(stream_decode(args.encoded_file,args.output,args.command)),flush=True)
