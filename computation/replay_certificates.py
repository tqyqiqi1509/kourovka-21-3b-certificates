#!/usr/bin/env python3
"""Deterministically replay every GAP witness; prove exact pair/parity coverage.

Usage: python3 replay_certificates.py PREFIX N [N ...]
Supported certificate formats: full permutation x; or a candidate index into a
shared even-permutation pool. Search randomness and the core sieve are absent.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import time

HERE=Path(__file__).resolve().parent
GAP=Path(os.environ.get('GAP_BINARY',HERE/'gap_local/usr/lib/x86_64-linux-gnu/gap/gap'))
ROOT=Path(os.environ.get('GAP_ROOT',HERE/'gap_local/usr/share/gap'))

def sha(path:Path)->str:
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()

def assignment(text:str,name:str):
    m=re.search(r'(?m)^'+re.escape(name)+r'\s*:=\s*(.*?);',text,re.S)
    if m is None:raise ValueError(f'missing GAP assignment: {name}')
    literal=re.sub(r'\[\s*(\d+)\s*\.\.\s*(\d+)\s*\]',lambda a:json.dumps(list(range(int(a.group(1)),int(a.group(2))+1))),m.group(1))
    return json.loads(literal)

def sign(a:list[int])->int:
    seen=[False]*len(a);cycles=0
    for i in range(len(a)):
        if not seen[i]:
            cycles+=1;j=i
            while not seen[j]:seen[j]=True;j=a[j]-1
    return -1 if (len(a)-cycles)%2 else 1

def valid_permutation(a:list[int],n:int)->bool:
    return len(a)==n and sorted(a)==list(range(1,n+1))

def replay(prefix:str,n:int)->None:
    start_all=time.monotonic()
    cert=HERE/f'{prefix}_{n}_certificates.jsonl'
    groupfile=HERE/f'{prefix}_{n}_groups.g'
    unresolved=HERE/f'{prefix}_{n}_unresolved.jsonl'
    poolfile=HERE/f'{prefix}_{n}_permutation_pool.jsonl'
    if unresolved.exists() and unresolved.stat().st_size:
        raise ValueError(f'unresolved pairs remain: {unresolved}')
    gt=groupfile.read_text()
    assert assignment(gt,'ExplicitDegree')==n
    gens=assignment(gt,'ExplicitGenerators')
    orders=assignment(gt,'ExplicitGroupOrders')
    covers=assignment(gt,'ExplicitCoverWords')
    k=len(gens);assert len(orders)==k==len(covers)
    for gs in gens:
        assert all(valid_permutation(g,n) for g in gs)
    hasodd=[any(sign(g)==-1 for g in gs) for gs in gens]
    e=sum(not b for b in hasodd)
    pair_count=k*(k+1)//2
    expected_count=pair_count+e*(e+1)//2
    flags=bytearray(pair_count)
    pool=[]
    if poolfile.exists():
        with poolfile.open() as f:
            for line in f:
                d=json.loads(line)
                assert d['index']==len(pool)+1
                assert valid_permutation(d['x'],n) and sign(d['x'])==1
                pool.append(d['x'])
    script=HERE/f'{prefix}_{n}_replay_generated.g'
    log=HERE/f'{prefix}_{n}_replay.log'
    count=0
    with script.open('w') as w:
        w.write('SizeScreen([1000,1000]);\nRead('+json.dumps(str(groupfile),ensure_ascii=False)+');\n')
        w.write('ReplayGroups:=List(ExplicitGenerators,gs->Group(List(gs,PermList)));;\n')
        w.write('for z in [1..Length(ReplayGroups)] do if Size(ReplayGroups[z])<>ExplicitGroupOrders[z] then Error("group order mismatch");fi;od;\n')
        w.write('ReplayCount:=0;;ReplayStart:=Runtime();;\n')
        w.write('ReplayIncrement:=function() ReplayCount:=ReplayCount+1; if ReplayCount mod 100000=0 then Print("REPLAY_PROGRESS count=",ReplayCount," ms=",Runtime()-ReplayStart,"\\n");fi;end;\n')
        if pool:
            w.write('ReplayPool:=List('+json.dumps(pool,separators=(',',':'))+',PermList);;\n')
            w.write('if not ForAll(ReplayPool,x->SignPerm(x)=1) then Error("pool parity");fi;\n')
            w.write('ReplayCache:=List(ReplayGroups,g->rec(even:=[],odd:=[]));;\n')
            w.write('VerifyPoolWitness:=function(i,j,parity,c) local x,cache,kx;\n')
            w.write(' x:=ReplayPool[c];if parity=-1 then x:=x*(1,2);fi;\n')
            w.write(' if SignPerm(x)<>parity then Error("wrong parity");fi;\n')
            w.write(' if parity=1 then cache:=ReplayCache[j].even;else cache:=ReplayCache[j].odd;fi;\n')
            w.write(' if c<=8 then if not IsBound(cache[c]) then cache[c]:=ReplayGroups[j]^x;fi;kx:=cache[c];else kx:=ReplayGroups[j]^x;fi;\n')
            w.write(' if Size(Intersection(ReplayGroups[i],kx))<>1 then Error("nontrivial intersection");fi;ReplayIncrement();end;\n')
        w.write('VerifyWitness:=function(i,j,parity,a) local x;\n')
        w.write(' if Length(a)<>ExplicitDegree or Set(a)<>[1..ExplicitDegree] then Error("invalid permutation");fi;\n')
        w.write(' x:=PermList(a);if SignPerm(x)<>parity then Error("wrong parity");fi;\n')
        w.write(' if Size(Intersection(ReplayGroups[i],ReplayGroups[j]^x))<>1 then Error("nontrivial intersection");fi;ReplayIncrement();end;\n')
        with cert.open() as f:
            for line_no,line in enumerate(f,1):
                d=json.loads(line)
                assert d['n']==n
                i,j,p=d['i'],d['j'],d['parity']
                if not (1<=i<=j<=k and p in (-1,1)):
                    raise ValueError(f'invalid pair/parity at {line_no}')
                if p==-1 and (hasodd[i-1] or hasodd[j-1]):
                    raise ValueError(f'unexpected odd record at {line_no}')
                off=k*(i-1)-(i-1)*(i-2)//2+(j-i)
                bit=1 if p==1 else 2
                if flags[off]&bit:raise ValueError(f'duplicate pair/parity at {line_no}')
                flags[off]|=bit
                count+=1
                if 'x' in d:
                    assert valid_permutation(d['x'],n) and sign(d['x'])==p
                    w.write(f'VerifyWitness({i},{j},{p},{json.dumps(d["x"],separators=(",",":"))});\n')
                else:
                    c=d['candidate'];assert 1<=c<=len(pool)
                    w.write(f'VerifyPoolWitness({i},{j},{p},{c});\n')
        if count!=expected_count or any(v==0 for v in flags):
            raise ValueError(f'incomplete coverage: {count}/{expected_count}')
        w.write(f'if ReplayCount<>{expected_count} then Error("count mismatch");fi;\n')
        w.write(f'Print("REPLAY_SUCCESS n={n} covers={k} witnesses=",ReplayCount," ms=",Runtime()-ReplayStart,"\\n");\nQUIT;\n')
    print(f'Prepared exact replay n={n}, groups={k}, witnesses={count}, pool={len(pool)}',flush=True)
    start=time.monotonic()
    with log.open('w') as f:
        p=subprocess.run([str(GAP),'-l',str(ROOT),'-A','-q','-b','-T','-m','128m','-o','1024m',str(script)],stdout=f,stderr=subprocess.STDOUT)
    text=log.read_text()
    ok=p.returncode==0 and f'REPLAY_SUCCESS n={n} covers={k} witnesses={count} ' in text and 'Error' not in text
    report={
        'degree':n,'cover_count':k,'unordered_pair_count':pair_count,
        'all_even_cover_count':e,'expected_witness_count':expected_count,
        'recorded_witness_count':count,'coverage_exact':True,
        'permutations_and_parities_checked':True,'deterministic_gap_replay_pass':ok,
        'wall_seconds_gap':round(time.monotonic()-start,3),
        'wall_seconds_total':round(time.monotonic()-start_all,3),
        'certificate_sha256':sha(cert),'explicit_groups_sha256':sha(groupfile),
        'replay_log_sha256':sha(log),'replay_program_sha256':sha(Path(__file__)),
        'randomness_used_in_replay':False,'core_sieve_used_in_replay':False,
        'shared_permutation_pool_size':len(pool),'gap_returncode':p.returncode,
    }
    if pool:report['permutation_pool_sha256']=sha(poolfile)
    (HERE/f'{prefix}_{n}_replay_audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report),flush=True)
    if not ok:raise RuntimeError('GAP replay failed: '+str(log))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('prefix');p.add_argument('degrees',type=int,nargs='+')
    a=p.parse_args()
    for n in a.degrees:replay(a.prefix,n)

