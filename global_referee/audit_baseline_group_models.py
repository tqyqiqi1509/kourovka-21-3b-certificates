#!/usr/bin/env python3
"""Audit actual finite group inputs against independently recognized forests.

Membership is checked directly by affine coordinate equations and the
imprimitive block recursion. This does not call the graph encoder, GAP, or
the permutation-intersection engines. The exact group orders themselves were
checked by the preserved GAP construction/replay; this audit checks their
theoretical values and the generators' membership in the named envelopes.
"""
from functools import lru_cache
from math import prod
from pathlib import Path
import json
import re
import time
from independent_structural_and_density_audit import (
    ROOT, HERE, PP, sha, words_through, word_order,
    pruned_independence_polynomial, literal_conflict, mul_poly
)

COORDINATES=json.loads((ROOT/'computation/affine_coordinates.json').read_text())

def assignment(text,name):
    raw=re.search(r'^'+name+r'\s*:=\s*(.*?);',text,re.M|re.S)
    assert raw is not None,name
    value=re.sub(r'\[\s*(\d+)\s*\.\.\s*(\d+)\s*\]',
                 lambda m:json.dumps(list(range(int(m[1]),int(m[2])+1))),raw[1])
    return json.loads(value)

@lru_cache(None)
def affine_member(p):
    q=len(p)
    if q<=4:
        return sorted(p)==list(range(q))
    char,d=PP[q]
    coords=tuple(map(tuple,COORDINATES[str(q)]))
    zero=coords.index((0,)*d)
    b=coords[p[zero]]
    images=[]
    for j in range(d):
        e=tuple(int(i==j) for i in range(d))
        e_index=coords.index(e)
        images.append(tuple((coords[p[e_index]][i]-b[i])%char for i in range(d)))
    return all(tuple((b[i]+sum(v[j]*images[j][i] for j in range(d)))%char
                     for i in range(d))==coords[p[k]] for k,v in enumerate(coords))

@lru_cache(None)
def word_member(w,p):
    if tuple(range(len(p)))==p:
        return True
    if not w:
        return False
    q=w[-1];m=prod(w[:-1])
    if len(p)!=m*q or sorted(p)!=list(range(m*q)):
        return False
    top=[]
    for block in range(q):
        dest={p[block*m+j]//m for j in range(m)}
        if len(dest)!=1:
            return False
        target=next(iter(dest));top.append(target)
        section=tuple(p[block*m+j]-target*m for j in range(m))
        if not word_member(w[:-1],section):
            return False
    return affine_member(tuple(top))

def sign(p):
    seen=set();cycles=0
    for i in range(len(p)):
        if i in seen:
            continue
        cycles+=1;j=i
        while j not in seen:
            seen.add(j);j=p[j]
    return -1 if (len(p)-cycles)%2 else 1

def main(output=None):
    start=time.monotonic();words=words_through(40)
    pruned=pruned_independence_polynomial(words,40)
    unpruned=[1]
    for w in words:
        unpruned=mul_poly(unpruned,[1]+[0]*(prod(w)-1)+[1],40)
    results=[]
    for n in range(21,40):
        prefix='affine' if n<=22 else 'affine_pruned' if n<=26 else 'graph'
        path=ROOT/f'computation/{prefix}_{n}_groups.g'
        text=path.read_text();assert assignment(text,'ExplicitDegree')==n
        covers=assignment(text,'ExplicitCoverWords')
        orders=assignment(text,'ExplicitGroupOrders')
        generators=assignment(text,'ExplicitGenerators')
        k=len(covers)
        assert len(orders)==len(generators)==k==(unpruned if n<=22 else pruned)[n]
        coverkeys=[];all_even=0;checks=0
        for c,order,gens in zip(covers,orders,generators):
            assert c==sorted(set(c)) and all(1<=i<=len(words) for i in c)
            ww=[words[i-1] for i in c]
            assert sum(map(prod,ww))==n
            if n>=23:
                assert all(not literal_conflict(u,v) for i,u in enumerate(ww) for v in ww[i+1:])
            assert order==prod(map(word_order,ww))
            coverkeys.append(tuple(c));has_odd=False
            for perm in gens:
                assert len(perm)==n and sorted(perm)==list(range(1,n+1))
                p=tuple(x-1 for x in perm)
                has_odd|=sign(p)==-1
                offset=0
                for w in ww:
                    d=prod(w)
                    section=tuple(p[offset+j]-offset for j in range(d))
                    assert word_member(w,section),(n,c,w,perm)
                    offset+=d
                checks+=1
            structural_odd=any(w and (PP[w[0]][0]%2 or w[0]<=4) for w in ww)
            assert bool(structural_odd)==has_odd
            all_even+=not has_odd
        assert len(set(coverkeys))==k
        auditpath=ROOT/f'computation/{prefix}_{n}_replay_audit.json'
        bound=False
        if auditpath.exists():
            audit=json.loads(auditpath.read_text())
            assert audit['explicit_groups_sha256']==sha(path)
            assert audit['cover_count']==k and audit['all_even_cover_count']==all_even
            bound=True
        row={'n':n,'complete_family_cover_count':k,'all_even_covers':all_even,
             'actual_generators_directly_recognized':checks,'every_order_formula_matches':True,
             'every_cover_unique_and_admissible':True,'cardinality_matches_independent_polynomial':True,
             'current_completed_replay_audit_hash_bound':bound,'group_file_sha256':sha(path)}
        results.append(row)
        print('BASELINE_MODEL_PASS',n,k,checks,'completed_replay_bound',bound,flush=True)
    report={'all_pass':True,'rows':results,'actual_generators_recognized':sum(r['actual_generators_directly_recognized'] for r in results),
            'method':'Direct affine equations and wreath block recursion; family cardinalities by independent-set polynomials',
            'no_group_intersections_rerun':True,'group_order_limit':'Actual group-order computations are the retained GAP construction/replay; this file checks generator membership and independently computes the expected orders.',
            'coordinates_sha256':sha(ROOT/'computation/affine_coordinates.json'),
            'source_sha256':sha(__file__),'independent_word_module_sha256':sha(HERE/'independent_structural_and_density_audit.py'),
            'seconds':time.monotonic()-start}
    (Path(output) if output is not None else HERE/'baseline_group_models_audit.json').write_text(json.dumps(report,indent=2)+'\n')

if __name__=='__main__':
    main()
