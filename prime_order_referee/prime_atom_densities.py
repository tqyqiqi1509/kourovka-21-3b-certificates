#!/usr/bin/env python3
"""Exact affine prime-order support counts and forest atom densities.

No group database or floating point is used. The formulas are justified in
PRIME_ORDER_MECHANISM.md. Child words may use full affine (nonsoluble) groups.
"""
from functools import lru_cache
from fractions import Fraction as Q
from math import prod, isqrt
from pathlib import Path
import argparse, json, sys

sys.set_int_max_str_digits(0)

def primes(B):
    return [p for p in range(2,B+1) if all(p % r for r in range(2,isqrt(p)+1))]

def prime_powers(B):
    out={}
    for p in primes(B):
        q=p; d=1
        while q<=B:
            out[q]=(p,d); q*=p; d+=1
    return out

@lru_cache(None)
def gl(d,p):
    return prod(p**d-p**i for i in range(d))

@lru_cache(None)
def partitions(n,maxpart):
    if n==0:return ((),)
    return tuple((k,)+rest for k in range(min(n,maxpart),0,-1)
                 for rest in partitions(n-k,k))

@lru_cache(None)
def multiplicative_order(a,p):
    assert a%p
    x=a%p;k=1
    while x!=1:
        x=x*a%p;k+=1
    return k

def conv(a,b,L):
    out=[Q(0)]*(L+1)
    for i,x in enumerate(a):
        for j,y in enumerate(b[:L+1-i]):out[i+j]+=x*y
    return out

@lru_cache(None)
def affine_prime_support(char,d,ell):
    """Return support->count for elements of exact prime order ell in AGL(d,char)."""
    q=char**d; out={}; GL=gl(d,char)
    def add(s,c):
        assert isinstance(c,int) and c>=0
        if c:
            assert s>0 and s%ell==0
            out[s]=out.get(s,0)+c
    if ell==char:
        for lam in partitions(d,ell):
            mult={j:lam.count(j) for j in set(lam)}
            columns=[sum(j>=i for j in lam) for i in range(1,lam[0]+1)]
            dim=sum(x*x for x in columns)
            rad=dim-sum(v*v for v in mult.values())
            cent=char**rad*prod(gl(v,char) for v in mult.values())
            assert GL%cent==0
            cls=GL//cent; blocks=len(lam); full=mult.get(ell,0)
            fix_b=char**(d-blocks)
            all_b=char**(d-full)
            if lam[0]>1:add(q-char**blocks,cls*fix_b)
            add(q,cls*(all_b-fix_b))
    else:
        f=multiplicative_order(char,ell)
        if f>d:return out
        k=(ell-1)//f; L=d//f
        a=[Q(1,gl(m,char**f)) for m in range(L+1)]
        poly=[Q(1)]+[Q(0)]*L
        for _ in range(k):poly=conv(poly,a,L)
        for h in range(1,L+1):
            m0=d-f*h
            count=Q(GL,gl(m0,char))*poly[h]*char**(d-m0)
            assert count.denominator==1
            add(q-char**m0,count.numerator)
    return out

def all_densities(B):
    pp=prime_powers(2*B); ps=primes(B)
    orders={q:q*gl(d,p) for q,(p,d) in pp.items()}
    O={1:1}; OW={1:()}
    for m in range(2,B//2+1):
        O[m],OW[m]=max((orders[q]*O[m//q]**q,OW[m//q]+(q,))
                       for q in pp if m%q==0)
    rows={ell:{} for ell in ps}; supports={}
    for q,(char,d) in sorted(pp.items()):
        supports[q]={}
        for ell in ps:
            if ell>q:break
            dist=affine_prime_support(char,d,ell)
            if dist:supports[q][ell]=dist
            for j,count in dist.items():
                if j>B:continue
                for m in range(1,B//j+1):
                    a=j*m
                    # Core cycles are encoded separately as individual cycles.
                    if m==1 and q in (2,3,4):continue
                    value=Q(count*O[m]**(j-j//ell),q*m)
                    rec={'a':a,'ell':ell,'numerator':value.numerator,
                         'denominator':value.denominator,'density':str(value),
                         'q':q,'char':char,'dimension':d,'j':j,'m':m,
                         'primitive_count':count,'child_order':O[m],
                         'child_word':OW[m]}
                    if a not in rows[ell] or value>Q(rows[ell][a]['numerator'],rows[ell][a]['denominator']):
                        rows[ell][a]=rec
    return {'max_weight':B,'noncore_prime_atom_densities':rows,
            'affine_prime_supports':supports,'child_orders':O,'child_words':OW}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--max-weight',type=int,default=200)
    args=ap.parse_args();x=all_densities(args.max_weight)
    p=Path(__file__).parent/f'prime_atom_densities_{args.max_weight}.json'
    p.write_text(json.dumps(x,indent=2)+'\n')
    for ell in (2,3,5,7):
        print(ell,[(a,v['density'],v['q'],v['m']) for a,v in sorted(x['noncore_prime_atom_densities'][ell].items()) if a<=40])
    print('written',p)

if __name__=='__main__':main()
