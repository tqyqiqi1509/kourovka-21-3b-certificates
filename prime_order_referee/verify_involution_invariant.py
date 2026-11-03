#!/usr/bin/env python3
"""Independent involution formulas and exact primitive polynomial audit.

The odd-characteristic eigenspace formula and characteristic-two rank formula
are independent of the general partition/semisimple code used for the density
table. Every polynomial inequality uses integers, with no numerical roots.
"""
from pathlib import Path
from fractions import Fraction as Q
from math import comb, prod
import json, hashlib
from prime_atom_densities import prime_powers, affine_prime_support, gl

HERE=Path(__file__).resolve().parent

def gaussian(n,r,p):
    x=prod((Q(p**(n-i)-1,p**(r-i)-1) for i in range(r)),start=Q(1))
    assert x.denominator==1
    return x.numerator

def independent_involutions(p,d):
    q=p**d;ans={}
    def add(s,c):
        if c:ans[s]=ans.get(s,0)+c
    if p%2:
        for r in range(1,d+1):
            den=gl(r,p)*gl(d-r,p)
            assert gl(d,p)%den==0
            add(q-p**(d-r),gl(d,p)//den*p**r)
    else:
        add(q,q-1)
        for r in range(1,d//2+1):
            linear=gaussian(d,r,2)*gaussian(d-r,r,2)*gl(r,2)
            add(q-2**(d-r),linear*2**r)
            add(q,linear*(2**(d-r)-2**r))
    return ans

def main():
    rows=[];failures=[];density_failures=[];support_rows=0
    for q,(p,d) in sorted(prime_powers(255).items()):
        inv=independent_involutions(p,d)
        assert inv==affine_prime_support(p,d,2),(q,inv,affine_prime_support(p,d,2))
        support_rows+=len(inv)
        f=[0]*(q+1);f[0]=1
        for s,c in inv.items():
            assert s%2==0;f[s//2]=c
            lhs=Q(c*s,2*q)
            if lhs**6>24**s:density_failures.append({'q':q,'support':s,'count':c,'lhs':str(lhs)})
        comparisons=[]
        for r in range(q+1):
            left=sum(f[k]*f[r-k] for k in range(r+1))+(3*q*gl(d,p) if r==q else 0)
            right=comb(q,r)*3**r
            if left>right:failures.append({'q':q,'coefficient_index':r,'left':str(left),'right':str(right)})
            if left:comparisons.append({'r':r,'left':str(left),'right':str(right),'passed':left<=right})
        rows.append({'q':q,'p':p,'d':d,'support_distribution':inv,
                     'involutions_including_identity':1+sum(inv.values()),
                     'coefficients_checked':q+1,'all_pass':not any(z['q']==q for z in failures),
                     'nonzero_left_coefficients':comparisons})
    assert not failures
    assert density_failures==[{'q':8,'support':4,'count':42,'lhs':'21/2'}]
    result={'primitive_groups_checked':len(rows),'independent_support_rows_matched':support_rows,
            'coefficients_checked':sum(r['coefficients_checked'] for r in rows),'failures':failures,
            'local_density_exceptions':density_failures,
            'independent_formula':'odd eigenspace classes; even square-zero rank and image/kernel flags',
            'verifier_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'rows':rows}
    (HERE/'involution_polynomial_invariant_audit.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='rows'}))

if __name__=='__main__':main()
