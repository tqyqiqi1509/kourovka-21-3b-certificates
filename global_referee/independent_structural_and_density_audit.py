#!/usr/bin/env python3
"""Fresh finite-input audit. Does not import production enumeration/count code.

The full pruned-cover counts are computed with independence polynomials on
the literal-prefix conflict forest, not by enumerating covers. Prime-order
semisimple classes are summed over partitions of multiplicities, not by the
production generating-function powers. No permutation intersections are rerun.
"""
from collections import Counter, defaultdict
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from math import comb, factorial, isqrt, prod
from pathlib import Path
import hashlib
import json
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def primes(bound):
    sieve = [True] * (bound + 1)
    sieve[:2] = [False, False]
    for p in range(2, isqrt(bound) + 1):
        if sieve[p]:
            for k in range(p*p, bound+1, p):
                sieve[k] = False
    return [p for p, yes in enumerate(sieve) if yes]

def powers(bound):
    ans = {}
    for p in primes(bound):
        for d in range(1, bound.bit_length()):
            q = p**d
    
            if q > bound:
                break
            ans[q] = (p, d)
    return ans

PP = powers(400)

@lru_cache(None)
def gl(d, p):
    return prod(p**d - p**i for i in range(d))

@lru_cache(None)
def partition(n, ceiling=None):
    if n == 0:
        return ((),)
    ceiling = min(n, n if ceiling is None else ceiling)
    return tuple((a,) + b for a in range(1, ceiling+1)
                 for b in partition(n-a, a))

@lru_cache(None)
def natural_words(n):
    """All ordered normalized factorizations, via the last factor."""
    if n == 1:
        return ((),)
    ans = []
    for last in sorted(PP):
        if last > n:
            break
        if n % last:
            continue
        for prefix in natural_words(n//last):
            if prefix and prefix[-1] == last == 2:
                continue
            ans.append(prefix+(last,))
    assert len(set(ans)) == len(ans)
    return tuple(sorted(ans))

def words_through(bound):
    return [w for n in range(1, bound+1) for w in natural_words(n)]

def mul_poly(a, b, bound):
    out = [0]*(min(bound, len(a)+len(b)-2)+1)
    for i, ai in enumerate(a):
        if not ai:
            continue
        for j, bj in enumerate(b[:len(out)-i]):
            if bj:
                out[i+j] += ai*bj
    return out

def add_poly(a, b):
    out = [0]*max(len(a), len(b))
    for i, ai in enumerate(a):
        out[i] += ai
    for i, bi in enumerate(b):
        out[i] += bi
    return out

def pruned_independence_polynomial(words, bound, minimum_degree=1):
    """Independent sets of the conflict forest w--(w,2/3), weighted z^degree."""
    present = {w for w in words if minimum_degree <= prod(w) <= bound}
    children = defaultdict(list)
    roots = []
    for w in present:
        parent = w[:-1] if w and w[-1] in (2,3) else None
        if parent in present:
            children[parent].append(w)
        else:
            roots.append(w)
    @lru_cache(None)
    def node(w):
        excluded, included = [1], [0]*prod(w)+[1]
        for child in children[w]:
            c0, c1 = node(child)
            excluded = mul_poly(excluded, add_poly(c0,c1), bound)
            included = mul_poly(included, c0, bound)
        return excluded, included
    result = [1]
    for root in roots:
        r0,r1 = node(root)
        result = mul_poly(result, add_poly(r0,r1), bound)
    return result+[0]*(bound+1-len(result))

def literal_conflict(u,v):
    return (u+(2,) == v or u+(3,) == v or
            v+(2,) == u or v+(3,) == u)

def residual_pruned(words,n):
    eligible = [i for i,w in enumerate(words) if (n+3)//4 < prod(w) <= n]
    by_degree = defaultdict(list)
    for i in eligible:
        by_degree[prod(words[i])].append(i)
    result = []
    for i in eligible:
        di = prod(words[i])
        if di == n:
            result.append([i])
        for j in eligible:
            if j <= i:
                continue
            dj = prod(words[j])
            if di+dj > n:
                break
            if literal_conflict(words[i],words[j]):
                continue
            if di+dj == n:
                result.append([i,j])
            for k in by_degree[n-di-dj]:
                if k > j and not literal_conflict(words[i],words[k]) and not literal_conflict(words[j],words[k]):
                    result.append([i,j,k])
    return sorted(result)

def audit_cover_partition():
    template = ROOT/'computation/lifting/lifted_residual_covers_through_63.json'
    data = json.loads(template.read_text())
    words = words_through(63)
    assert data['words'] == [list(w) for w in words]
    polys = {b:pruned_independence_polynomial(words,63,b) for b in range(1,18)}
    rows = []
    for n in range(40,64):
        auditpath = ROOT/f'computation/lifting/inheritance_partition_{n}.json'
        expected = json.loads(auditpath.read_text())
        residual = residual_pruned(words,n)
        assert residual == data['covers_0based'][str(n)]
        bound = (n+3)//4
        deps = {str(n-b):polys[b][n]-polys[b+1][n] for b in range(1,bound+1)}
        deps = {d:c for d,c in deps.items() if c}
        assert expected['full_cover_count'] == polys[1][n]
        assert expected['residual_cover_count'] == len(residual) == polys[bound+1][n]
        assert expected['inherited_by_dependency'] == deps
        assert expected['inherited_cover_count'] == sum(deps.values())
        rows.append({'n':n,'full_covers':polys[1][n], 'residual_covers':len(residual),
                     'residual_pairs':len(residual)*(len(residual)+1)//2,
                     'every_dependency_count_matched':True,
                     'partition_record_sha256':sha(auditpath)})
    old = json.loads((ROOT/'structural_reduction/affine_covers_pruned_through_40.json').read_text())
    assert all(len(covers) == polys[1][int(n)] for n,covers in old['covers_0based'].items())
    return {'method':'Conflict-forest independence polynomials, plus at-most-three-factor residual enumeration',
            'all_pass':True,'rows':rows,'total_residual_pairs':sum(r['residual_pairs'] for r in rows),
            'template_sha256':sha(template)}

def degree_patterns(n):
    b=(n+3)//4
    all_patterns=[(n,)]
    all_patterns += [(i,n-i) for i in range(b+1,n//2+1)]
    all_patterns += [(i,j,n-i-j) for i in range(b+1,n//3+1)
                     for j in range(i,(n-i)//2+1)]
    def separated(ds):
        return any(21<=sum(d for k,d in enumerate(ds) if mask&(1<<k))<=n-21
                   for mask in range(1,2**len(ds)-1))
    return [ds for ds in all_patterns if not separated(ds)]

def audit_extended_residuals():
    path=ROOT/'prime_order_referee/lifting_extended_unpruned_templates.json'
    data=json.loads(path.read_text())
    words=words_through(127);indices={w:i for i,w in enumerate(words)}
    assert data['words']==[list(w) for w in words]
    rows=[]
    for n in range(64,128):
        pats=degree_patterns(n)
        assert [list(p) for p in pats]==data['degree_patterns'][str(n)]
        result=set()
        for ds in pats:
            for ww in product(*(natural_words(d) for d in ds)):
                result.add(tuple(sorted(indices[w] for w in ww)))
        expected={tuple(c) for c in data['covers_0based'][str(n)]}
        assert result==expected and len(expected)==len(data['covers_0based'][str(n)])
        k=len(result)
        rows.append({'n':n,'degree_patterns':pats,'residual_covers':k,'pairs':k*(k+1)//2})
    return {'all_pass':True,'scope':'Structural inputs only, not intersection replay or induction closure',
            'rows':rows,'residual_covers':sum(r['residual_covers'] for r in rows),
            'residual_pairs':sum(r['pairs'] for r in rows),'template_sha256':sha(path)}

def affine_prime_counts(characteristic,d,p):
    """Partition-of-multiplicities class sum, distinct from production powers."""
    q=characteristic**d; ans=defaultdict(int)
    if characteristic==p:
        for lam in partition(d,p):
            multiplicities=Counter(lam)
            centralizer=F(p**sum(sum(x>=j for x in lam)**2 for j in range(1,max(lam)+1)))
            for multiplicity in multiplicities.values():
                for k in range(1,multiplicity+1):
                    centralizer*=1-F(1,p**k)
            assert centralizer.denominator==1 and gl(d,p)%centralizer.numerator==0
            classes=gl(d,p)//centralizer.numerator
            fixed_translation_count=p**(d-len(lam))
            allowed_translation_count=p**(d-multiplicities[p])
            support=q-p**len(lam)
            if support:
                ans[support]+=classes*fixed_translation_count
            if allowed_translation_count>fixed_translation_count:
                ans[q]+=classes*(allowed_translation_count-fixed_translation_count)
    else:
        f=next(k for k in range(1,p) if pow(characteristic,k,p)==1)
        h=(p-1)//f
        for total in range(1,d//f+1):
            m0=d-f*total
            support=q-characteristic**m0
            for lam in partition(total):
                if len(lam)>h:
                    continue
                ways=F(prod(range(h-len(lam)+1,h+1)),prod(factorial(c) for c in Counter(lam).values()))
                centralizer=gl(m0,characteristic)*prod(gl(m,characteristic**f) for m in lam)
                count=ways*F(gl(d,characteristic)*characteristic**(d-m0),centralizer)
                assert count.denominator==1
                ans[support]+=count.numerator
    assert all(s>0 and s%p==0 and c>0 for s,c in ans.items())
    return dict(ans)

def word_order(w):
    value=1
    for q in w:
        p,d=PP[q]
        value=value**q*q*gl(d,p)
    return value

def audit_density_table():
    path=ROOT/'prime_order_referee/prime_atom_densities_200.json'
    data=json.loads(path.read_text())
    ps=primes(200)
    distributions={}
    compared=0
    for q,(char,d) in PP.items():
        distributions[q]={}
        for p in ps:
            if p>q:
                break
            dist=affine_prime_counts(char,d,p)
            if dist:
                distributions[q][p]=dist
                compared+=len(dist)
        expected={int(p):{int(s):c for s,c in dist.items()}
                  for p,dist in data['affine_prime_supports'][str(q)].items()}
        assert distributions[q]==expected,(q,distributions[q],expected)
    orders={m:max(word_order(w) for w in natural_words(m)) for m in range(1,101)}
    assert orders=={int(m):v for m,v in data['child_orders'].items()}
    densities={p:{} for p in ps}
    for q,rows in distributions.items():
        for p,dist in rows.items():
            for j,count in dist.items():
                for m in range(1,200//j+1):
                    if m==1 and q in (2,3,4):
                        continue
                    a=j*m
                    v=F(count*orders[m]**(j-j//p),q*m)
                    densities[p][a]=max(v,densities[p].get(a,F(0)))
    expected={int(p):{int(a):F(v['numerator'],v['denominator']) for a,v in rows.items()}
              for p,rows in data['noncore_prime_atom_densities'].items()}
    assert densities==expected
    inv_comparisons=0
    for q,(char,d) in PP.items():
        if q>=256:
            continue
        coefficients={0:1, **{s//2:c for s,c in distributions[q][2].items()}}
        for k in range(q+1):
            left=sum(c*coefficients.get(k-j,0) for j,c in coefficients.items())
            if k==q:
                left+=3*q*gl(d,char)
            assert left<=comb(q,k)*3**k,(q,k)
            inv_comparisons+=1
    return {'all_pass':True,'affine_prime_support_rows':compared,
            'all_atom_density_rows':sum(len(v) for v in densities.values()),
            'child_maximum_orders_compared':len(orders),'involution_primitive_coefficient_checks':inv_comparisons,
            'semisimple_method':'integer partitions of multiplicity multisets, times falling factorial assignment counts',
            'unipotent_method':'centralizer dimension times product (1-p^-k)',
            'child_order_method':'direct maximum over all normalized ordered factorization words',
            'density_file_sha256':sha(path)}

def main(output=None):
    start=time.monotonic()
    result={'auditor':'global_referee fresh independent implementation','started_unix':time.time()}
    result['finite_40_63_partition']=audit_cover_partition()
    print('COVER_PARTITION_PASS',result['finite_40_63_partition']['total_residual_pairs'],flush=True)
    result['finite_64_127_structure']=audit_extended_residuals()
    print('EXTENDED_STRUCTURE_PASS',result['finite_64_127_structure']['residual_pairs'],flush=True)
    result['N128_density_and_invariant']=audit_density_table()
    print('DENSITIES_AND_INVARIANT_PASS',json.dumps(result['N128_density_and_invariant']),flush=True)
    result['seconds']=time.monotonic()-start
    result['source_sha256']=sha(__file__)
    result['all_checks_passed']=True
    (Path(output) if output is not None else HERE/'independent_structural_and_density_audit.json').write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':
    main()
