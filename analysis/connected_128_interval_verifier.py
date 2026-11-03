#!/usr/bin/env python3
"""Independent interval-arithmetic replay of the connected-involution model.

Uses unnormalized, actual atom weights at n=128, instead of the normalized
Fraction implementation. Every addition, product, division and logarithmic
coefficient subtraction has lower and upper endpoints. No imported certificate
coefficient is used. Exact Fraction is used only for input constants/tail sums.
"""
from decimal import Decimal as D, Context, ROUND_FLOOR, ROUND_CEILING
from fractions import Fraction as Q
from math import factorial
from pathlib import Path
import hashlib, json, time

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
LOW=Context(prec=90,rounding=ROUND_FLOOR)
HIGH=Context(prec=90,rounding=ROUND_CEILING)

class I:
    __slots__=('lo','hi')
    def __init__(self,lo=0,hi=None):
        self.lo=D(lo); self.hi=D(lo if hi is None else hi)
    @staticmethod
    def rational(q):
        q=Q(q)
        return I(LOW.divide(D(q.numerator),D(q.denominator)),
                 HIGH.divide(D(q.numerator),D(q.denominator)))
    def __add__(self,other):
        if not isinstance(other,I):other=I(other)
        return I(LOW.add(self.lo,other.lo),HIGH.add(self.hi,other.hi))
    __radd__=__add__
    def __mul__(self,other):
        if not isinstance(other,I):other=I(other)
        assert self.lo>=0 and other.lo>=0
        return I(LOW.multiply(self.lo,other.lo),HIGH.multiply(self.hi,other.hi))
    __rmul__=__mul__
    def __truediv__(self,other):
        if not isinstance(other,I):other=I(other)
        assert self.lo>=0 and other.lo>0
        return I(LOW.divide(self.lo,other.hi),HIGH.divide(self.hi,other.lo))
    def positive_difference(self,other):
        lo=LOW.subtract(self.lo,other.hi)
        hi=HIGH.subtract(self.hi,other.lo)
        assert hi>=0
        return I(max(D(0),lo),hi)
    def record(self):return {'lower':str(self.lo),'upper':str(self.hi)}

def exponential_coefficients(weights,limit):
    # From A'=f'A, with unnormalized atom counts as coefficients of f.
    a=[I(1)]+[I() for _ in range(limit)]
    for j in range(1,limit+1):
        a[j]=sum((a[j-k]*weights[k]*k for k in range(1,j+1)),I())/j
    return a

def exact_record(q):
    q=Q(q)
    return {'numerator':str(q.numerator),'denominator':str(q.denominator),
            'interval':I.rational(q).record()}

def main():
    start=time.monotonic()
    n=128; limit=24; bins=30; cap=Q(2048); tau=Q(39,25); rho=Q(5,4)
    raw=ROOT/'prime_order_referee'/'prime_atom_densities_200.json'
    data=json.loads(raw.read_text())
    densities={int(s)//2:Q(v['numerator'],v['denominator'])
               for s,v in data['noncore_prime_atom_densities']['2'].items()}
    fac=[factorial(k) for k in range(201)]
    # Independent exact checks of LLL, conditional inflation, and inverse core cap.
    q=tau/Q(n*(n-1))
    penalty=36*n*q/(1-q)
    lt=2*sum(((tau-1)/(tau+1))**(2*k+1)/Q(2*k+1) for k in range(2))
    lr=2*sum(((rho-1)/(rho+1))**(2*k+1)/Q(2*k+1) for k in range(2))
    lc=11*2*(Q(1,3)+Q(1,3)**3/3)
    inverse=Q(9,2)*n*n*q/(1-q)
    assert penalty<lt and penalty/2<lr and inverse<lc
    first=Q(18**25*fac[25],fac[50])
    invol_tail=cap*first/(1-Q(9,51))
    a=Q(3,5); k=200
    odd_tail=cap*a**k*(Q(k*k)/(1-a)+2*k*a/(1-a)**2+a*(1+a)/(1-a)**3)
    assert Q(729*8,200**2)<a**3 and odd_tail<Q(1,10**36)
    # Independently rebuild every odd-prime finite row with raw atom weights.
    odd_total=I(); odd_count=0
    for p,rows in sorted(data['noncore_prime_atom_densities'].items(),key=lambda t:int(t[0])):
        p=int(p)
        if p==2 or p>=200:continue
        for r in range(2 if p==3 else 1,200//p+1):
            s=p*r
            if s>=200:continue
            m=max(n,s)
            w=[I() for _ in range(r+1)]
            for support,v in rows.items():
                j=int(support)//p
                if j<=r:w[j]=I.rational(m*Q(v['numerator'],v['denominator']))
            if p==3:w[1]=I(2*m)
            count=exponential_coefficients(w,r)[r]
            factor=Q(p**r*fac[r],(p-1)*(fac[m]//fac[m-s]))*min(rho**s,cap)
            odd_total+=count*count*I.rational(factor)
            odd_count+=1
    assert odd_count==267 and odd_total.hi<D('0.013')
    common_rest=odd_total+I.rational(invol_tail+odd_tail)

    H=[]; core=[]
    for i in range(bins):
        left=Q(i,20); right=Q(i+1,20)
        weights=[I() for _ in range(limit+1)]
        for j,d in densities.items():
            if j<=limit:weights[j]=I.rational(n*d)
        weights[2]=I.rational(n*(Q(21,4)-Q(7,2)*left))
        H.append(exponential_coefficients(weights,limit))
        core.append([I.rational((n*right)**j/Q(fac[j])) for j in range(limit+1)])
    factors=[I()]+[I.rational(Q(2**r,fac[n]//fac[n-2*r])*min(rho**(2*r),cap))
                    for r in range(1,limit+1)]
    expected={}
    for line in (ROOT/'prime_order_referee'/'connected_128_density_boxes.jsonl').read_text().splitlines():
        row=json.loads(line); rec=row['connected_involution_upper']
        expected[row['i'],row['j']]=Q(int(rec['numerator']),int(rec['denominator']))
    output=HERE/'connected_128_interval_boxes.jsonl'
    worst=D(0); worst_pair=None; empty=None; rows_checked=0
    with output.open('w') as out:
        for i in range(bins):
            for j in range(i,bins):
                all_graphs=[I(1)]+[I() for _ in range(limit)]
                connected=[I() for _ in range(limit+1)]
                total=I()
                for r in range(1,limit+1):
                    terms=[]
                    # u,v are the noncore port counts, independently of the
                    # other verifier's core-count indexing.
                    for u in range(r+1):
                        for v in range(r-u,r+1):
                            if H[i][u].hi==0 or H[j][v].hi==0:continue
                            terms.append(H[i][u]*H[j][v]*core[i][r-u]*core[j][r-v]
                                         *fac[u]*fac[v]/fac[u+v-r])
                    all_graphs[r]=sum(terms,I())
                    disconnected=sum(((r-k)*all_graphs[k]*connected[r-k]
                                      for k in range(1,r)),I())
                    connected[r]=(all_graphs[r]*r).positive_difference(disconnected)/r
                    total+=connected[r]*factors[r]
                exact=expected[i,j]
                exact_interval=I.rational(exact)
                assert total.lo<=exact_interval.lo<=exact_interval.hi<=total.hi
                complete=total+common_rest
                assert complete.hi<D('0.999')
                if complete.hi>worst:worst=complete.hi;worst_pair=[i,j]
                if i==j==0:
                    empty=complete.record()
                    assert complete.hi<D('0.5')
                out.write(json.dumps({'i':i,'j':j,'involution':total.record(),
                                      'complete':complete.record(),'exact_fraction_enclosed':True})+'\n')
                rows_checked+=1
    assert rows_checked==465
    summary={'method':'independent unnormalized polynomial model, 90-digit two-sided Decimal intervals',
             'N':n,'density_box_pairs':rows_checked,'odd_prime_rows':odd_count,
             'all_fraction_results_enclosed':True,'worst_complete_upper':str(worst),
             'worst_pair':worst_pair,'empty_core_interval':empty,
             'odd_finite_interval':odd_total.record(),'involution_tail':exact_record(invol_tail),
             'odd_tail':exact_record(odd_tail),'lll_penalty':exact_record(penalty),
             'inverse_core_exponent':exact_record(inverse),
             'all_complete_less_than':'999/1000','empty_core_less_than':'1/2',
             'density_file_sha256':hashlib.sha256(raw.read_bytes()).hexdigest(),
             'interval_boxes_sha256':hashlib.sha256(output.read_bytes()).hexdigest(),
             'verifier_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
             'seconds':time.monotonic()-start}
    (HERE/'connected_128_interval_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2),flush=True)

if __name__=='__main__':main()
