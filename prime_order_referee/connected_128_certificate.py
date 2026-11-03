#!/usr/bin/env python3
"""Exact rational certificate for the candidate uniform cutoff n>=128.

Every probability bound, logarithm comparison, formal logarithm coefficient,
and final comparison is computed with fractions.Fraction and integer
factorials. Decimal is used only to print readable (upward) summaries.
The mathematical model and monotonicity lemmas must accompany this file.
"""
from fractions import Fraction as Q
from decimal import Decimal as D, localcontext, ROUND_CEILING
from math import factorial, prod
from pathlib import Path
import json, hashlib, time, sys
sys.set_int_max_str_digits(0)

HERE=Path(__file__).resolve().parent
N=128; R=24; L=30; ODD_S=200
TAU=Q(39,25); RHO=Q(5,4); CAP=Q(2048)
DATA_PATH=HERE/'prime_atom_densities_200.json'
DATA=json.loads(DATA_PATH.read_text())
DENS={int(p):{int(a)//int(p):Q(v['numerator'],v['denominator'])
             for a,v in rows.items()}
      for p,rows in DATA['noncore_prime_atom_densities'].items()}
FAC=[factorial(r) for r in range(max(ODD_S,N)+1)]

def show(x):
    with localcontext() as c:
        c.prec=45;c.rounding=ROUND_CEILING
        return str(D(x.numerator)/D(x.denominator))

def record(x):
    return {'numerator':str(x.numerator),'denominator':str(x.denominator),'decimal_upper':show(x)}

def falling(n,s):
    assert 0<=s<=n
    return FAC[n]//FAC[n-s]

def exp_coeff(f,R):
    a=[Q(0)]*(R+1);a[0]=Q(1)
    for r in range(1,R+1):
        a[r]=sum((k*f[k]*a[r-k] for k in range(1,r+1)),Q(0))/r
    return a

def odd_coeff(ell,n,r):
    f=[Q(0)]*(r+1)
    for k,v in DENS[ell].items():
        if k<=r:f[k]=v/Q(n**(k-1))
    if ell==3:f[1]=Q(2)
    return exp_coeff(f,r)[r]

def main():
    start=time.monotonic()
    # Uniform LLL uses a logarithmic criterion, not the weaker Bernoulli one.
    q=TAU/Q(N*(N-1))
    lll_penalty=36*N*q/(1-q)
    t=(TAU-1)/(TAU+1)
    log_tau_lower=2*(t+t**3/3)
    assert lll_penalty<log_tau_lower
    u=(RHO-1)/(RHO+1)
    log_rho_lower=2*(u+u**3/3)
    assert lll_penalty/2<log_rho_lower
    inverse_core_exponent=Q(9,2)*N*N*q/(1-q)
    log_cap_lower=11*Q(56,81)
    assert inverse_core_exponent<log_cap_lower

    # The p=2 polynomial invariant gives a completely uniform infinite tail.
    r0=R+1
    invol_tail=CAP*Q(18**r0*FAC[r0],FAC[2*r0])/(1-Q(9,2*r0+1))
    assert invol_tail<Q(31,10**6)
    # For s>=200, every odd-prime row is below s*(3/5)^s.
    assert Q(729*8,ODD_S**2)<Q(3,5)**3
    a=Q(3,5)
    odd_tail=CAP*a**ODD_S*(Q(ODD_S**2)/(1-a)+2*ODD_S*a/(1-a)**2+a*(1+a)/(1-a)**3)
    assert odd_tail<Q(1,10**36)

    # All odd prime rows of support below 200; (p=3,r=1) forces a core hit.
    odd_sum=Q(0);odd_rows=[]
    for ell in sorted(DENS):
        if ell==2 or ell>=ODD_S:continue
        for r in range(2 if ell==3 else 1,(ODD_S-1)//ell+1):
            s=ell*r;n=max(N,s)
            aa=odd_coeff(ell,n,r)
            value=aa*aa*Q(ell**r*FAC[r]*n**(2*r),(ell-1)*falling(n,s))*min(RHO**s,CAP)
            odd_sum+=value
            odd_rows.append({'p':ell,'r':r,'s':s,'evaluation_n':n,'upper':record(value)})
    assert odd_sum<Q(13,1000)
    (HERE/'connected_128_odd_prime_rows.json').write_text(json.dumps(odd_rows,indent=2)+'\n')

    # Build the positive labelled-atom configuration model in every density box.
    # P[i][r][c] includes the factorial of the r-c noncore ports.
    P=[]
    for i in range(L):
        low=Q(3*i,2*L);high=Q(3*(i+1),2*L)
        f=[Q(0)]*(R+1)
        for k,v in DENS[2].items():
            if k<=R:f[k]=v/Q(N**(k-1))
        f[2]=(Q(21,4)-Q(7,2)*low)/N
        bad=exp_coeff(f,R)
        core=[high**c/FAC[c] for c in range(R+1)]
        P.append([[core[c]*bad[r-c]*FAC[r-c] for c in range(r+1)] for r in range(R+1)])
    W=[Q(0)]+[Q(2**r*N**(2*r),falling(N,2*r))*min(RHO**(2*r),CAP) for r in range(1,R+1)]
    rest=odd_sum+invol_tail+odd_tail
    worst=Q(0);worst_pair=None;empty=None;count=0;worst_coeff=None
    cert_path=HERE/'connected_128_density_boxes.jsonl'
    with cert_path.open('w') as out:
        for i in range(L):
            for j in range(i,L):
                Z=[Q(1)]+[Q(0)]*R
                C=[Q(0)]*(R+1)
                total=Q(0)
                for r in range(1,R+1):
                    Z[r]=sum((P[i][r][c]*P[j][r][d]/FAC[r-c-d]
                              for c in range(r+1) for d in range(r-c+1)
                              if P[i][r][c] and P[j][r][d]),Q(0))
                    C[r]=Z[r]-sum((k*C[k]*Z[r-k] for k in range(1,r)),Q(0))/r
                    assert C[r]>=0,(i,j,r)
                    total+=C[r]*W[r]
                assert total<Q(985,1000),(i,j,show(total))
                assert total+rest<Q(999,1000),(i,j,show(total+rest))
                if i==j==0:
                    empty=total+rest
                    assert empty<Q(1,2)
                if total>worst:
                    worst=total;worst_pair=[i,j];worst_coeff=C
                out.write(json.dumps({'i':i,'j':j,'t_i':[str(Q(3*i,2*L)),str(Q(3*(i+1),2*L))],
                          't_j':[str(Q(3*j,2*L)),str(Q(3*(j+1),2*L))],
                          'connected_involution_upper':record(total),
                          'complete_conditional_upper':record(total+rest),'passed':True},separators=(',',':'))+'\n')
                count+=1
            print(json.dumps({'stage':'boxes','row':i,'checked':count,'worst_upper':show(worst),'seconds':round(time.monotonic()-start,2)}),flush=True)
    assert count==L*(L+1)//2
    summary={'N':N,'involution_connected_cycle_cutoff':R,'density_bins':L,'odd_prime_support_cutoff':ODD_S,
      'arithmetic':'exact Fraction and integer factorials; decimals are upward-rounded display only',
      'tau':str(TAU),'rho':str(RHO),'inverse_core_cap':str(CAP),
      'lll_log_penalty_upper':record(lll_penalty),'log_tau_lower':record(log_tau_lower),
      'log_rho_lower':record(log_rho_lower),'inverse_core_exponent_upper':record(inverse_core_exponent),
      'log_cap_lower':record(log_cap_lower),'involution_infinite_tail':record(invol_tail),
      'odd_prime_finite_sum':record(odd_sum),'odd_prime_infinite_tail':record(odd_tail),
      'odd_prime_rows_checked':len(odd_rows),'density_box_pairs_checked':count,
      'max_connected_involution_upper':record(worst),'worst_pair':worst_pair,
      'max_complete_conditional_upper':record(worst+rest),'empty_core_upper':record(empty),
      'all_conditional_bad_less_than':'999/1000','empty_core_bad_less_than':'1/2',
      'density_file_sha256':hashlib.sha256(DATA_PATH.read_bytes()).hexdigest(),
      'certificate_file_sha256':hashlib.sha256(cert_path.read_bytes()).hexdigest(),
      'verifier_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
      'seconds':time.monotonic()-start,
      'scope':'Every integer n>=128, assuming the accompanying connected-model, affine-count, and prime-polynomial-invariant proofs.'}
    (HERE/'connected_128_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    (HERE/'connected_128_worst_coefficients.json').write_text(json.dumps([record(x) for x in worst_coeff],indent=2)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if k not in ('max_connected_involution_upper','max_complete_conditional_upper','empty_core_upper')},default=str),flush=True)
    print('FINAL_MAX',show(worst+rest),'EMPTY_CORE',show(empty),flush=True)

if __name__=='__main__':main()
