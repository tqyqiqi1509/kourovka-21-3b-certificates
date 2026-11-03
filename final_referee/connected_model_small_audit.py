"""Independent port-matching enumeration versus the positive-model formal log.

This checks the semantic combinatorial interface, not the 128-box calculation.
No source code from the connected certificate is imported.
"""
from fractions import Fraction as Q
from itertools import permutations
from math import factorial, prod
from pathlib import Path
import json

R=5
# First weight is a core degree-one vertex. The others are noncore types.
LEFT=[(1,True,Q(3,2)),(2,False,Q(5,4)),(3,False,Q(2,3))]
RIGHT=[(1,True,Q(7,5)),(2,False,Q(11,6)),(3,False,Q(1,2))]

def profiles(types,r):
    def recurse(k,remaining,counts):
        if k==len(types):
            if remaining==0: yield counts
            return
        degree=types[k][0]
        for count in range(remaining//degree+1):
            yield from recurse(k+1,remaining-count*degree,counts+[count])
    for counts in recurse(0,r,[]):
        ports=[]; weight=Q(1); vertex=0
        for (degree,core,q),count in zip(types,counts):
            weight*=q**count/factorial(count)
            for _ in range(count):
                ports.extend([(vertex,core)]*degree);vertex+=1
        yield ports,vertex,weight

all_coeff=[Q(1)]; connected=[Q(0)]
matchings_examined=0
for r in range(1,R+1):
    total=Q(0); conn=Q(0)
    for lp,lv,lw in profiles(LEFT,r):
        for rp,rv,rw in profiles(RIGHT,r):
            all_count=0; conn_count=0
            for perm in permutations(range(r)):
                matchings_examined+=1
                if any(lp[i][1] and rp[perm[i]][1] for i in range(r)): continue
                adj=[set() for _ in range(lv+rv)]
                for i,j in enumerate(perm):
                    a=lp[i][0];b=lv+rp[j][0]
                    adj[a].add(b);adj[b].add(a)
                seen={0};stack=[0]
                while stack:
                    a=stack.pop()
                    for b in adj[a]-seen:
                        seen.add(b);stack.append(b)
                all_count+=1
                if len(seen)==lv+rv:
                    assert lv+rv<=r+1
                    conn_count+=1
            total+=lw*rw*all_count
            conn+=lw*rw*conn_count
    all_coeff.append(total);connected.append(conn)

def coefficient_exp_bad(types,r):
    return sum((weight for ports,v,weight in profiles([t for t in types if not t[1]],r)),Q(0))
H=[coefficient_exp_bad(LEFT,r) for r in range(R+1)]
K=[coefficient_exp_bad(RIGHT,r) for r in range(R+1)]
a=LEFT[0][2];b=RIGHT[0][2]
Z=[Q(1)]
for r in range(1,R+1):
    Z.append(sum((H[r-c]*K[r-d]*a**c/factorial(c)*b**d/factorial(d)
                  *factorial(r-c)*factorial(r-d)/factorial(r-c-d)
                  for c in range(r+1) for d in range(r-c+1)),Q(0)))
assert Z==all_coeff
formal=[Q(0)]
for r in range(1,R+1):
    formal.append(Z[r]-sum((k*formal[k]*Z[r-k] for k in range(1,r)),Q(0))/r)
assert formal==connected
assert all(v>=0 for v in connected)
report={'status':'passed','maximum_edges':R,'port_matchings_examined':matchings_examined,
        'method':'Direct enumeration of port bijections; BFS tests connectedness; no formal log used in the direct count.',
        'all_coefficients':[str(v) for v in all_coeff],
        'connected_coefficients':[str(v) for v in connected],
        'formal_log_agrees_exactly':True}
path=Path(__file__).with_suffix('.json')
path.write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
