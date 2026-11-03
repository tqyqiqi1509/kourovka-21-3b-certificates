# Exact certificate search for the complete affine-envelope cover.
# Randomness is used only to propose a permutation. Each accepted certificate
# is checked by the exact permutation-group intersection algorithm.
# Parameters may be pre-bound by a driver: CoverInput, OutputPrefix,
# DegreesToRun, MaxProposals, MaxIntersections, Seed, OnlyDiagonal.
if not IsBound(CoverInput) then CoverInput := "../structural_reduction/affine_covers_through_24.g"; fi;
if not IsBound(OutputPrefix) then OutputPrefix := "affine"; fi;
if not IsBound(DegreesToRun) then DegreesToRun := [21..24]; fi;
if not IsBound(MaxProposals) then MaxProposals := 200000; fi;
if not IsBound(MaxIntersections) then MaxIntersections := 3000; fi;
if not IsBound(Seed) then Seed := 21320261102; fi;
if not IsBound(OnlyDiagonal) then OnlyDiagonal := false; fi;
SizeScreen([1000,1000]);
Read(CoverInput);
AffinePrimitive := function(q)
 local f,p,d,vectors,lin,gens,m,e,i,g,ord;
 if q<=4 then return SymmetricGroup(q);fi;
 f:=FactorsInt(q);p:=f[1];d:=Length(f);
 if not ForAll(f,x->x=p) then Error("factor is not a prime power");fi;
 vectors:=Tuples(AsList(GF(p)),d);
 lin:=GL(d,p);
 gens:=List(GeneratorsOfGroup(lin),m->PermList(List(vectors,v->Position(vectors,v*m))));
 for i in [1..d] do
  e:=List([1..d],j->Zero(GF(p)));e[i]:=One(GF(p));
  Add(gens,PermList(List(vectors,v->Position(vectors,v+e))));
 od;
 g:=Group(gens);
 ord:=q*Product([0..d-1],i->q-p^i);
 if Size(g)<>ord or not IsTransitive(g,[1..q]) then Error("AGL construction failed");fi;
 return g;
end;
PrimitiveCache:=[];
for q in Set(Concatenation(AffineWords)) do PrimitiveCache[q]:=AffinePrimitive(q);od;
WordGroup := function(w)
 local g,i,deg,expected;
 if Length(w)=0 then return Group(());fi;
 g:=PrimitiveCache[w[1]];deg:=w[1];
 for i in [2..Length(w)] do
  expected:=Size(g)^w[i]*Size(PrimitiveCache[w[i]]);
  g:=WreathProduct(g,PrimitiveCache[w[i]]);deg:=deg*w[i];
  if Size(g)<>expected or LargestMovedPoint(g)<>deg or not IsTransitive(g,[1..deg]) then Error("wreath construction failed");fi;
 od;
 return g;
end;
WordGroups:=List(AffineWords,WordGroup);
CoverGroup := function(c,n)
 local offset,gens,k,d,g,perm,images,j,res,expected;
 offset:=0;gens:=[];expected:=1;
 for k in c do
  d:=AffineWordDegrees[k];g:=WordGroups[k];expected:=expected*Size(g);
  for perm in GeneratorsOfGroup(g) do
   images:=[1..n];
   for j in [1..d] do images[offset+j]:=offset+j^perm;od;
   Add(gens,PermList(images));
  od;
  offset:=offset+d;
 od;
 if offset<>n then Error("cover degree mismatch");fi;
 if Length(gens)=0 then gens:=[()];fi;
 res:=Group(gens);
 if Size(res)<>expected then Error("direct product order mismatch");fi;
 return res;
end;
CoreData := function(g,n)
 local trans,a,b,blocks,labels,k;
 trans:=[];
 for a in [1..n-1] do for b in [a+1..n] do
  if (a,b) in g then Add(trans,(a,b));fi;
 od;od;
 if Length(trans)=0 then trans:=[()];fi;
 blocks:=Orbits(Group(trans),[1..n]);labels:=List([1..n],i->0);
 for k in [1..Length(blocks)] do for a in blocks[k] do labels[a]:=k;od;od;
 return rec(labels:=labels,blocks:=Filtered(blocks,b->Length(b)>1));
end;
CoreDisjoint := function(a,b,x)
 local block,seen,j,l;
 for block in b.blocks do
  seen:=[];
  for j in block do
   l:=a.labels[j^x];
   if IsBound(seen[l]) then return false;fi;
   seen[l]:=true;
  od;
 od;
 return true;
end;
