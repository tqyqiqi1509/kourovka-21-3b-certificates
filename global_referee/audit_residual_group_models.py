#!/usr/bin/env python3
"""Direct affine/wreath recognition for all actual residual group generators.

Targets the new finite fields and deeper words in degrees 40 through 127.
Does not invoke the graph encoder, GAP, nauty, Traces, or any intersection test.
"""
from pathlib import Path
from math import prod
import json,time
import audit_baseline_group_models as model
from independent_structural_and_density_audit import ROOT,HERE,PP,sha,word_order,words_through

def main(output=None):
    started=time.monotonic();rows=[];all_factors=set();new_factors=set();coordinate_bindings={}
    for dirname,bound,lo,hi,template_name in (
        ('lifting',63,40,63,'lifted_residual_covers_through_63.json'),
        ('hyperplane_lifting',127,64,127,'extended_residual_covers_through_127.json')):
        directory=ROOT/'computation'/dirname
        coordinates=directory/f'affine_coordinates_through_{bound}.json'
        model.COORDINATES=json.loads(coordinates.read_text())
        model.affine_member.cache_clear();model.word_member.cache_clear()
        coordinate_bindings[str(coordinates.relative_to(ROOT))]=sha(coordinates)
        for q_text,coords in model.COORDINATES.items():
            q=int(q_text);p,d=PP[q]
            assert len(coords)==q and len(set(map(tuple,coords)))==q
            assert all(len(v)==d and all(0<=x<p for x in v) for v in coords)
        template=directory/template_name
        data=json.loads(template.read_text());words=[tuple(w) for w in data['words']]
        assert words==words_through(bound)
        for n in range(lo,hi+1):
            path=directory/f'graph_{n}_groups.g'
            text=path.read_text();assert model.assignment(text,'ExplicitDegree')==n
            covers=model.assignment(text,'ExplicitCoverWords')
            assert covers==[[i+1 for i in c] for c in data['covers_0based'][str(n)]]
            orders=model.assignment(text,'ExplicitGroupOrders')
            generators=model.assignment(text,'ExplicitGenerators')
            assert len(covers)==len(orders)==len(generators)
            checks=0;odd_flags=[];max_depth=0
            for c,order,gens in zip(covers,orders,generators):
                ww=[words[i-1] for i in c]
                assert sum(map(prod,ww))==n and order==prod(map(word_order,ww))
                hasodd=False
                for w in ww:
                    all_factors.update(w)
                    new_factors.update(q for q in w if q>40)
                    max_depth=max(max_depth,len(w))
                for perm in gens:
                    assert len(perm)==n and sorted(perm)==list(range(1,n+1))
                    p=tuple(x-1 for x in perm);hasodd|=model.sign(p)==-1
                    offset=0
                    for w in ww:
                        d=prod(w);section=tuple(p[offset+j]-offset for j in range(d))
                        assert model.word_member(w,section),(n,c,w,perm)
                        offset+=d
                    checks+=1
                assert hasodd==bool(any(w and (PP[w[0]][0]%2 or w[0]<=4) for w in ww))
                odd_flags.append(hasodd)
            ga=directory/f'graph_{n}_generator_audit.json'
            graph_audit_bound=False
            if ga.exists():
                G=json.loads(ga.read_text())
                assert G['has_odd_generator']==odd_flags and G['groups_checked']==len(covers)
                assert G['explicit_groups_sha256']==sha(path)
                graph_audit_bound=True
            rows.append({'n':n,'residual_covers':len(covers),'actual_generators_directly_recognized':checks,
                         'max_word_depth':max_depth,'all_even_covers':sum(not o for o in odd_flags),
                         'all_expected_group_orders_match':True,'all_theoretical_parities_match':True,
                         'graph_generator_audit_already_present_and_bound':graph_audit_bound,
                         'group_file':str(path.relative_to(ROOT)),'group_file_sha256':sha(path),
                         'template_sha256':sha(template),'coordinates_sha256':sha(coordinates)})
        print('RESIDUAL_MODELS_PASS',lo,hi,sum(r['actual_generators_directly_recognized'] for r in rows if lo<=r['n']<=hi),flush=True)
    result={'all_pass':True,'degrees':list(range(40,128)),'records':rows,
            'actual_generators_directly_recognized':sum(r['actual_generators_directly_recognized'] for r in rows),
            'residual_cover_count':sum(r['residual_covers'] for r in rows),
            'primitive_degrees_checked':sorted(all_factors),'primitive_degrees_above_40_checked':sorted(new_factors),
            'largest_word_depth':max(r['max_word_depth'] for r in rows),
            'coordinate_files':coordinate_bindings,'membership_engine_sha256':sha(HERE/'audit_baseline_group_models.py'),
            'independent_word_module_sha256':sha(HERE/'independent_structural_and_density_audit.py'),
            'source_sha256':sha(__file__),'seconds':time.monotonic()-started,
            'scope':'Actual generator membership, theoretical order values and parity. This is not a graph-intersection replay or a finite-theorem closure claim.'}
    (Path(output) if output is not None else HERE/'residual_group_models_audit.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='records'},indent=2),flush=True)

if __name__=='__main__':
    main()
