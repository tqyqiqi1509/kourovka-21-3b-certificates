from pathlib import Path
import sys,json,time
R=Path(__file__).resolve().parent
sys.path[:0]=[str(R/'computation'),str(R/'structural_reduction'),str(R/'global_referee')]
from replay_certificates import assignment,sign
from independent_structural_and_density_audit import words_through
import affine_graphs as old
import affine_hyperplane_graphs as new

def main(out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True);reports=[]
    for row in json.loads((R/'CERTIFICATE_INDEX.json').read_text())['finite_degrees']:
        if row['method']!='graph':continue
        n=row['degree'];bound=40 if n<40 else 63 if n<64 else 127
        directory=R/'computation'
        if n>=40:directory/= 'lifting' if n<64 else 'hyperplane_lifting'
        coord=directory/('affine_coordinates.json' if n<40 else f'affine_coordinates_through_{bound}.json')
        coordinates=json.loads(coord.read_text());words=words_through(bound)
        data=(R/row['groups']).read_text();covers=assignment(data,'ExplicitCoverWords');gens=assignment(data,'ExplicitGenerators')
        module=old if n<64 else new
        count=0
        with (R/row['graph']).open() as f:
            first=f.readline().strip();assert first=='KOUROVKA_GRAPH_COVERS_V1',first
            assert list(map(int,f.readline().split()))==[n,len(covers)]
            for index,(cover,gg) in enumerate(zip(covers,gens),1):
                g=module.cover_graph(words,[i-1 for i in cover],coordinates)
                header=list(map(int,f.readline().split()))
                assert header==[index,len(g.colors),len(g.edges),int(any(sign(p)==-1 for p in gg))],(n,header)
                assert list(map(int,f.readline().split()))==g.colors
                edges={tuple(map(int,f.readline().split())) for _ in range(len(g.edges))}
                assert edges==g.edges
                # Every actual generator extends to this exact colored graph.
                for p in gg:
                    ok,why=module.verify_permutation(g,p);assert ok,(n,why)
                    count+=1
            assert not f.read().strip()
        reports.append({'degree':n,'graphs':len(covers),'actual_generators':count,'all_pass':True})
        print('graph-model pass',n,count,flush=True)
    result={'all_pass':True,'degrees':len(reports),'actual_generators':sum(r['actual_generators'] for r in reports),'records':reports}
    (out/'graph_models.json').write_text(json.dumps(result,indent=2))
if __name__=='__main__':main(sys.argv[1] if len(sys.argv)>1 else R/'verification_results')

