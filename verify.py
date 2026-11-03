#!/usr/bin/env python3
"""Public verification entry points; stored records are not acceptance oracles."""
from pathlib import Path
from fractions import Fraction
from concurrent.futures import ThreadPoolExecutor
import argparse,hashlib,importlib.util,json,os,shutil,subprocess,sys,tempfile,time
ROOT=Path(__file__).resolve().parent
sys.path[:0]=[str(ROOT/'computation'),str(ROOT/'structural_reduction'),str(ROOT/'global_referee')]

def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()
def save(p,x):Path(p).write_text(json.dumps(x,indent=2),encoding='utf-8')
def rows():return read(ROOT/'CERTIFICATE_INDEX.json')['finite_degrees']
def integrity(out):
    manifest=read(ROOT/'MANIFEST.json')
    for r in manifest['files']:
        rel=Path(r['path']);assert not rel.is_absolute() and '..' not in rel.parts
        p=ROOT/rel;assert p.stat().st_size==r['bytes'] and sha(p)==r['sha256'],str(rel)
    z={'files_checked':len(manifest['files']),'all_pass':True,'mathematical_computation':False}
    save(out/'integrity.json',z);print(z)

def coverage(out):
    from replay_certificates import assignment,sign,valid_permutation
    import certificate_codec as codec
    reports=[]
    for r in rows():
        n=r['degree'];text=(ROOT/r['groups']).read_text()
        gens=assignment(text,'ExplicitGenerators');k=len(gens)
        assert assignment(text,'ExplicitDegree')==n and k==r['covers']
        odd=[]
        for gg in gens:
            assert all(valid_permutation(g,n) for g in gg)
            odd.append(any(sign(g)==-1 for g in gg))
        e=sum(not x for x in odd);expected=k*(k+1)//2+e*(e+1)//2
        assert e==r['even_covers'] and expected==r['witnesses']
        cert=ROOT/r['certificate']
        if cert.suffix=='.kc4':
            with cert.open('rb') as f:h=codec.read_header(f)
            assert h['has_odd_generator']==odd and h['degree']==n
            z=codec.stream_decode(cert);assert z['records']==expected
        else:
            pool=[]
            if r.get('permutation_pool'):
                for line in (ROOT/r['permutation_pool']).read_text().splitlines():
                    p=json.loads(line);assert p['index']==len(pool)+1
                    assert valid_permutation(p['x'],n) and sign(p['x'])==1
                    pool.append(p['x'])
            flags=bytearray(k*(k+1)//2);count=0
            qf=(ROOT/r['queries']).open() if r.get('queries') else None
            if qf:assert int(next(qf))==expected
            with cert.open() as f:
                for line in f:
                    a=json.loads(line);i,j,p=a['i'],a['j'],a['parity']
                    assert a['n']==n and 1<=i<=j<=k and p in (-1,1)
                    assert p==1 or not(odd[i-1] or odd[j-1])
                    offset=k*(i-1)-(i-1)*(i-2)//2+j-i;bit=1 if p==1 else 2
                    assert not flags[offset]&bit
                    flags[offset]|=bit;count+=1
                    if 'x' in a:assert valid_permutation(a['x'],n) and sign(a['x'])==p
                    else:
                        c=a['candidate'];assert 1<=c<=len(pool)
                        if qf:assert list(map(int,next(qf).split()))==[i,j,p,c]
            if qf:assert qf.read().strip()=='';qf.close()
            assert count==expected
            for i in range(k):
                for j in range(i,k):
                    offset=k*i-i*(i-1)//2+j-i
                    assert flags[offset]==(3 if not(odd[i] or odd[j]) else 1)
        reports.append({'degree':n,'covers':k,'even_covers':e,'witnesses':expected,'all_slots_present':True})
        if n%10==0:print('coverage through',n,flush=True)
    assert sum(z['witnesses'] for z in reports)==38972586
    save(out/'coverage.json',{'all_pass':True,'witnesses':38972586,'degrees':107,'records':reports})

def models(out):
    import independent_structural_and_density_audit as a
    import audit_baseline_group_models as b
    import audit_residual_group_models as c
    z={'finite_40_63':a.audit_cover_partition(),'finite_64_127':a.audit_extended_residuals(),'densities':a.audit_density_table()}
    b.main(out/'baseline_models.json');c.main(out/'residual_models.json')
    z['actual_generators']=read(out/'baseline_models.json')['actual_generators_recognized']+read(out/'residual_models.json')['actual_generators_directly_recognized']
    assert z['actual_generators']==703569
    claims={}
    for rel in ['lifting/INDUCTION_CLOSED_THROUGH_63.json','hyperplane_lifting/INDUCTION_CLOSED_THROUGH_127.json']:
        summary=read(ROOT/'computation'/rel)
        for row in summary['closed_degree_records']:
            n=row['degree']
            if n in claims:assert claims[n]==row
            claims[n]=row
    assert set(claims)==set(range(21,128))
    edges=0
    for n,row in claims.items():
        path=ROOT/'computation'/row['file'];assert sha(path)==row['sha256']
        if n<40:continue
        claim=read(path);part=ROOT/'computation'/claim['inheritance_partition_file']
        assert sha(part)==claim['inheritance_partition_sha256']
        partition=read(part)
        deps=claim['dependencies']
        assert [d['degree'] for d in deps]==partition['required_dependencies' if n<64 else 'conservative_required_dependencies']
        for d in deps:
            assert 21<=d['degree']<n and d==claims[d['degree']]
            edges+=1
    assert edges==5086
    z['finite_claims']=107;z['strictly_decreasing_dependencies']=edges
    z['all_pass']=True;save(out/'models.json',z)

def math(out):
    work=out/'math_work';work.mkdir(exist_ok=True)
    for name in ['prime_order_referee','analysis','final_referee']:
        dest=work/name;dest.mkdir(exist_ok=True)
        for p in (ROOT/name).iterdir():
            if p.is_file() and (p.suffix=='.py' or p.name=='prime_atom_densities_200.json'):
                shutil.copy2(p,dest/p.name)
    def run(rel):
        with (out/(Path(rel).stem+'.log')).open('wb') as log:
            p=subprocess.run([sys.executable,str(work/rel)],stdout=log,stderr=subprocess.STDOUT)
        assert p.returncode==0,rel
    run('prime_order_referee/prime_atom_densities.py')
    assert read(work/'prime_order_referee/prime_atom_densities_200.json')==read(ROOT/'prime_order_referee/prime_atom_densities_200.json')
    run('prime_order_referee/verify_involution_invariant.py')
    run('prime_order_referee/connected_128_certificate.py')
    run('analysis/connected_128_interval_verifier.py')
    run('final_referee/connected_model_small_audit.py')
    s=read(work/'prime_order_referee/connected_128_summary.json')
    def q(x):return Fraction(int(x['numerator']),int(x['denominator']))
    assert s['density_box_pairs_checked']==465 and s['odd_prime_rows_checked']==267
    assert q(s['max_complete_conditional_upper'])<Fraction(999,1000)
    assert q(s['empty_core_upper'])<Fraction(1,2)
    for key in ['max_complete_conditional_upper','empty_core_upper']:
        assert q(s[key])==q(read(ROOT/'prime_order_referee/connected_128_summary.json')[key])
    save(out/'math_summary.json',{'all_pass':True,'density_data_regenerated':True,'primitive_coefficients':7251,'density_box_pairs':465,'odd_prime_rows':267,'exact_summary':s,'interval_summary':read(work/'analysis/connected_128_interval_summary.json')})

def graph(out,degrees,jobs):
    import certificate_codec as codec
    engine=ROOT/'software_referee/runtime_293/graph_engine'
    def command(args):
        if os.name!='nt':return [str(engine),*map(str,args)]
        def wp(p):
            p=str(Path(p).resolve()).replace('\\','/')
            return '/mnt/'+p[0].lower()+p[2:]
        return ['wsl.exe','--exec',wp(engine),'check',*[wp(p) for p in args[1:5]],str(args[5])]
    selected=[r for r in rows() if r['method']=='graph' and (not degrees or r['degree'] in degrees)]
    assert selected
    def one(r):
        n=r['degree'];q=ROOT/r['queries'] if r['queries'] else out/f'queries_{n}.txt'
        if r['queries'] is None:
            assert not q.exists(),'Use a fresh output directory for a new replay.'
            codec.stream_decode(ROOT/r['certificate'],q,'queries')
        log=out/f'graph_{n}.log'
        with (out/f'graph_{n}_driver.log').open('wb') as f:
            p=subprocess.run(command(['check',ROOT/r['graph'],ROOT/r['pool'],q,log,'traces']),stdout=f,stderr=f)
        assert p.returncode==0,n
        text=log.read_text();c=r['witnesses']
        assert f'queries={c} pass={c} fail=0 ' in text and 'nauty_version=29300' in text
        z={'degree':n,'queries':c,'passed':True,'engine_sha256':sha(engine),'graph_sha256':sha(ROOT/r['graph']),'pool_sha256':sha(ROOT/r['pool']),'query_sha256':sha(q),'log_sha256':sha(log)}
        save(out/f'graph_{n}.json',z);print('graph pass',n,c,flush=True);return z
    with ThreadPoolExecutor(max_workers=jobs) as pool:z=list(pool.map(one,selected))
    save(out/'graph_summary.json',{'all_pass':True,'degrees':len(z),'queries':sum(a['queries'] for a in z),'records':z})

def gap(out,binary,library):
    assert binary and library,'Supply --gap-binary and --gap-root.'
    work=out/'gap_work';work.mkdir(exist_ok=True)
    shutil.copy2(ROOT/'computation/replay_certificates.py',work/'replay_certificates.py')
    for r in rows()[:6]:
        for key in ['groups','certificate']:
            p=ROOT/r[key];shutil.copy2(p,work/p.name)
    env=os.environ.copy();env.update(GAP_BINARY=str(Path(binary).resolve()),GAP_ROOT=str(Path(library).resolve()))
    for prefix,ns in [('affine',[21,22]),('affine_pruned',[23,24,25,26])]:
        with (out/(prefix+'_driver.log')).open('wb') as log:
            p=subprocess.run([sys.executable,str(work/'replay_certificates.py'),prefix,*map(str,ns)],env=env,stdout=log,stderr=log)
        assert p.returncode==0,prefix
    z=[read(work/(('affine' if n<=22 else 'affine_pruned')+f'_{n}_replay_audit.json')) for n in range(21,27)]
    assert all(a['deterministic_gap_replay_pass'] for a in z)
    assert sum(a['recorded_witness_count'] for a in z)==301061
    save(out/'gap_summary.json',{'all_pass':True,'witnesses':301061,'records':z})

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode',choices=['integrity','coverage','models','math','graph','gap'])
    p.add_argument('--output',type=Path,default=ROOT/'verification_results')
    p.add_argument('--degrees',type=int,nargs='+');p.add_argument('--jobs',type=int,default=4)
    p.add_argument('--gap-binary');p.add_argument('--gap-root');a=p.parse_args()
    a.output.mkdir(parents=True,exist_ok=True)
    if a.mode=='graph':graph(a.output,a.degrees,a.jobs)
    elif a.mode=='gap':gap(a.output,a.gap_binary,a.gap_root)
    else:globals()[a.mode](a.output)
if __name__=='__main__':main()

