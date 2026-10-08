"""Run the fixed 3 x 5 matrix, retaining all commands, failures, and artifacts."""
import argparse
import datetime as dt
import hashlib
import json
import platform
from pathlib import Path
import subprocess
import time

ROOT = Path(__file__).resolve().parents[2]
MODELS = {'wine_glass24': '0.0214378565', 'camel2': '0.001', 'shark13': '0.003106642746'}
MODES = ['cgal-adaptive', 'cgal-adaptive-pmap', 'cgal-adaptive-radius-cap',
         'cgal-adaptive-radius-recovery', 'cgal-adaptive-radius']

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def git(*args):
    return subprocess.check_output(['git', *args], cwd=ROOT, text=True).strip()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--exe', type=Path, default=ROOT/'build/Release/rar_cgal.exe')
    ap.add_argument('--results', type=Path, default=Path(__file__).parent/'results'/dt.datetime.now().strftime('%Y%m%d_%H%M%S'))
    args = ap.parse_args()
    args.exe = args.exe.resolve(); out = args.results.resolve()
    out.mkdir(parents=True, exist_ok=False)  # never overwrite an earlier run
    status = git('status','--porcelain')
    tracked = git('status','--porcelain','--untracked-files=no')
    if tracked:
        raise RuntimeError('Commit tracked source changes before comparison runs: '+tracked)
    manifest = {
        'experiment': {'id': out.name, 'title': 'Radius mechanism ablation', 'status': 'running',
                       'hypothesis': 'Recovery reduces wine over-refinement; Cap explains shark improvement; camel changes little'},
        'code': {'branch': git('branch','--show-current'), 'commit': git('rev-parse','HEAD'),
                 'dirty': bool(status), 'tracked_dirty': bool(tracked), 'status': status,
                 'baseline_commit': '5f460e724eb8e94e84f5308880499089523c8eaf'},
        'build': {'executable': str(args.exe), 'sha256': sha(args.exe),
                  'generator': 'Visual Studio 17 2022', 'build_type': 'Release',
                  'CGAL': '6.1.1-I-900', 'triplet': 'x64-windows'},
        'data': {'inputs': []},
        'design': {'min_edge': .001, 'max_edge': .5, 'iterations': 5, 'relax_steps': 3,
                   'preserve_features': False, 'projection': True, 'Chen': False,
                   'samples_per_direction': 80000, 'reference_seed': 731, 'output_seed': 343},
        'execution': {'started_at': dt.datetime.now().isoformat(), 'host': platform.platform(),
                      'artifact_root': str(out)}, 'runs': []}
    for name, eps in MODELS.items():
        source = ROOT/'data'/f'{name}.ply'
        manifest['data']['inputs'].append({'path': str(source), 'sha256': sha(source) if source.exists() else None,
                                          'epsilon': eps, 'provenance': 'existing repository data; no substitution'})
        for mode in MODES:
            cell = out/name/mode; cell.mkdir(parents=True)
            cmd = [str(args.exe), str(source), str(cell/'output.ply'), '--field', mode, '--epsilon', eps,
                   '--min-edge', '0.001', '--max-edge', '0.5', '--iterations', '5', '--relax-steps', '3',
                   '--no-preserve-features', '--no-export-field', '--radius-diagnostics', str(cell/'initial.csv')]
            (cell/'command.json').write_text(json.dumps(cmd,indent=2),encoding='utf-8')
            (cell/'command.txt').write_text(subprocess.list2cmdline(cmd),encoding='utf-8')
            record = {'model': name, 'mode': mode, 'command': cmd, 'path': str(cell)}
            print('START', name, mode, flush=True)
            start=time.perf_counter()
            try:
                with (cell/'stdout.log').open('w',encoding='utf-8') as stdout, (cell/'stderr.log').open('w',encoding='utf-8') as stderr:
                    proc=subprocess.run(cmd,stdout=stdout,stderr=stderr,cwd=ROOT,check=False)
                record.update(returncode=proc.returncode,status='success' if proc.returncode==0 and (cell/'output.ply').exists() else 'failed')
                if record['status']=='failed': record['failure']=(cell/'stderr.log').read_text(encoding='utf-8') or 'Nonzero exit or missing output'
            except Exception as exc:
                record.update(status='failed', failure=repr(exc))
            record['wall_seconds']=time.perf_counter()-start
            if (cell/'output.ply').exists(): record['output_sha256']=sha(cell/'output.ply')
            (cell/'status.json').write_text(json.dumps(record,indent=2),encoding='utf-8')
            manifest['runs'].append(record)
            (out/'manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
            print('DONE',name,mode,record['status'],round(record['wall_seconds'],2),flush=True)
    manifest['experiment']['status']='complete'
    manifest['execution']['finished_at']=dt.datetime.now().isoformat()
    manifest['results']={'status':'success' if all(r['status']=='success' for r in manifest['runs']) else 'partial',
                         'failures':[r for r in manifest['runs'] if r['status']!='success']}
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    print('RESULTS',out,flush=True)

if __name__=='__main__': main()
