"""Compare current legacy-mode outputs against a separately built baseline executable."""
import argparse
import json
import subprocess
from pathlib import Path
from run import sha

def main():
    ap=argparse.ArgumentParser();ap.add_argument('results',type=Path);ap.add_argument('baseline_exe',type=Path)
    args=ap.parse_args();root=args.results.resolve();exe=args.baseline_exe.resolve()
    manifest=json.loads((root/'manifest.json').read_text(encoding='utf-8'))
    dest=root/'legacy_regression';dest.mkdir(exist_ok=False)
    records=[]
    for run in manifest['runs']:
        if run['mode'] not in ['cgal-adaptive','cgal-adaptive-radius']: continue
        cell=dest/run['model']/run['mode'];cell.mkdir(parents=True)
        command=run['command'][:];command[0]=str(exe);command[2]=str(cell/'output.ply')
        index=command.index('--radius-diagnostics');del command[index:index+2]
        (cell/'command.json').write_text(json.dumps(command,indent=2),encoding='utf-8')
        with (cell/'stdout.log').open('w') as stdout,(cell/'stderr.log').open('w') as stderr:
            proc=subprocess.run(command,stdout=stdout,stderr=stderr)
        record={'model':run['model'],'mode':run['mode'],'returncode':proc.returncode,
                'baseline_commit':manifest['code']['baseline_commit'],'baseline_exe_sha256':sha(exe),
                'baseline_sha256':sha(cell/'output.ply') if proc.returncode==0 else None,
                'current_sha256':run.get('output_sha256')}
        record['identical']=record['baseline_sha256']==record['current_sha256'] and proc.returncode==0
        records.append(record)
        (dest/'comparison.json').write_text(json.dumps(records,indent=2),encoding='utf-8')
        print(record,flush=True)
    if not all(r['identical'] for r in records): raise SystemExit('Legacy regression mismatch')

if __name__=='__main__':main()
