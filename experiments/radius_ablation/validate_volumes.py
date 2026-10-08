"""Gate relative volume on CGAL self-intersection and enclosed-volume checks."""
import argparse
import json
import subprocess
from pathlib import Path
import pandas as pd
from run import sha

def main():
    ap=argparse.ArgumentParser();ap.add_argument('results',type=Path);ap.add_argument('validator',type=Path)
    args=ap.parse_args();root=args.results.resolve();exe=args.validator.resolve()
    manifest=json.loads((root/'manifest.json').read_text(encoding='utf-8'))
    fields=json.loads((root/'field_summary.json').read_text(encoding='utf-8'))
    table=pd.read_csv(root/'summary.csv');records=[]
    def validate(path):
        command=[str(exe),str(path)];proc=subprocess.run(command,capture_output=True,text=True)
        record={'path':str(path),'command':command,'returncode':proc.returncode,
                'stdout':proc.stdout,'stderr':proc.stderr}
        if proc.returncode==0: record.update(json.loads(proc.stdout))
        else: record.update(bounds_volume=False,validation_error='CGAL validation failed')
        records.append(record)
        (root/'volume_validation.json').write_text(json.dumps({'executable_sha256':sha(exe),'records':records},indent=2),encoding='utf-8')
        print(path,record.get('self_intersects'),record['bounds_volume'],flush=True)
        return record
    for source in manifest['data']['inputs']:
        name=Path(source['path']).stem;reference=validate(source['path'])
        fields[name]['input']['cgal_volume_validation']=reference
        for run in manifest['runs']:
            if run['model']!=name or run['status']!='success':continue
            cell=Path(run['path']);checked=validate(cell/'output.ply')
            row=json.loads((cell/'metrics.json').read_text(encoding='utf-8'))
            valid=reference['bounds_volume'] and checked['bounds_volume'] and fields[name]['input']['volume']>0 and row['vol']>0
            row.update(volume_valid=bool(valid), self_intersects=checked.get('self_intersects'),
                       input_self_intersects=reference.get('self_intersects'),bounds_volume=checked['bounds_volume'],
                       relative_volume_change=row['vol']/fields[name]['input']['volume']-1 if valid else None)
            (cell/'metrics.json').write_text(json.dumps(row,indent=2),encoding='utf-8')
            mask=(table.model==name)&(table['mode']==run['mode'])
            for col in ['volume_valid','self_intersects','input_self_intersects','bounds_volume','relative_volume_change']:
                table.loc[mask,col]=row[col]
    table.to_csv(root/'summary.csv',index=False)
    (root/'field_summary.json').write_text(json.dumps(fields,indent=2),encoding='utf-8')
    if any(r['returncode'] for r in records):raise SystemExit('Volume validation failure; see volume_validation.json')

if __name__=='__main__':main()
