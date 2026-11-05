"""Copy this completed study into an existing publication checkout.

No deletion, no git commit/push. Raw local logs remain untouched. Public logs
omit machine/license header lines; numeric solver and field records remain.
"""
from pathlib import Path
import argparse,hashlib,json,re,shutil

ROOT=Path(__file__).resolve().parents[1]
def sha(b):return hashlib.sha256(b).hexdigest()
def main(destination):
    dest=Path(destination).resolve()
    if dest==ROOT or ROOT in dest.parents:raise ValueError('Publication destination must be outside source study')
    dest.mkdir(parents=True,exist_ok=True);records=[]
    for source in sorted(ROOT.rglob('*')):
        if not source.is_file():continue
        rel=source.relative_to(ROOT)
        if any(x in ('__pycache__','.git') for x in rel.parts):continue
        if source.suffix in ('.pyc','.bak') or source.name.startswith('.report_pdf'):continue
        if rel.parts[0]=='qa' and (source.suffix=='.png' or 'partial' in source.name):continue
        if source.name=='PUBLIC_COPY_MANIFEST.json':continue
        if rel.parts[:2]==('fdtd','cases'):
            case=ROOT/'fdtd/cases'/rel.parts[2]
            if not (case/'DONE.json').exists():
                if not (case/'CANCELLED.json').exists():continue
                if source.name not in ('CANCELLED.json','settings.json','engine_stdout.log','engine_stderr.log','build_this_case.lsf','run_this_case.lsf'):continue
        payload=source.read_bytes();redacted=False
        if source.suffix=='.log':
            text=payload.decode('utf8',errors='replace')
            clean=[]
            for line in text.splitlines(keepends=True):
                if re.search(r'license|flexnet|host name|hostname|user name|computer name|machine name|license_server',line,re.I):
                    clean.append('[Machine/license header omitted from public copy]\n');redacted=True
                else:clean.append(line)
            payload=''.join(clean).encode('utf8')
        target=dest/rel;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(payload)
        records.append(dict(path=rel.as_posix(),bytes=len(payload),sha256=sha(payload),
            source_sha256=sha(source.read_bytes()),log_header_redacted=redacted))
    (dest/'PUBLIC_COPY_MANIFEST.json').write_text(json.dumps(dict(files=records,
        policy='Only machine/license log header lines and transient build files are excluded; scientific arrays retained unchanged',
        file_count=len(records),total_bytes=sum(x['bytes'] for x in records)),indent=2),encoding='utf8')
    print(json.dumps(dict(destination=str(dest),files=len(records),total_MiB=sum(x['bytes'] for x in records)/1024**2,
        redacted_logs=sum(x['log_header_redacted'] for x in records))))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--destination',required=True);a=p.parse_args();main(a.destination)
