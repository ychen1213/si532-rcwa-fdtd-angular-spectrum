"""Verify publication bytes against the package manifest AND staged Git blobs."""
from pathlib import Path
import argparse,hashlib,json,subprocess

def main(repository,subdirectory):
    repo=Path(repository).resolve();folder=(repo/subdirectory).resolve()
    if repo not in folder.parents:raise ValueError('Study must be inside repository')
    manifest=json.loads((folder/'PUBLIC_COPY_MANIFEST.json').read_text(encoding='utf8'))
    raw=subprocess.check_output(['git','ls-files','--stage','-z','--',subdirectory],cwd=repo)
    index={}
    for row in raw.split(b'\0'):
        if not row:continue
        header,path=row.split(b'\t',1);mode,sha,stage=header.split()
        if stage!=b'0':raise RuntimeError('Unmerged index entry')
        index[path.decode('utf8')]=sha.decode('ascii')
    for item in manifest['files']:
        path=folder/item['path'];payload=path.read_bytes()
        if hashlib.sha256(payload).hexdigest()!=item['sha256']:raise RuntimeError('Manifest changed: '+item['path'])
        name=path.relative_to(repo).as_posix()
        digest=hashlib.sha1(b'blob '+str(len(payload)).encode()+b'\0'+payload).hexdigest()
        if index.get(name)!=digest:raise RuntimeError('Ignored, filtered or changed Git blob: '+name)
    print(json.dumps(dict(manifest_files=len(manifest['files']),indexed_study_files=len(index),
        manifest_SHA256_verified=True,exact_Git_blob_bytes_verified=True,
        line_ending_conversion=False,total_bytes=manifest['total_bytes'])))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--repository',required=True);p.add_argument('--subdirectory',required=True);a=p.parse_args();main(a.repository,a.subdirectory)
