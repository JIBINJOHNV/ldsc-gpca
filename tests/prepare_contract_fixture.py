"""Capture deterministic preparation products for before/after source comparison.

Run with PYTHONPATH pointing to each checkout:
  python tests/prepare_contract_fixture.py /tmp/preparation-parity before
  python tests/prepare_contract_fixture.py /tmp/preparation-parity after
Then compare before.json and after.json (only output paths are normalized).
Worker-attempt timings are excluded; all tables, QC, manifests and settings stay.
"""
import contextlib
import io
import json
import sys
from pathlib import Path
from ldsc_gpca import prepare
root=Path(sys.argv[1]).resolve();root.mkdir(exist_ok=True, parents=True)
tag=sys.argv[2]
source=root/'source.vcf'
header=['##fileformat=VCFv4.2']+[f'##contig=<ID={i}>' for i in range(1,23)]
header += [f'##FORMAT=<ID={name},Number=1,Type=Float,Description="fixture">' for name in ['AF','ES','SE','LP','NEF','SI']]
header += ['#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample']
records=[f'{i}\t1\trs{i}\tA\tG\t.\tPASS\t.\tAF:ES:SE:LP:NEF:SI\t.2:.1:.05:4:1000:.95' for i in range(1,23)]
records += [f'1\t{i}\t{ident}\tA\tG\t.\tPASS\t.\tAF:ES:SE:LP:NEF:SI\t{vals}' for i,ident,vals in [
(2,'duplicate','.2:.1:.05:400:2000:.9'),(3,'duplicate','.2:.1:.05:1:2000:.9'),
(4,'missing_info','.2:.1:.05:4:1000:.'),(5,'bad_se','.2:.1:0:4:1000:.9')]]
source.write_text('\n'.join(header+records)+'\n')
hm=root/'hm3';hm.write_text('SNP\n'+'\n'.join(f'rs{i}' for i in range(1,23))+'\nduplicate\n')
coords=root/'coords';coords.write_text('SNP\n'+'\n'.join(f'{i}_1_A_G' for i in range(1,23))+'\n1_2_A_G\n1_3_A_G\n')
manifest=root/'manifest.csv';manifest.write_text('traitname,vcf_files,N\n002,source.vcf,\n001,source.vcf,\n')
result={}
for label,kwargs in [('default',{}),('split',{'splitby_chr':'split'}),
('legacy',{'write_munge_inputs':True,'hapmap_file':coords,'munge_id_source':'chr_pos_ref_alt'}),
('ldsc',{'mode':'ldsc','hapmap_file':hm}),('both',{'mode':'both','hapmap_file':hm})]:
    out=root/tag/label
    with contextlib.redirect_stdout(io.StringIO()):
        prepare.prepare_inputs(manifest,out,**{'splitby_chr':'nosplit','prepare_workers':1,**kwargs})
    for path in sorted(out.rglob('*')):
        if path.is_file() and 'Worker_Attempts' not in path.name:
            result[f'{label}/{path.relative_to(out)}']=path.read_text().replace(str(out),'<OUT>')
manifest.write_text('traitname,vcf_files,N\n002,source.vcf,45000\n001,source.vcf,45000\n')
out=root/tag/'both_override'
with contextlib.redirect_stdout(io.StringIO()):
    prepare.prepare_inputs(manifest,out,splitby_chr='nosplit',prepare_workers=1,mode='both',hapmap_file=hm)
for path in sorted(out.rglob('*')):
    if path.is_file() and 'Worker_Attempts' not in path.name:
        result[f'both_override/{path.relative_to(out)}']=path.read_text().replace(str(out),'<OUT>')
(root/(tag+'.json')).write_text(json.dumps(result,sort_keys=True,indent=2))
print(tag,len(result),'tables and audit files captured')
