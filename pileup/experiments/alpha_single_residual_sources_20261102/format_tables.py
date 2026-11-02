"""Normalize CSV line endings with lossless row comparison and retained originals."""
from pathlib import Path
import csv
import hashlib
import io
import json
import shutil

exp=Path(__file__).resolve().parent
root=exp.parents[2]
rec=root/'jupyter/run_records'/exp.name/'report_format_before'
revisions=exp/'report-revisions'
rec.mkdir(exist_ok=False)
revisions.mkdir(exist_ok=False)
for name in ('report-source-freeze.json','report-output-freeze.json','delivery-freeze.json','validation.json'):
    shutil.copyfile(exp/name,revisions/('v1-'+name))
rows=[]
for p in sorted((exp/'tables').glob('*.csv')):
    old=p.read_bytes()
    new=old.replace(b'\r\n',b'\n')
    assert old!=new and b'\r' not in new
    assert list(csv.reader(io.StringIO(old.decode())))==list(csv.reader(io.StringIO(new.decode())))
    (rec/p.name).write_bytes(old)
    p.write_bytes(new)
    rows.append(dict(path=str(p.relative_to(root)),old_sha256=hashlib.sha256(old).hexdigest(),
        new_sha256=hashlib.sha256(new).hexdigest(),all_parsed_rows_equal=True,
        original=str((rec/p.name).relative_to(root))))
source=root/'pileup/analysis/alpha_single_residual_sources_report.py'
new_source=exp/'source/pileup/analysis/alpha_single_residual_sources_report_v2.py'
shutil.copyfile(source,new_source)
wrapper=exp/'run_report.py'
shutil.copyfile(wrapper,revisions/'v1-run_report.py')
wrapper.write_text(wrapper.read_text().replace('from pileup.analysis.alpha_single_residual_sources_report import main',
    'from pileup.analysis.alpha_single_residual_sources_report_v2 import main'))
freeze=json.loads((exp/'report-source-freeze.json').read_text())
for p in (source,new_source,wrapper,Path(__file__).resolve()):
    freeze[str(p.relative_to(root))]=hashlib.sha256(p.read_bytes()).hexdigest()
(exp/'report-source-freeze.json').write_text(json.dumps(freeze,indent=2)+'\n')
freeze=json.loads((exp/'report-output-freeze.json').read_text())
for row in rows:freeze[row['path']]=row['new_sha256']
(exp/'report-output-freeze.json').write_text(json.dumps(freeze,indent=2)+'\n')
(exp/'report-format-revision.json').write_text(json.dumps(dict(
    reason='Git staged whitespace check flagged CSV CRLF; export now explicitly uses LF',
    numerical_change=False,source_snapshot_v1_retained=True,
    executed_display_cells='original frozen report v1; figures and notebook outputs unchanged',
    actual_table_rewrite='format_tables.py; parsed rows exactly equal',rows=rows),ensure_ascii=False,indent=2)+'\n')
print('Three CSV files normalized to LF; all parsed rows identical; v1 sources and hashes retained')
