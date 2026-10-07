"""Thin dbCAN 5.2.9 runner and explicit CAZy evidence conversion.

The upstream tool owns search thresholds, domain overlap and recommendations.
Kolach imports those results and never makes CAZy families vote for a KO.
"""
from collections import defaultdict
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
from urllib.request import urlopen

from .evaluate import DATA_DIR, read_tsv
from .specialise import parse_bool, read_annotations, write_tsv

VERSION = '5.2.9'
RELEASE = 'db_v5-2-9_5-5-2026'
BASE_URL = f'https://dbcan.s3.us-west-2.amazonaws.com/{RELEASE}'
DATABASE_FILES = ('CAZy.dmnd', 'dbCAN.hmm', 'dbCAN-sub.hmm', 'fam-substrate-mapping.tsv')
DEFINITIONS = DATA_DIR / 'cazy_definitions.tsv'
LOW_MATCHING = {'GH63', 'GH49', 'GH127', 'CE2', 'GH23', 'PL30', 'PL8', 'PL33', 'GH88', 'GH67', 'GH98', 'CE4'}
DBCAN_COLUMNS = ['dbcan_id', 'dbcan_sub_id', 'dbcan_sub_ec', 'dbcan_sub_substrate',
                 'dbcan_family_ec',
                 'dbcan_tools', 'dbcan_recommended', 'dbcan_hmm_candidates',
                 'dbcan_sub_candidates', 'dbcan_diamond_candidates']


def file_hash(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def download_dbcan(database_dir):
    """Download four versioned upstream files, publish only complete downloads.

    Existing files are reused only if their recorded release and hash match.
    Hashes identify downloaded bytes; upstream does not publish checksums here.
    """
    directory = Path(database_dir)
    directory.mkdir(parents=True, exist_ok=True)
    manifest_path = directory / 'kolach_dbcan_database.json'
    old = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
    manifest = {'dbcan_version': VERSION, 'release': RELEASE, 'files': {}}
    for name in DATABASE_FILES:
        target = directory / name
        if target.exists():
            recorded = old.get('files', {}).get(name, {})
            if old.get('release') != RELEASE or file_hash(target) != recorded.get('sha256'):
                raise ValueError(f'Unverified existing database file: {target}; choose a new database directory')
        else:
            with tempfile.NamedTemporaryFile(dir=directory, delete=False) as handle:
                temporary = Path(handle.name)
                try:
                    with urlopen(f'{BASE_URL}/{name}', timeout=120) as response:
                        shutil.copyfileobj(response, handle)
                    if handle.tell() == 0:
                        raise ValueError(f'Empty database download: {name}')
                except BaseException:
                    handle.close()
                    temporary.unlink(missing_ok=True)
                    raise
            temporary.replace(target)
        manifest['files'][name] = {'url': f'{BASE_URL}/{name}', 'sha256': file_hash(target), 'bytes': target.stat().st_size}
        # Save after each file so interrupted downloads can resume safely.
        manifest_path.write_text(json.dumps(manifest, indent=2) + '\n', encoding='utf-8')
    return manifest_path


def run_dbcan(protein_fasta, database_dir, output_dir, threads=1):
    """Run the released protein annotation command with all three methods."""
    if not isinstance(threads, int) or threads < 1:
        raise ValueError('threads must be a positive integer')
    executable = shutil.which('run_dbcan')
    if executable is None:
        raise ValueError('Install dbcan==5.2.9 in the annotation environment (see docs/dbcan.md)')
    version = subprocess.check_output([executable, 'version'], text=True).strip()
    if not re.fullmatch(r'(?:dbCAN version:\s*)?5\.2\.9', version):
        raise ValueError(f'Expected dbCAN {VERSION}; executable reports {version!r}')
    db = Path(database_dir).resolve()
    for name in DATABASE_FILES:
        if not (db / name).is_file():
            raise ValueError(f'Missing dbCAN reference: {db / name}')
    out = Path(output_dir).resolve()
    out.mkdir(parents=True, exist_ok=True)
    # Every search runs in a fresh directory. Previous raw results cannot leak
    # into an empty/failed rerun. Publish outputs only after successful execution.
    with tempfile.TemporaryDirectory(dir=out) as temporary:
        command = [executable, 'CAZyme_annotation', '--input_raw_data', str(Path(protein_fasta).resolve()),
                   '--mode', 'protein', '--db_dir', str(db), '--output_dir', temporary,
                   '--methods', 'diamond,hmm,dbCANsub', '--threads', str(threads)]
        subprocess.run(command, check=True)
        if not (Path(temporary) / 'overview.tsv').is_file():
            raise ValueError('dbCAN succeeded without producing overview.tsv')
        for source in Path(temporary).iterdir():
            if source.is_file():
                source.replace(out / source.name)
    overview = out / 'overview.tsv'
    if not overview.is_file():
        raise ValueError('dbCAN succeeded without producing overview.tsv')
    provenance = {'dbcan_version': version, 'command': command,
                  'protein_fasta': {'path': str(Path(protein_fasta).resolve()), 'sha256': file_hash(protein_fasta)},
                  'references': {name: file_hash(db / name) for name in DATABASE_FILES},
                  'overview_sha256': file_hash(overview)}
    (out / 'kolach_dbcan_run.json').write_text(json.dumps(provenance, indent=2) + '\n', encoding='utf-8')
    return overview


def families(text):
    """Extract whole dbCAN family/subfamily names, dropping domain coordinates."""
    if text in ('', '-', 'N/A'):
        return set()
    result = set()
    for token in re.split(r'[+|;]', text):
        token = re.sub(r'\([^)]*\)', '', token).strip()
        if not re.fullmatch(r'(?:GH|GT|CBM|AA|CE|PL)\d+(?:_[A-Za-z0-9]+)*', token):
            raise ValueError(f'Invalid dbCAN family: {token!r}')
        result.add(token)
    return result


def read_overview(path):
    columns, rows = read_tsv(path)
    required = {'Gene ID', 'EC#', 'dbCAN_hmm', 'dbCAN_sub', 'DIAMOND', '#ofTools', 'Recommend Results'}
    if not required <= set(columns):
        raise ValueError(f'dbCAN 5 overview requires columns: {sorted(required)}')
    annotations = {}
    for row in rows:
        gene = row['Gene ID']
        if not gene or gene in annotations:
            raise ValueError(f'Empty or duplicate dbCAN gene ID: {gene!r}')
        hmm = families(row['dbCAN_hmm'])
        sub = families(row['dbCAN_sub'])
        diamond = families(row['DIAMOND'])
        tools = int(row['#ofTools'])
        if tools != sum(bool(hits) for hits in (hmm, sub, diamond)):
            raise ValueError(f'dbCAN tool count disagrees with evidence: {gene}')
        recommended = families(row['Recommend Results'])
        if tools < 2 and recommended:
            raise ValueError(f'dbCAN recommendation has fewer than two tools: {gene}')
        ecs = set()
        pairs = set()
        sub_domains = row['dbCAN_sub'].split('+') if sub else []
        ec_domains = row['EC#'].split('|') if sub else []
        # dbCAN writes both lists from the same ordered sub_results rows. Keep
        # that association: no Cartesian product of every family with every EC.
        if len(sub_domains) != len(ec_domains):
            raise ValueError(f'dbCAN subfamily/EC domain counts disagree: {gene}')
        for domain, text in zip(sub_domains, ec_domains):
            names = families(domain)
            for entry in text.split(';'):
                ec = entry.split(':')[0]
                if ec in ('', '-'):
                    continue
                if not re.fullmatch(r'\d+\.(?:\d+|-)\.(?:\d+|-)\.(?:\d+|-)', ec):
                    raise ValueError(f'Invalid dbCAN EC: {ec!r}')
                if names <= recommended:
                    ecs.add(ec)
                    if len(names) == 1:
                        pairs.add(next(iter(names)).split('_')[0] + ';EC:' + ec)
        confident = tools >= 2 and bool(recommended)
        annotations[gene] = {
            'gene_id': gene, 'dbcan_id': ';'.join(sorted({name.split('_')[0] for name in recommended})) if confident else '-',
            'dbcan_sub_id': ';'.join(sorted(name for name in recommended if '_' in name)) if confident else '-',
            'dbcan_sub_ec': ';'.join(sorted(ecs)) if confident else '-',
            'dbcan_sub_substrate': (row.get('Substrate') or '-') if sub <= recommended else '-',
            'dbcan_family_ec': '|'.join(sorted(pairs)) if confident else '-',
            'dbcan_tools': tools, 'dbcan_recommended': confident,
            'dbcan_hmm_candidates': row['dbCAN_hmm'], 'dbcan_sub_candidates': row['dbCAN_sub'],
            'dbcan_diamond_candidates': row['DIAMOND'],
        }
        if not confident:
            annotations[gene]['dbcan_sub_substrate'] = '-'
    return annotations


def merge_overview(annotation_table, overview, output_file):
    """Attach orthogonal CAZy columns by exact gene ID, retaining the inventory."""
    columns, rows = read_tsv(annotation_table)
    if 'gene_id' not in columns or len({row['gene_id'] for row in rows}) != len(rows):
        raise ValueError('dbCAN join requires globally unique gene_id values')
    if set(columns) & set(DBCAN_COLUMNS):
        raise ValueError('Annotation table already contains dbCAN columns')
    calls = read_overview(overview)
    unknown = set(calls) - {row['gene_id'] for row in rows}
    if unknown:
        raise ValueError(f'dbCAN genes absent from annotation inventory: {sorted(unknown)[:5]}')
    empty = dict.fromkeys(DBCAN_COLUMNS, '-')
    empty.update(dbcan_tools=0, dbcan_recommended=False)
    merged = [{**row, **{column: calls.get(row['gene_id'], empty)[column] for column in DBCAN_COLUMNS}} for row in rows]
    metadata = {'overview': {'path': str(Path(overview).resolve()), 'sha256': file_hash(overview)},
                'annotations': {'path': str(Path(annotation_table).resolve()), 'sha256': file_hash(annotation_table)},
                'import_policy': 'upstream recommendations; at least two tools; globally unique exact gene IDs',
                'recommended_genes': sum(bool(call['dbcan_recommended']) for call in calls.values())}
    write_tsv(output_file, columns + DBCAN_COLUMNS, merged)
    (Path(output_file).parent / 'kolach_dbcan_import.json').write_text(json.dumps(metadata, indent=2) + '\n', encoding='utf-8')
    return Path(output_file)


def cazy_calls(features, definitions_file=DEFINITIONS):
    """DRAM1 all-cleavage calls plus AnnoGuild's unique-family signature gate."""
    _, definitions = read_tsv(definitions_file)
    groups = defaultdict(list)
    substrates = defaultdict(set)
    for row in definitions:
        required = {feature.strip() for feature in row['function_ids'].split(',')}
        groups[row['function_name']].append(required)
        for feature in required:
            substrates[feature].add(row['function_name'])
    signatures = {feature for feature, names in substrates.items() if len(names) == 1} - LOW_MATCHING
    calls = []
    details = []
    for genome, observed in features.items():
        row = {'genome': genome}
        for substrate, steps in groups.items():
            hits = set().union(*steps) & observed
            signature_hits = hits & signatures
            # DRAM1 requires every cleavage row. Its Boolean product is 0/1,
            # and AnnoGuild's >=0.6 gate therefore has the same interpretation.
            present = all(step & observed for step in steps) and bool(signature_hits)
            name = 'CAZy-' + re.sub(r'[-\s()]', '_', substrate).strip()
            row[name] = present
            details.append({'genome': genome, 'pathway': name, 'called': present,
                            'cleavage_steps': len(steps), 'steps_present': sum(bool(step & observed) for step in steps),
                            'signature_features': ';'.join(sorted(signature_hits)), 'supporting_features': ';'.join(sorted(hits))})
        calls.append(row)
    return calls, details


def make_cazy_product(annotation_table, output_dir, genome_id=None, genome_column=None, gene_genome_map=None):
    """Distill confirmed family/KO evidence into the existing metabolism input."""
    inventory = read_annotations(annotation_table, genome_id=genome_id, genome_column=genome_column, gene_genome_map=gene_genome_map)
    columns, rows = read_tsv(annotation_table)
    if len({row['gene_id'] for row in rows}) != len(rows):
        raise ValueError('CAZy distillation requires globally unique gene IDs')
    if not {'dbcan_id', 'dbcan_recommended'} <= set(columns):
        raise ValueError('CAZy distillation requires a merged dbCAN annotation table')
    calls = {row['gene_id']: row for row in rows}
    features = {}
    for genome, genes in inventory.items():
        observed = set().union(*genes.values()) if genes else set()
        for gene in genes:
            row = calls[gene]
            if parse_bool(row['dbcan_recommended']):
                observed.update(families(row['dbcan_id']))
        features[genome] = observed
    product, details = cazy_calls(features)
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    write_tsv(out / 'kolach_cazy_product.tsv', list(product[0]) if product else ['genome'], product)
    write_tsv(out / 'kolach_cazy_pathways.tsv', ['genome', 'pathway', 'called', 'cleavage_steps', 'steps_present', 'signature_features', 'supporting_features'], details)
    metadata = {'target_dbcan_version': VERSION, 'download_reference_release': RELEASE,
                'definition_source': 'dylancronin/AnnoGuild@291b232907e6568597178e21384ff0f964993a82/resources/DRAM_CAZy_definitions.tsv',
                'definition_sha256': file_hash(DEFINITIONS), 'annotation_sha256': file_hash(annotation_table),
                'genome_id': genome_id, 'genome_column': genome_column,
                'genome_map_sha256': file_hash(gene_genome_map) if gene_genome_map else None,
                'policy': 'dbCAN recommendations with at least two tools; DRAM1 all-cleavage calls; AnnoGuild signature gate'}
    (out / 'kolach_cazy_metadata.json').write_text(json.dumps(metadata, indent=2) + '\n', encoding='utf-8')
    return out / 'kolach_cazy_product.tsv'
