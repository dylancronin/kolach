"""DRAM2 expressions, evaluated directly without graphs or Python eval.

Each immutable node is (kind, value, children). Boolean operators require
explicit grouping. Commas are steps, not OR (unlike the AnnoGuild grammar).
Upstream tables are imported unchanged; this module owns only their execution.
"""
from collections import Counter
import csv
from dataclasses import dataclass
import json
from pathlib import Path
import operator
import re

from .evaluate import read_tsv
from .specialise import read_annotations, write_tsv

DATA = Path(__file__).parent / 'data/dram2'
OPS = {name: getattr(operator, name) for name in ('gt', 'ge', 'lt', 'le', 'eq', 'ne')}


@dataclass(frozen=True)
class Node:
    kind: str
    value: str = ''
    children: tuple = ()


def parse_rule(text):
    tokens = re.findall(r'->|`(?:[^`\\]|\\.)*`|"(?:[^"\\]|\\.)*"|[A-Za-z0-9_.:;\-]+|\S', text)
    position = 0

    def take(expected=None):
        nonlocal position
        if position == len(tokens) or expected is not None and tokens[position] != expected:
            raise ValueError(f'Expected {expected or "token"} in {text!r}')
        token = tokens[position]
        position += 1
        return token

    def peek():
        return tokens[position] if position < len(tokens) else ''

    def atom():
        token = take()
        if token == '[':
            node = steps()
            take(']')
            return node
        if token == '@':
            return Node('alias', name(take()))
        if token.startswith('"'):
            return Node('literal', json.loads(token))
        if not re.fullmatch(r'[A-Za-z0-9_.:;\-]+|`(?:[^`\\]|\\.)*`', token):
            raise ValueError(f'Unexpected {token!r} in {text!r}')
        if peek() == '(':
            take('(')
            args = []
            if peek() != ')':
                args.append(expression())
                while peek() == ',':
                    take(',')
                    args.append(expression())
            take(')')
            return Node('call', token, tuple(args))
        if peek() == '@':
            take('@')
            # Upstream currently ignores database qualifiers during evaluation.
            token = take()
        return Node('name', name(token))

    def pipe():
        parts = [atom()]
        while peek() == '->':
            take('->')
            parts.append(atom())
        return parts[0] if len(parts) == 1 else Node('pipe', children=tuple(parts))

    def expression():
        parts = [pipe()]
        separator = peek()
        if separator not in ('&', '|'):
            return parts[0]
        while peek() == separator:
            take(separator)
            parts.append(pipe())
        if peek() in ('&', '|'):
            raise ValueError(f'Mixed AND/OR requires square brackets: {text!r}')
        return Node('and' if separator == '&' else 'or', children=tuple(parts))

    def steps():
        parts = [expression()]
        while peek() == ',':
            take(',')
            parts.append(expression())
        return parts[0] if len(parts) == 1 else Node('steps', children=tuple(parts))

    def name(token):
        return token[1:-1].replace(r'\`', '`').replace('\\\\', '\\') if token.startswith('`') else token

    result = steps()
    if position != len(tokens):
        raise ValueError(f'Trailing tokens in {text!r}')
    return result


def load_rules(path, common=DATA / 'common.tsv'):
    """Expand aliases once; reject missing aliases and cycles.

    Last definition wins, matching DRAM2's dictionaries. Some ecosystem sheets
    repeat aliases or output names; preserve every metadata row in the report.
    """
    def read_rules(source):
        with open(source, encoding='utf-8-sig', newline='') as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            if not reader.fieldnames or not {'name', 'rule'} <= set(reader.fieldnames):
                raise ValueError(f'Rule TSV requires name and rule columns: {source}')
            # Upstream loads curated tables with truncate_ragged_lines=True.
            # A few ecosystem tables contain extra trailing metadata cells.
            return [{key: value or '' for key, value in row.items() if key is not None}
                    for row in reader]
    rows = read_rules(path)
    shared = read_rules(common)
    definitions = {}
    outputs = []
    output_nodes = {}
    for row in rows + shared:
        if not row.get('rule', '').strip():
            continue
        node = parse_rule(row['rule'])
        alias = row.get('alias', row.get('name', '').replace(' ', ''))
        if alias:
            definitions[alias] = node
        if row.get('name'):
            output_nodes[row['name']] = node
            outputs.append(row)
    cache = {}

    def expand(node, stack=()):
        if node.kind == 'alias':
            if node.value in stack:
                raise ValueError(f'Alias cycle: {stack + (node.value,)}')
            if node.value not in definitions:
                raise ValueError(f'Undefined alias: {node.value}')
            if node.value not in cache:
                cache[node.value] = expand(definitions[node.value], stack + (node.value,))
            return cache[node.value]
        return Node(node.kind, node.value, tuple(expand(child, stack) for child in node.children))

    return [(row, expand(output_nodes[row['name']])) for row in outputs]


def evaluate(node, counts, rows=()):
    """Return a Boolean/count, row mask, or pathway coverage dictionary.

    Upstream casts every at_least step to bool, including COUNTS mode. Preserve
    that behavior explicitly; gene multiplicity does not increase its score.
    """
    kind, fn, args = node.kind, node.value, node.children
    ev = lambda child: evaluate(child, counts, rows)
    if kind == 'name':
        return counts.get(fn.upper(), 0)
    if kind in ('and', 'steps'):
        return all(ev(child) for child in args)
    if kind == 'or':
        return any(ev(child) for child in args)
    if kind == 'pipe':
        mask = [True] * len(rows)
        for child in args[:-1]:
            part = evaluate(child, counts, rows)
            if not isinstance(part, list):
                raise ValueError('Only row filters can precede ->')
            mask = [a and b for a, b in zip(mask, part)]
        return evaluate(args[-1], counts, [row for row, keep in zip(rows, mask) if keep])
    if kind != 'call':
        raise ValueError(f'Expected an evaluable expression: {node}')
    value = lambda i: args[i].value
    arity = {'not': 1, 'tax': 1, 'percent': 2, 'at_least': 3,
             'column_count_values': 5, 'column_sum_values': 3,
             'filter_contains': 2, 'filter_compare': 3}
    if fn not in arity and fn not in ('path_steps', 'path_subunits'):
        raise ValueError(f'Unsupported DRAM2 function: {fn}')
    if fn in arity and len(args) != arity[fn]:
        raise ValueError(f'{fn} requires {arity[fn]} arguments')
    if fn == 'not':
        result = ev(args[0])
        return [not item for item in result] if isinstance(result, list) else not result
    if fn == 'tax':
        return any(re.search(value(0), row.get('taxonomy') or '') for row in rows)
    if fn in ('percent', 'at_least'):
        group = args[1 if fn == 'percent' else 2]
        if group.kind not in ('steps', 'or'):
            raise ValueError(f'{fn} requires steps or OR alternatives')
        hits = sum(bool(ev(child)) for child in group.children)
        threshold = int(value(0))
        if fn == 'percent':
            return hits / len(group.children) >= threshold / 100
        if value(1).upper() not in ('PRESENCE', 'COUNTS'):
            raise ValueError('at_least mode must be PRESENCE or COUNTS')
        return hits >= threshold
    if fn in ('filter_contains', 'filter_compare'):
        if not any(value(0) in row for row in rows):
            if fn == 'filter_contains':
                return [False] * len(rows)
            raise ValueError(f'Missing filter column: {value(0)}')
        if fn == 'filter_contains':
            return [bool(re.search(value(1), row.get(value(0)) or '')) for row in rows]
        return [row.get(value(0)) not in (None, '') and OPS[value(1)](float(row[value(0)]), float(value(2))) for row in rows]
    if fn in ('column_count_values', 'column_sum_values'):
        if not any(value(0) in row for row in rows):
            return False
        numbers = [float(row[value(0)]) for row in rows if row.get(value(0)) not in (None, '')]
        if fn == 'column_sum_values':
            return OPS[value(1)](sum(numbers), float(value(2)))
        hits = sum(OPS[value(1)](number, float(value(2))) for number in numbers)
        return OPS[value(3)](hits, float(value(4)))
    if fn in ('path_steps', 'path_subunits'):
        results = []
        for child in args:
            parts = child.children if fn == 'path_subunits' and child.kind == 'and' else (child,)
            results.extend(bool(ev(part)) for part in parts)
        unit = 'steps' if fn == 'path_steps' else 'subunits'
        return {unit: len(results), unit + '_present': sum(results),
                'coverage_percentage': sum(results) / len(results) if results else 0.0}
    raise ValueError(f'Unhandled DRAM2 function: {fn}')


def annotation_features(row, ko_column='accepted_ko'):
    """Use accepted KOs and explicit supplementary IDs, never KO candidates.

    Supplementary columns follow DRAM2's separators. Family+EC joins are made
    on one gene (the upstream adapter explodes these columns on one gene too).
    """
    features = set()
    separators = {'dbcan_id': ';', 'dbcan_sub_id': ';', 'dbcan_sub_substrate': ';',
                  'fegenie_id': None, 'sulfur_id': None, 'camper_id': None,
                  'methyl_id': ',', 'merops_family': ';', 'peptidase_family': ';'}
    if row.get(ko_column, '') not in ('', '-', 'None', 'nan'):
        features.update(re.split(r'[,;\s]+', row[ko_column]))
    for column, separator in separators.items():
        raw = row.get(column) or ''
        if raw in ('', '-'):
            continue
        parts = raw.split(separator) if separator else [raw]
        features.update(part.strip().split(' ')[0] for part in parts)
    for column in ('pfam_hits', 'pfam_id'):
        features.update(re.findall(r'\[(PF\d{5})\.\d+\]', row.get(column) or ''))
    for column in ('camper_EC', 'kegg_EC', 'kofam_EC'):
        features.update(re.findall(r'\[EC:(\d+\.\d+\.\d+\.\d+)\]', row.get(column) or ''))
    features.update('EC:' + ec for ec in re.findall(r'\(EC\s+([\d.\-]+)\)', row.get('cazy_hits') or ''))
    features.update('EC:' + ec for ec in re.findall(r'([\d.\-]+)', row.get('cazy_subfam_ec') or ''))
    ecs = {part.split(':')[0] for part in (row.get('dbcan_sub_ec') or '').split(';') if part and part != '-'}
    features.update('EC:' + ec for ec in ecs)
    # Family+EC substrate rules require correlated evidence. The dbCAN importer
    # supplies selected-domain pairs; unrelated domains must not be cross-joined.
    for pair in (row.get('dbcan_family_ec') or '').split('|'):
        if pair in ('', '-'):
            continue
        if not re.fullmatch(r'(?:GH|GT|CBM|AA|CE|PL)\d+;EC:\d+\.(?:\d+|-)\.(?:\d+|-)\.(?:\d+|-)', pair):
            raise ValueError(f'Invalid correlated CAZy family/EC feature: {pair!r}')
        features.add(pair)
    return {feature.upper() for feature in features if feature and feature != '-'}


def run_dram2(annotation_table, output_dir, genome_id=None, genome_column=None,
              gene_genome_map=None, ko_column='accepted_ko', rules_file=None,
              common_file=None, system='default'):
    inventory = read_annotations(annotation_table, ko_column, genome_id, genome_column, gene_genome_map)
    columns, annotations = read_tsv(annotation_table)
    grouped = {genome: [] for genome in inventory}
    gene_genomes = {}
    for genome, genes in inventory.items():
        for gene in genes:
            gene_genomes.setdefault(gene, []).append(genome)
    for row in annotations:
        matches = gene_genomes[row['gene_id']] if genome_id is not None or gene_genome_map is not None else [row[genome_column or 'genome']]
        if len(matches) != 1:
            raise ValueError(f'Ambiguous gene/genome assignment: {row["gene_id"]}')
        grouped[matches[0]].append(row)
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    systems = {'default': 'rules.tsv', 'ag': 'ag_rules.tsv', 'bgc': 'rules_bgc.tsv',
               'eng_sys': 'rules_eng_sys.tsv', 'gut': 'rules_gut.tsv', 'marine': 'rules_marine.tsv'}
    tables = {'traits': DATA / 'traits.tsv', 'product': rules_file or DATA / systems[system]}
    results = []
    for category, path in tables.items():
        rules = load_rules(path, common_file or DATA / 'common.tsv')
        for genome, rows in grouped.items():
            counts = Counter(feature for row in rows for feature in annotation_features(row, ko_column))
            for metadata, node in rules:
                options = json.loads(metadata.get('rule_options') or '{}').get('include_when', {})
                if not set(options.get('all_columns_present', [])) <= set(columns) or not set(options.get('all_columns_absent', [])).isdisjoint(columns):
                    continue
                try:
                    result = evaluate(node, counts, rows)
                except (ValueError, KeyError) as error:
                    raise ValueError(f'{category}/{metadata["name"]}/{genome}: {error}') from error
                coverage = result if isinstance(result, dict) else {'present': bool(result)}
                results.append({'genome': genome, 'category': category, 'name': metadata['name'],
                                'group': metadata.get('group', metadata.get('topic_ecosystem', '')),
                                'present': coverage.get('present', ''),
                                'coverage_percentage': coverage.get('coverage_percentage', ''),
                                'units': coverage.get('steps', coverage.get('subunits', '')),
                                'units_present': coverage.get('steps_present', coverage.get('subunits_present', ''))})
    write_tsv(out / 'kolach_dram2.tsv', ['genome', 'category', 'name', 'group', 'present', 'coverage_percentage', 'units', 'units_present'], results)
    import hashlib
    inputs = {'annotations': annotation_table, 'common': common_file or DATA / 'common.tsv', **tables}
    if gene_genome_map:
        inputs['genome_map'] = gene_genome_map
    provenance = {'sources': json.loads((DATA / 'manifest.json').read_text()), 'system': system,
                  'ko_column': ko_column, 'genome_id': genome_id, 'genome_column': genome_column,
                  'policy': 'accepted KOs plus explicitly supplied supplementary evidence; absent IDs have zero counts',
                  'inputs': {name: {'path': str(Path(path).resolve()), 'sha256': hashlib.sha256(Path(path).read_bytes()).hexdigest()} for name, path in inputs.items()}}
    (out / 'kolach_dram2_metadata.json').write_text(json.dumps(provenance, indent=2) + '\n', encoding='utf-8')
    return out / 'kolach_dram2.tsv'
