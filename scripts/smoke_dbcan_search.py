"""Exercise real dbCAN searches and the Snakemake pipeline with tiny references.

Linux: install dbcan==5.2.9, snakemake, and diamond, then PYTHONPATH=src python
scripts/smoke_dbcan_search.py. Synthetic sequences test plumbing, not biology.
"""
from pathlib import Path
import random
import subprocess
import sys
import tempfile

import pyhmmer
from kolach.metabolism.dbcan import read_overview
from kolach.metabolism.evaluate import read_tsv


def main():
    with tempfile.TemporaryDirectory() as temporary:
        root = Path(temporary)
        db = root / 'database/dbcan'
        db.mkdir(parents=True)
        rng = random.Random(723)
        sequences = {family: ''.join(rng.choices('ACDEFGHIKLMNPQRSTVWY', k=300)) for family in ('GH13', 'GH15')}
        fasta = root / 'proteins.faa'
        fasta.write_text(''.join(f'>{family}\n{sequence}\n' for family, sequence in sequences.items()) + '>empty\n' + ''.join(rng.choices('ACDEFGHIKLMNPQRSTVWY', k=300)) + '\n')
        reference = root / 'reference.faa'
        reference.write_text(''.join(f'>ref_{family}|{family}\n{sequence}\n' for family, sequence in sequences.items()))
        subprocess.run(['diamond', 'makedb', '--in', str(reference), '--db', str(db / 'CAZy')], check=True)
        alphabet = pyhmmer.easel.Alphabet.amino()
        builder = pyhmmer.plan7.Builder(alphabet)
        background = pyhmmer.plan7.Background(alphabet)
        for filename, subfamily in (('dbCAN.hmm', False), ('dbCAN-sub.hmm', True)):
            with open(db / filename, 'wb') as handle:
                for family, sequence in sequences.items():
                    model_name = f'{family}_e1.hmm|{family}:1|3.2.1.1:1' if subfamily else family + '.hmm'
                    msa = pyhmmer.easel.TextMSA(name=model_name.encode(), sequences=[pyhmmer.easel.TextSequence(name=b'reference', sequence=sequence)])
                    model, _, _ = builder.build_msa(msa.digitize(alphabet), background)
                    model.write(handle)
        (db / 'fam-substrate-mapping.tsv').write_text('substrate\tunused\tfamily\tunused\tec\nstarch\t-\tGH13\t-\t3.2.1.1\nstarch\t-\tGH15\t-\t3.2.1.1\n')
        output = root / 'out'
        command = [sys.executable, '-m', 'kolach.cli', 'annotate', '--protein-fasta', str(fasta),
                   '--database-dir', str(root / 'database'), '--databases', 'dbcan',
                   '--output-dir', str(output), '--threads', '2', '--add-specialisations', '--genome-id', '001']
        subprocess.run(command, check=True)
        calls = read_overview(output / 'dbcan/overview.tsv')
        assert calls['GH13']['dbcan_tools'] == 3, calls
        assert calls['GH15']['dbcan_tools'] == 3, calls
        assert not calls['empty']['dbcan_recommended'], calls
        _, annotations = read_tsv(output / 'kolach_annotations.tsv')
        assert len(annotations) == 3 and all(row['accepted_ko'] == '-' for row in annotations)
        _, product = read_tsv(output / 'kolach_cazy_product.tsv')
        assert product[0]['genome'] == '001' and product[0]['CAZy-Starch'] == 'True', product
        print('PASS: real HMM, subfamily HMM, DIAMOND, integration and metabolic distillation')


if __name__ == '__main__':
    main()
