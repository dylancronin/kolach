"""GraftM marker refinement: placement lookup, gating, parsing, and bundles.

No GraftM install is needed: a stand-in `graftM` executable writes placements.
"""
import csv
import io
import json
import os
from pathlib import Path
import tarfile
import tempfile
import unittest
from unittest.mock import patch

from kolach.methods import graftm

FIXTURES = Path(__file__).resolve().parent / "fixtures"

FAKE_GRAFTM = """#!/usr/bin/env python3
# Stand-in for `graftM graft`: place each query by the placement named in
# $GRAFTM_PLACEMENTS ({package: {sequence_sha: placement}}); others get no hit.
import hashlib, json, os, pathlib, sys
args = sys.argv[1:]
if pathlib.Path(sys.argv[0]).name != "graftM" or not args or args[0] != "graft":
    raise SystemExit(0)
query = pathlib.Path(args[args.index("--forward") + 1])
package = pathlib.Path(args[args.index("--graftm_package") + 1]).name[:-len(".gpkg")]
out = pathlib.Path(args[args.index("--output_directory") + 1]) / query.stem
out.mkdir(parents=True, exist_ok=True)
wanted = json.loads(os.environ["GRAFTM_PLACEMENTS"]).get(package, {})
records, name = {}, None
for line in query.read_text().splitlines():
    if line.startswith(">"):
        name = line[1:]
        records[name] = ""
    else:
        records[name] += line.strip()
with open(out / f"{query.stem}_read_tax.tsv", "w") as handle:
    for name, sequence in records.items():
        placement = wanted.get(sequence)
        if placement:
            handle.write(f"{name}\\t{placement}\\n")
"""


def fake_packages(directory):
    """Package folders whose taxonomies contain exactly the placement table's clades."""
    directory = Path(directory)
    lineages = {}
    for package, placement in graftm.load_placements():
        lineages.setdefault(package, []).append(placement.split("; "))
    for package, chains in lineages.items():
        refpkg = directory / f"{package}.gpkg" / f"{package}.gpkg.refpkg"
        refpkg.mkdir(parents=True)
        (refpkg.parent / "CONTENTS.json").write_text("{}")
        rows, ids = [], {}
        for chain in sorted(chains, key=len):
            tax_id = "|".join(chain)
            parent = "|".join(chain[:-1]) or tax_id
            ids[tax_id] = parent
            rows.append([tax_id, parent, chain[-1]])
        with open(refpkg / f"{package}_taxonomy.csv", "w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(["tax_id", "parent_id", "rank", "tax_name"])
            for tax_id, parent, name in rows:
                writer.writerow([tax_id, parent, "rank", name])
    return directory


def install_fake_tools(bindir):
    bindir = Path(bindir)
    bindir.mkdir(parents=True, exist_ok=True)
    for name in graftm.EXECUTABLES:
        path = bindir / name
        path.write_text(FAKE_GRAFTM)
        path.chmod(0o755)
    return bindir


class TestPlacementTable(unittest.TestCase):
    def setUp(self):
        self.table = graftm.load_placements()

    def test_every_package_and_marker_is_covered(self):
        self.assertEqual({package for package, _ in self.table}, set(graftm.PACKAGES))
        markers = {marker for marker, _ in self.table.values() if marker}
        # norZ is now reachable (qNor); every encoded marker is reachable.
        self.assertEqual(markers, set(graftm.MARKER_FEATURES))
        for package in graftm.PACKAGES:
            self.assertEqual(self.table[package, "Root"], ("", "root_placement"))

    def test_agrees_with_temporal_paper_lookup_except_documented_fixes(self):
        # tests/fixtures/graftm_reference_lookup.tsv is the 69-row tribble from
        # temporal_paper_202505/GraftM-analysis/update_annotations.R.
        with open(FIXTURES / "graftm_reference_lookup.tsv") as handle:
            reference = dict(row for row in csv.reader(handle, delimiter="\t"))
        self.assertEqual(len(reference), 69)
        differences = {}
        for (package, placement), (marker, reason) in self.table.items():
            if placement in reference and reference[placement] != marker:
                differences[package, placement] = (reference[placement], marker or reason)
        expected = {(package, "Root"): ("norB", "root_placement") for package in graftm.PACKAGES}
        for clade in ("qNor_Alphaproteobacteria", "qNor_Mixed", "qNor_Mixed_2", "qNor_Planctomycetota"):
            expected["norZ_norB", f"Root; {clade}"] = ("norB", "norZ")
        self.assertEqual(differences, expected)
        placements = {placement for _, placement in self.table}
        self.assertTrue(set(reference) <= placements)

    def test_annoguild_substring_bugs_are_fixed(self):
        self.assertEqual(self.table["dsrA", "Root; rdsrA; Alphaproteobacteria"][0], "rdsrA")
        self.assertEqual(self.table["dsrB", "Root; rdsrB; Gammaproteobacteria_2"][0], "rdsrB")
        self.assertEqual(self.table["mmoX_sMMO", "Root; Homologues_mmoX"], ("", "outgroup_homologue"))
        self.assertEqual(self.table["norZ_norB", "Root; Nod"], ("", "nitric_oxide_dismutase"))
        self.assertEqual(self.table["mmoX_sMMO", "Root; likely_mmoX; Rhizobiales_mmoX; HYP_mmoX"][0], "mmoX")

    def test_interpret(self):
        self.assertEqual(graftm.interpret("amoA_pmoA", "Root; amoA_Mixed", self.table)["feature"], "K109441")
        self.assertEqual(graftm.interpret("amoA_pmoA", "Root; Nonexistent", self.table)["unresolved_reason"],
                         "placement_not_in_lookup")
        self.assertEqual(graftm.interpret("narG_nxrA", "Root; unknown_3", self.table)["status"], "unresolved")


class TestGraftmRun(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.packages = fake_packages(self.root / "db")
        self.bin = install_fake_tools(self.root / "bin")
        self.fasta = self.root / "proteins.faa"
        self.fasta.write_text(">001 desc\nMAAA\n>g2\nMCCC\n>g3\nMDDD\n>g4\nMEEE\n>g5\nMFFF\n")
        self.annotations = self.root / "kolach_annotations.tsv"
        self.annotations.write_text(
            "gene_id\taccepted_ko\talternative_kos\n"
            "001\tK10944\t-\n"        # accepted pmoA/amoA KO -> amoA
            "g2\tK00001\tK11180\n"   # dsrA only as an alternative -> rdsrA
            "g3\tK04561\t-\n"        # norB KO, placed in Nod -> unresolved
            "g4\tK16157\t-\n"        # mmoX KO, no GraftM hit -> not_placed
            "g5\t-\t-\n")             # no gate KO: never placed
        self.placements = {"amoA_pmoA": {"MAAA": "Root; amoA_Mixed"},
                           "dsrA": {"MCCC": "Root; rdsrA; Chlorobi"},
                           "norZ_norB": {"MDDD": "Root; Nod"}}

    def tearDown(self):
        self.tmp.cleanup()

    def run_graftm(self):
        env = {"PATH": f"{self.bin}{os.pathsep}{os.environ['PATH']}", "GRAFTM_PLACEMENTS": json.dumps(self.placements)}
        out = self.root / "out" / "kolach_graftm_markers.tsv"
        out.parent.mkdir(exist_ok=True)
        with patch.dict(os.environ, env):
            graftm.check_dependencies()
            graftm.annotate(self.fasta, self.annotations, self.packages, out, threads=2)
        with open(out) as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    def test_gating_placement_and_statuses(self):
        rows = self.run_graftm()
        got = {(row["gene_id"], row["package"]): row for row in rows}
        self.assertEqual(set(got), {("001", "amoA_pmoA"), ("g2", "dsrA"), ("g3", "norZ_norB"), ("g4", "mmoX_sMMO")})
        self.assertEqual((got["001", "amoA_pmoA"]["marker"], got["001", "amoA_pmoA"]["feature"]), ("amoA", "K109441"))
        self.assertEqual(got["001", "amoA_pmoA"]["gate_source"], "accepted")
        self.assertEqual((got["g2", "dsrA"]["marker"], got["g2", "dsrA"]["gate_source"]), ("rdsrA", "alternative"))
        self.assertEqual((got["g3", "norZ_norB"]["status"], got["g3", "norZ_norB"]["unresolved_reason"]),
                         ("excluded", "nitric_oxide_dismutase"))
        self.assertEqual(got["g4", "mmoX_sMMO"]["status"], "not_placed")
        self.assertFalse((self.root / "out" / ".graftm_work").exists())
        self.assertTrue((self.root / "out" / "kolach_graftm_markers.tsv.log").is_file())

    def test_one_function_per_protein(self):
        # dsrA accepted, dsrB alternative: the accepted package wins.
        self.annotations.write_text("gene_id\taccepted_ko\talternative_kos\n"
                                    "001\tK11180\tK11181\n"
                                    "g2\t-\tK00370;K00371\n")  # narG and narH, both alternative
        self.placements = {"dsrA": {"MAAA": "Root; dsrA; Chloroflexota"},
                           "dsrB": {"MAAA": "Root; dsrB; Chloroflexota"},
                           "narG_nxrA": {"MCCC": "Root; NarG_Mixed"},
                           "narH_nxrB": {"MCCC": "Root; NarH_Proteobacteria"}}
        got = {(row["gene_id"], row["package"]): row for row in self.run_graftm()}
        self.assertEqual(got["001", "dsrA"]["status"], "refined")
        self.assertEqual((got["001", "dsrB"]["status"], got["001", "dsrB"]["feature"]), ("superseded", ""))
        self.assertEqual({got["g2", p]["unresolved_reason"] for p in ("narG_nxrA", "narH_nxrB")},
                         {"competing_packages"})

    def test_logs_are_kept_when_graftm_fails(self):
        fail = self.bin / "graftM"
        fail.write_text("#!/bin/sh\necho boom\nexit 3\n")
        with self.assertRaises(Exception):
            self.run_graftm()
        log = self.root / "out" / "kolach_graftm_markers.tsv.log"
        self.assertIn("boom", log.read_text())

    def test_missing_package_clade_is_rejected(self):
        taxonomy = self.packages / "dsrA.gpkg" / "dsrA.gpkg.refpkg" / "dsrA_taxonomy.csv"
        taxonomy.write_text(taxonomy.read_text() + "Root|New,Root,rank,New\n")
        with self.assertRaisesRegex(ValueError, "unmapped"):
            graftm.validate_database(self.packages)

    def test_gate_requires_alternative_column(self):
        self.annotations.write_text("gene_id\taccepted_ko\n001\tK10944\n")
        with self.assertRaisesRegex(ValueError, "alternative_kos"):
            graftm.gate_genes(self.annotations)


class TestBundle(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.packages = fake_packages(self.root / "src")
        self.sources = {name: self.packages / f"{name}.gpkg" for name in graftm.PACKAGES}

    def tearDown(self):
        self.tmp.cleanup()

    def test_reproducible_bundle_download_and_checksum(self):
        first = graftm.build_bundle(self.sources, self.root / "a.tar.gz")
        second = graftm.build_bundle(self.sources, self.root / "b.tar.gz")
        self.assertEqual(first, second)
        url = (self.root / "a.tar.gz").as_uri()
        target = self.root / "installed"
        graftm.download_graftm(target, url=url, sha256=first)
        self.assertTrue((target / ".download_complete").exists())
        self.assertEqual(json.loads((target / "manifest.json").read_text())["sha256"], first)
        with self.assertRaisesRegex(ValueError, "checksum"):
            graftm.download_graftm(self.root / "bad", url=url, sha256="0" * 64)
        self.assertFalse((self.root / "bad" / ".download_complete").exists())

    def test_unsafe_archive_is_rejected(self):
        archive = self.root / "evil.tar.gz"
        with tarfile.open(archive, "w:gz") as handle:
            info = tarfile.TarInfo("../escape.txt")
            info.size = 1
            handle.addfile(info, io.BytesIO(b"x"))
        with self.assertRaisesRegex(ValueError, "Unsafe"):
            graftm.install_bundle(archive, self.root / "target")
        self.assertFalse((self.root / "escape.txt").exists())

    def test_symlinks_are_not_bundled(self):
        link = self.packages / "dsrA.gpkg" / "link"
        try:
            link.symlink_to(self.packages)
        except OSError:
            self.skipTest("symlinks unsupported")
        with self.assertRaisesRegex(ValueError, "Symlinks"):
            graftm.build_bundle(self.sources, self.root / "c.tar.gz")


if __name__ == "__main__":
    unittest.main()
