"""Build the reproducible GraftM package bundle published for `kolach download`.

Usage:
  python scripts/build_graftm_bundle.py --temporal-dir <temporal_paper_202505/GraftM-analysis> \
      --mmox-package <mmoX_sMMO.gpkg> --output kolach-graftm-packages-v1.tar.gz

Sources:
- Eight EMERGE packages from the temporal paper (GraftM-analysis/<name>/<name>.gpkg).
- mmoX_sMMO from the public GraftM package collection
  (https://data.ace.uq.edu.au/public/graftm/7/7.60.mmoX_sMMO.gpkg.tar.gz), as used by AnnoGuild.

The bundle is byte-for-byte reproducible from identical inputs. The script
validates every package's clades against the placement table, then prints the
SHA256 to record in `kolach.methods.graftm.BUNDLE_SHA256`.
"""
import argparse
from pathlib import Path
import tempfile

from kolach.methods.graftm import PACKAGES, build_bundle, install_bundle


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--temporal-dir", type=Path, required=True)
    parser.add_argument("--mmox-package", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    sources = {name: args.temporal_dir / name / f"{name}.gpkg" for name in PACKAGES if name != "mmoX_sMMO"}
    sources["mmoX_sMMO"] = args.mmox_package
    digest = build_bundle(sources, args.output)
    with tempfile.TemporaryDirectory() as tmp:
        install_bundle(args.output, tmp)  # Validates clades against the placement table.
    print(f"{args.output}\t{digest}")


if __name__ == "__main__":
    main()
