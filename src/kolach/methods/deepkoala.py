"""DeepKOALA search and annotation wrapper for kolach.

This module delegates inference and processing directly to the official
DeepKOALA package (zhaoxi120/deepkoala), managing database downloads from
GenomeNet and linking them for upstream model loading.
"""
import csv
import hashlib
import json
from pathlib import Path
import re
import shutil
import tempfile
from urllib.error import HTTPError
from urllib.request import urlopen

REFERENCE = "https://github.com/zhaoxi120/deepkoala"
BASE_URL = "https://www.genome.jp/ftp/db/deepkoala/"


def get_available_releases():
    """Retrieve available release dates (YYYYMM) from GenomeNet."""
    pattern = re.compile(r'href="(\d{6})/"')
    with urlopen(BASE_URL, timeout=60) as response:
        html = response.read().decode("utf-8")
    releases = sorted(set(pattern.findall(html)))
    if not releases:
        raise ValueError("Could not find any DeepKOALA releases at GenomeNet.")
    return releases


def setup_resources_link(database_dir, release):
    """Ensure deepkoala's expected resources directory points to our downloaded models."""
    try:
        import deepkoala
        deepkoala_pkg_dir = Path(deepkoala.__file__).resolve().parent
        resources_dir = deepkoala_pkg_dir.parent / "resources"
        resources_dir.mkdir(parents=True, exist_ok=True)

        source_dir = Path(database_dir) / release
        target_link = resources_dir / release

        if target_link.is_symlink() or target_link.exists():
            if target_link.resolve() != source_dir.resolve():
                if target_link.is_symlink():
                    target_link.unlink()
                else:
                    shutil.rmtree(target_link)
                target_link.symlink_to(source_dir, target_is_directory=True)
        else:
            target_link.symlink_to(source_dir, target_is_directory=True)
    except Exception:
        # Fallback if site-packages is read-only; infer.py may already find it or symlink existed
        pass


def download_file(url, destination):
    """Download a file with SHA256 calculation."""
    digest = hashlib.sha256()
    destination = Path(destination)
    with urlopen(url, timeout=180) as response, open(destination, "wb") as out:
        while block := response.read(1024 * 1024):
            out.write(block)
            digest.update(block)
    return digest.hexdigest()


def download_deepkoala(directory, release="latest"):
    """Download DeepKOALA model weights and configs from GenomeNet."""
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)

    available = get_available_releases()
    if release == "latest":
        selected_release = available[-1]
    elif release in available:
        selected_release = release
    else:
        raise ValueError(
            f"Release '{release}' not found on GenomeNet. Available releases: {', '.join(available)}"
        )

    release_dir = directory / selected_release
    release_dir.mkdir(parents=True, exist_ok=True)
    release_url = f"{BASE_URL}{selected_release}/"

    manifest = {
        "method_reference": REFERENCE,
        "release": selected_release,
        "downloads": {}
    }

    # Core files to fetch
    file_specs = [
        ("evaluation.csv", ["evaluation.csv"]),
        ("ko_config_full.json", ["ko_config_full.json"]),
        ("weights_full.pt", ["weights_full.pt"]),
        ("ko_config_frag.json", ["ko_config_fragment.json", "ko_config_frag.json"]),
        ("weights_frag.pt", ["weights_fragment.pt", "weights_frag.pt"]),
    ]

    for target_name, candidates in file_specs:
        downloaded = False
        target_path = release_dir / target_name
        for cand in candidates:
            url = f"{release_url}{cand}"
            try:
                sha256 = download_file(url, target_path)
                manifest["downloads"][target_name] = {"url": url, "sha256": sha256}
                downloaded = True

                # If downloaded as frag/fragment, ensure both aliases exist
                if "frag" in target_name:
                    alt_name = target_name.replace("_frag", "_fragment")
                    alt_path = release_dir / alt_name
                    if not alt_path.exists():
                        try:
                            alt_path.symlink_to(target_path.name)
                        except Exception:
                            shutil.copy2(target_path, alt_path)
                break
            except HTTPError as e:
                if e.code == 404:
                    continue
                raise

        if not downloaded:
            raise RuntimeError(f"Failed to download required DeepKOALA file: {target_name} from {release_url}")

    (release_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    setup_resources_link(directory, selected_release)


def annotate(args, output_file):
    """Run upstream DeepKOALA inference directly."""
    try:
        from deepkoala.infer import inference
    except ImportError as e:
        raise ImportError(
            "The 'deepkoala' package is not installed. "
            "Please install it using: pip install git+https://github.com/zhaoxi120/deepkoala.git"
        ) from e

    db_dir = Path(args.database_dir).expanduser().resolve() / "deepkoala"
    release = getattr(args, "deepkoala_release", "latest")

    # If release is latest, find highest YYYYMM in db_dir
    pat = re.compile(r"^\d{6}$")
    local_releases = [p.name for p in db_dir.iterdir() if p.is_dir() and pat.match(p.name)]
    if not local_releases:
        raise FileNotFoundError(
            f"No DeepKOALA model directory found in {db_dir}. Run 'kolach download --databases deepkoala' first."
        )

    resolved_release = max(local_releases) if release == "latest" else release
    setup_resources_link(db_dir, resolved_release)

    output_path = Path(output_file).expanduser().resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    threads = getattr(args, "threads", 1)
    num_workers = min(max(threads, 1), 4)

    with tempfile.NamedTemporaryFile(prefix="deepkoala_out_", suffix=".csv", delete=False) as tmp:
        tmp_csv = tmp.name

    try:
        inference(
            input_path=str(Path(args.protein_fasta).expanduser().resolve()),
            output_path=tmp_csv,
            model=getattr(args, "deepkoala_model", "full"),
            date=resolved_release,
            batch_size=getattr(args, "batch_size", 64),
            num_workers=num_workers,
            detail=getattr(args, "detail", True),
            device=getattr(args, "device", "auto"),
        )

        # Convert CSV product to TSV format for kolach standard
        with open(tmp_csv, "r", encoding="utf-8") as in_f, open(output_path, "w", encoding="utf-8", newline="") as out_f:
            reader = csv.reader(in_f)
            writer = csv.writer(out_f, delimiter="\t", lineterminator="\n")
            for row in reader:
                writer.writerow(row)

    finally:
        Path(tmp_csv).unlink(missing_ok=True)
