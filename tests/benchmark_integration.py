"""Run old/new integration in fresh processes; report time and process peak RSS.

    PYTHONPATH=src python tests/benchmark_integration.py --genes 20000 --repeats 3

The larger fixture is deterministic and exercises all three tools. This script
also verifies full DataFrame and TSV parity for that fixture before measuring.
"""
import argparse
import ctypes
import json
from pathlib import Path
import statistics
import subprocess
import sys
import tempfile
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))


def peak_rss_mib():
    if sys.platform != "win32":
        import resource
        peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        return peak / (1024 * 1024) if sys.platform == "darwin" else peak / 1024

    from ctypes import wintypes

    class MemoryCounters(ctypes.Structure):
        _fields_ = [
            ("cb", wintypes.DWORD), ("PageFaultCount", wintypes.DWORD),
            ("PeakWorkingSetSize", ctypes.c_size_t), ("WorkingSetSize", ctypes.c_size_t),
            ("QuotaPeakPagedPoolUsage", ctypes.c_size_t), ("QuotaPagedPoolUsage", ctypes.c_size_t),
            ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t), ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
            ("PagefileUsage", ctypes.c_size_t), ("PeakPagefileUsage", ctypes.c_size_t),
        ]

    kernel = ctypes.WinDLL("kernel32", use_last_error=True)
    kernel.GetCurrentProcess.restype = wintypes.HANDLE
    psapi = ctypes.WinDLL("psapi", use_last_error=True)
    psapi.GetProcessMemoryInfo.argtypes = [wintypes.HANDLE, ctypes.POINTER(MemoryCounters), wintypes.DWORD]
    psapi.GetProcessMemoryInfo.restype = wintypes.BOOL
    counters = MemoryCounters()
    counters.cb = ctypes.sizeof(counters)
    if not psapi.GetProcessMemoryInfo(kernel.GetCurrentProcess(), ctypes.byref(counters), counters.cb):
        raise ctypes.WinError(ctypes.get_last_error())
    return counters.PeakWorkingSetSize / (1024 * 1024)


def fixture(directory, genes):
    paths = {tool: directory / (tool + "_annotations.tsv") for tool in ("kofam", "deepkoala", "eggnog")}
    with paths["kofam"].open("w") as kofam, paths["deepkoala"].open("w") as deepkoala, paths["eggnog"].open("w") as eggnog:
        kofam.write("gene_id\tko\tassignment\tbit_score\te_value\tthreshold\tscore_type\n")
        deepkoala.write("name\tpredict_label\tprobability\tthreshold\tannotate\n")
        eggnog.write("#query\tkegg_ko\tscore\tevalue\n")
        for index in range(genes):
            gene = f"gene_{index:06}"
            ko = f"K{index % 500 + 1:05}"
            other = f"K{(index + 1) % 500 + 1:05}"
            assignment = ("threshold", "rescued", "below")[index % 3]
            kofam.write(f"{gene}\t{ko}\t{assignment}\t100\t1e-10\t80\tfull\n")
            deepkoala.write(f"{gene}\t{ko if index % 4 else other}\t{0.9 if index % 2 else 0.4}\t0.5\t-\n")
            eggnog.write(f"{gene}\t{ko},{other}\t{100 if index % 5 else 40}\t0\n")
            if index % 7 == 0:
                eggnog.write(f"{gene}\t{other}\t90\t1e-9\n")
    return paths


def worker(version, directory):
    from kolach import integrate
    from tests.reference import integrate_db5a5ce
    module = integrate if version == "new" else integrate_db5a5ce
    paths = {tool + "_tsv": directory / (tool + "_annotations.tsv") for tool in ("kofam", "deepkoala", "eggnog")}
    start = time.perf_counter()
    result = module.integrate_annotations(**paths)
    elapsed = time.perf_counter() - start
    print(json.dumps({"seconds": elapsed, "peak_rss_mib": peak_rss_mib(), "genes": len(result)}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genes", type=int, default=20000)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--worker", choices=("old", "new"))
    parser.add_argument("--directory", type=Path)
    args = parser.parse_args()
    if args.worker:
        worker(args.worker, args.directory)
        return
    if args.genes < 1 or args.repeats < 1:
        parser.error("--genes and --repeats must be positive")

    from tests.test_integration_parity import PipelineParity

    report = {}
    with tempfile.TemporaryDirectory() as directory:
        directory = Path(directory)
        paths = fixture(directory, args.genes)
        PipelineParity().compare(**{tool + "_tsv": path for tool, path in paths.items()})
        datasets = [("examples", ROOT / "examples"), ("synthetic", directory)]
        for label, dataset in datasets:
            samples = {"old": [], "new": []}
            # Alternate versions to reduce bias from machine load or caching.
            for repeat in range(args.repeats):
                order = ("old", "new") if repeat % 2 == 0 else ("new", "old")
                for version in order:
                    command = [sys.executable, __file__, "--worker", version, "--directory", str(dataset)]
                    completed = subprocess.run(command, capture_output=True, text=True, check=True)
                    samples[version].append(json.loads(completed.stdout))
            report[label] = {
                version: {
                    "median_seconds": statistics.median(sample["seconds"] for sample in values),
                    "median_peak_rss_mib": statistics.median(sample["peak_rss_mib"] for sample in values),
                    "genes": values[0]["genes"], "samples": values,
                }
                for version, values in samples.items()
            }
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
