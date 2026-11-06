"""Step 1: download the 762 genome assemblies from the University of Liverpool dataset.

Source: https://datacat.liverpool.ac.uk/3008/  (licence CC BY 4.0), file genomes.zip (about 1.2 GB).
Each genome is fetched from inside the zip on its own (HTTP byte range) and checked against the zip's own
size and CRC32 checksum, so a dropped connection only costs one file. Genomes already on disk are skipped.

Writes: genomes/<id>.fna   (762 files, about 3.9 GB; not stored in the repository)
Usage:  python scripts/01_download_genomes.py [--jobs 10] [--only ID ID ...]
"""
import argparse, io, os, struct, threading, time, urllib.request, zipfile, zlib
from concurrent.futures import ThreadPoolExecutor

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
URL = "https://datacat.liverpool.ac.uk/3008/7/genomes.zip"
OUT = os.path.join(ROOT, "genomes")


def get(url, headers=None, timeout=120, tries=30):
    for attempt in range(tries):
        try:
            return urllib.request.urlopen(urllib.request.Request(url, headers=headers or {}), timeout=timeout).read()
        except Exception as ex:
            print("retry", attempt + 1, ex, flush=True); time.sleep(min(30, 2 + attempt))
    raise IOError("download failed: " + url)


def zip_entries(size):
    """Read the zip's table of contents from the last 400 kB of the file."""
    tail = 400_000
    data = get(URL, {"Range": f"bytes={size - tail}-{size - 1}"}, timeout=90)

    class Tail(io.RawIOBase):          # a file-like object that only holds the end of the zip
        def __init__(s): s.pos = 0
        def seekable(s): return True
        def readable(s): return True
        def tell(s): return s.pos
        def seek(s, off, wh=0):
            s.pos = off if wh == 0 else s.pos + off if wh == 1 else size + off; return s.pos
        def read(s, n=-1):
            if n < 0: n = size - s.pos
            start = s.pos - (size - tail)
            if start < 0: raise IOError("outside the downloaded tail")
            out = data[start:start + n]; s.pos += len(out); return out

    z = zipfile.ZipFile(Tail())
    return [(os.path.basename(i.filename)[:-4], i.file_size, i.compress_size, i.header_offset, i.CRC)
            for i in z.infolist() if i.filename.endswith(".fna") and not i.filename.startswith("__MACOSX")]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", type=int, default=10, help="parallel downloads")
    ap.add_argument("--only", nargs="*", help="download only these genome IDs (for testing)")
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    size = int(urllib.request.urlopen(urllib.request.Request(URL, method="HEAD"), timeout=60).headers["Content-Length"])
    entries = zip_entries(size)
    print("genomes in the zip:", len(entries), flush=True)
    if args.only: entries = [e for e in entries if e[0] in set(args.only)]
    todo = [e for e in entries if not os.path.exists(os.path.join(OUT, e[0] + ".fna"))]
    print("to download:", len(todo), "| compressed MB:", round(sum(e[2] for e in todo) / 1e6), flush=True)
    lock = threading.Lock(); done = [0]; failed = []

    def fetch(e):
        gid, usize, csize, off, crc = e
        dest = os.path.join(OUT, gid + ".fna")
        for attempt in range(60):
            try:
                req = urllib.request.Request(URL, headers={"Range": f"bytes={off}-{min(size - 1, off + 30 + 400 + csize + 64)}"})
                buf = urllib.request.urlopen(req, timeout=120).read()
                nlen, xlen = struct.unpack("<HH", buf[26:30])          # local file header: name and extra-field lengths
                raw = buf[30 + nlen + xlen: 30 + nlen + xlen + csize]
                if len(raw) != csize: raise IOError("short read")
                txt = zlib.decompress(raw, -15)
                if len(txt) != usize or zlib.crc32(txt) != crc: raise IOError("size/CRC mismatch")
                open(dest + ".part", "wb").write(txt); os.replace(dest + ".part", dest)
                with lock:
                    done[0] += 1
                    if done[0] % 10 == 0 or done[0] == len(todo): print(f"{done[0]}/{len(todo)} downloaded", flush=True)
                return
            except Exception:
                time.sleep(min(30, 2 + attempt))
        failed.append(gid); print("FAILED", gid, flush=True)

    with ThreadPoolExecutor(args.jobs) as ex: list(ex.map(fetch, todo))
    have = sum(os.path.exists(os.path.join(OUT, e[0] + ".fna")) for e in entries)
    print(f"finished: {have}/{len(entries)} genomes on disk", "| failed: " + ", ".join(failed) if failed else "", flush=True)


if __name__ == "__main__":
    main()
