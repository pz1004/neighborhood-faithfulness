#!/usr/bin/env python3
"""Fetch only the preprocessed TEST-set file of THINGS-EEG2 (Gifford et al. 2022) for each
subject, without downloading the training set.

OSF (node anp5v, "preprocessed EEG data", 17 posterior channels) distributes one zip per
subject, about 1.07 GB, holding both preprocessed_eeg_training.npy and
preprocessed_eeg_test.npy. This script reads each zip's central directory with HTTP range
requests, then downloads only the test member's bytes, checks the CRC-32 recorded in the
zip, and writes data/things_eeg2/sub-XX/preprocessed_eeg_test.npy. Nothing is executed or
unpickled here; loading is done separately by a restricted unpickler (reach/things.py).

    python3 fetch_things_eeg2_test.py --list          # central directories only (a few KB)
    python3 fetch_things_eeg2_test.py [--subjects 1,2,...]
"""
import argparse, hashlib, json, os, struct, sys, time, urllib.request, zlib

OSF_FILES = "https://api.osf.io/v2/nodes/anp5v/files/osfstorage/?page[size]=100"
OUT = "data/things_eeg2"
MEMBER = "preprocessed_eeg_test.npy"


def http_range(url, start, end, tries=5):
    for t in range(tries):
        try:
            req = urllib.request.Request(url, headers={"Range": f"bytes={start}-{end}"})
            with urllib.request.urlopen(req, timeout=300) as r:
                if r.status != 206:
                    raise IOError(f"expected 206, got {r.status}")
                data = r.read()
            if len(data) != end - start + 1:
                raise IOError(f"short read {len(data)} of {end - start + 1}")
            return data
        except Exception as e:                       # transient network errors: retry
            if t == tries - 1:
                raise
            print(f"   retry {t + 1} after {e}", flush=True)
            time.sleep(5 * (t + 1))


def subject_zips():
    with urllib.request.urlopen(OSF_FILES, timeout=60) as r:
        d = json.load(r)
    z = {}
    for it in d["data"]:
        a = it["attributes"]
        name = a["name"]
        if a["kind"] == "file" and name.endswith(".zip") and "63_channels" not in name:
            z[int(name[4:6])] = dict(name=name, size=a["size"], url=it["links"]["download"])
    return dict(sorted(z.items()))


def central_directory(url, size):
    tail = http_range(url, max(0, size - 65_557), size - 1)
    k = tail.rfind(b"PK\x05\x06")
    if k < 0:
        raise IOError("no end-of-central-directory record")
    _, _, _, _, n_ent, cd_size, cd_off, _ = struct.unpack("<4sHHHHIIH", tail[k:k + 22])
    cd = http_range(url, cd_off, cd_off + cd_size - 1)
    members, p = [], 0
    for _ in range(n_ent):
        (sig, _, _, flag, method, _, _, crc, csize, usize, nlen, xlen, clen, _, _, _,
         off) = struct.unpack("<4sHHHHHHIIIHHHHHII", cd[p:p + 46])
        assert sig == b"PK\x01\x02"
        name = cd[p + 46:p + 46 + nlen].decode()
        members.append(dict(name=name, method=method, crc=crc, csize=csize, usize=usize,
                            offset=off, flag=flag))
        p += 46 + nlen + xlen + clen
    return members


def fetch_member(url, m, dest):
    head = http_range(url, m["offset"], m["offset"] + 29)
    sig, _, _, method, _, _, _, _, _, nlen, xlen = struct.unpack("<4sHHHHHIIIHH", head)
    assert sig == b"PK\x03\x04" and method == m["method"]
    start = m["offset"] + 30 + nlen + xlen
    raw = bytearray()
    step = 64 * 1024 * 1024
    for s in range(start, start + m["csize"], step):
        raw += http_range(url, s, min(s + step, start + m["csize"]) - 1)
    if m["method"] == 0:
        data = bytes(raw)
    elif m["method"] == 8:
        data = zlib.decompress(bytes(raw), -15)
    else:
        raise IOError(f"unsupported compression method {m['method']}")
    if len(data) != m["usize"] or (zlib.crc32(data) & 0xFFFFFFFF) != m["crc"]:
        raise IOError("size or CRC-32 mismatch")
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    with open(dest + ".part", "wb") as f:
        f.write(data)
    os.replace(dest + ".part", dest)
    return hashlib.sha256(data).hexdigest(), len(raw)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--subjects", default="")
    a = ap.parse_args()
    zips = subject_zips()
    want = [int(s) for s in a.subjects.split(",") if s] or list(zips)
    manifest = []
    for s in want:
        z = zips[s]
        members = central_directory(z["url"], z["size"])
        test = [m for m in members if m["name"].endswith(MEMBER)]
        print(f"sub-{s:02d}: {z['name']} {z['size'] / 1e9:.2f} GB; members: " +
              "; ".join(f"{m['name']} ({m['usize'] / 1e6:.1f} MB, method {m['method']})"
                        for m in members if not m["name"].endswith("/")), flush=True)
        if a.list:
            continue
        assert len(test) == 1, "expected exactly one test member"
        dest = os.path.join(OUT, f"sub-{s:02d}", MEMBER)
        t0 = time.time()
        sha, nbytes = fetch_member(z["url"], test[0], dest)
        print(f"   fetched {nbytes / 1e6:.1f} MB in {time.time() - t0:.0f}s, CRC-32 ok, "
              f"sha256 {sha}", flush=True)
        manifest.append(dict(subject=s, zip=z["name"], zip_size=z["size"], member=test[0]["name"],
                             bytes=test[0]["usize"], crc32=test[0]["crc"], sha256=sha,
                             source=z["url"], fetched=time.strftime("%Y-%m-%d %H:%M:%S")))
    if manifest:
        os.makedirs(OUT, exist_ok=True)
        json.dump(manifest, open(os.path.join(OUT, "manifest.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
