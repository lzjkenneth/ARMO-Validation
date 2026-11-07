"""Read small windows of NISAR HDF5 products straight from ASF over HTTP range requests (no full download).
Auth: Earthdata login in a netrc file readable only by you (path in EARTHDATA_NETRC, default ~/.netrc; never printed). The signed S3 URL that ASF redirects to
supports byte ranges, so h5py only fetches the chunks it needs (about 1-10 MB per window instead of 1-7 GB per file).
Run reads inside multiprocessing workers: h5py on a Python file object can crash at interpreter exit otherwise."""
import io, os, subprocess, tempfile, time, urllib.request
NETRC = os.environ.get("EARTHDATA_NETRC", os.path.expanduser("~/.netrc"))      # Earthdata login; never printed
COOKIES = os.path.join(os.environ.get("ARMO_SCRATCH", os.path.join(tempfile.gettempdir(), "armo_nisar")), "ed_cookies")

class HTTPRange(io.RawIOBase):
    def __init__(self, url, block=1 << 20):
        self.url, self.block, self.pos, self.cache, self.nread = url, block, 0, {}, 0
        r = urllib.request.urlopen(urllib.request.Request(url, headers={"Range": "bytes=0-0"}), timeout=60)
        self.size = int(r.headers["Content-Range"].split("/")[-1])
    def readable(self): return True
    def seekable(self): return True
    def tell(self): return self.pos
    def seek(self, off, whence=0):
        self.pos = off if whence == 0 else self.pos + off if whence == 1 else self.size + off; return self.pos
    def _blk(self, i):
        if i not in self.cache:
            a = i * self.block; b = min(self.size, a + self.block) - 1
            for k in range(5):
                try:
                    self.cache[i] = urllib.request.urlopen(urllib.request.Request(self.url, headers={"Range": f"bytes={a}-{b}"}), timeout=120).read(); break
                except Exception:
                    if k == 4: raise
                    time.sleep(2 + 3 * k)
            self.nread += len(self.cache[i])
        return self.cache[i]
    def readinto(self, buf):
        n = min(len(buf), self.size - self.pos)
        if n <= 0: return 0
        out = bytearray(); p = self.pos
        while len(out) < n:
            i = p // self.block; o = p - i * self.block; chunk = self._blk(i)[o:o + n - len(out)]; out += chunk; p += len(chunk)
        buf[:n] = out; self.pos += n; return n

def signed_url(collection, granule):
    os.makedirs(os.path.dirname(COOKIES), exist_ok=True)
    u = f"https://nisar.asf.earthdatacloud.nasa.gov/NISAR/{collection}/{granule}/{granule}.h5"
    r = subprocess.run(["curl", "-s", "-L", "--netrc-file", NETRC, "-c", COOKIES, "-b", COOKIES, "-r", "0-0", "-o", "/dev/null",
                        "-w", "%{url_effective}", u], capture_output=True, text=True, timeout=120)
    return r.stdout.strip()
