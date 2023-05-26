from tqdm import tqdm
from urllib.request import urlretrieve
from pathlib import Path

def _tqdm_urlretrieve_hook(t: tqdm):
    """From https://github.com/tqdm/tqdm/blob/master/examples/tqdm_wget.py ."""
    last_b = [0]

    def update_to(b=1, bsize=1, tsize=None):
        """
        b  : int, optional
            Number of blocks transferred so far [default: 1].
        bsize  : int, optional
            Size of each block (in tqdm units) [default: 1].
        tsize  : int, optional
            Total size (in tqdm units). If [default: None] or -1,
            remains unchanged.
        """
        if tsize not in (None, -1):
            t.total = tsize
        displayed = t.update((b - last_b[0]) * bsize)
        last_b[0] = b
        return displayed

    return update_to

def cache_files(file_paths: list, cache_dir: str, remote_prefix: str):
    for url in file_paths:
        out_path = Path(cache_dir) / url.removeprefix(remote_prefix).lstrip('/')
        out_path.parent.mkdir(parents=True, exist_ok=True)
        if not out_path.exists():
            with tqdm(unit='B', unit_scale=True, unit_divisor=1024, miniters=1, desc=out_path.name) as t:
                urlretrieve(url, out_path.absolute(), reporthook=_tqdm_urlretrieve_hook(t))
