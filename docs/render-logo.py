"""Render an SVG to an exact-size PNG.

Headless Chrome reserves part of --window-size for browser chrome, so the
viewport is shorter than asked for and an exactly-sized image gets clipped.
Measure the reservation once, render into a window tall enough that the
viewport is the size we want, then crop the surplus rows off.
"""
import pathlib, struct, subprocess, sys, tempfile, zlib

CHROME = ["google-chrome", "--headless", "--disable-gpu", "--no-sandbox",
          "--hide-scrollbars", "--default-background-color=00000000",
          "--force-device-scale-factor=1"]


def decode(path):
    data = path.read_bytes()
    i, idat = 8, b""
    while i < len(data):
        length = struct.unpack(">I", data[i:i + 4])[0]
        kind, chunk = data[i + 4:i + 8], data[i + 8:i + 8 + length]
        if kind == b"IHDR":
            w, h, depth, colour = struct.unpack(">IIBB", chunk[:10])
            if depth != 8 or colour != 6:
                raise SystemExit(f"expected 8-bit RGBA, got depth {depth} colour {colour}")
        elif kind == b"IDAT":
            idat += chunk
        i += 12 + length
    raw, stride = zlib.decompress(idat), w * 4
    rows, prev, pos = [], bytearray(stride), 0
    for _ in range(h):
        f = raw[pos]; pos += 1
        line = bytearray(raw[pos:pos + stride]); pos += stride
        for x in range(stride):
            a = line[x - 4] if x >= 4 else 0
            b = prev[x]
            c = prev[x - 4] if x >= 4 else 0
            if f == 1: line[x] = (line[x] + a) & 255
            elif f == 2: line[x] = (line[x] + b) & 255
            elif f == 3: line[x] = (line[x] + (a + b) // 2) & 255
            elif f == 4:
                p = a + b - c
                pa, pb, pc = abs(p - a), abs(p - b), abs(p - c)
                line[x] = (line[x] + (a if pa <= pb and pa <= pc else b if pb <= pc else c)) & 255
        prev = line
        rows.append(bytes(line))
    return w, h, rows


def encode(path, w, h, rows):
    raw = b"".join(b"\x00" + r[:w * 4] for r in rows[:h])
    def chunk(kind, payload):
        return (struct.pack(">I", len(payload)) + kind + payload
                + struct.pack(">I", zlib.crc32(kind + payload) & 0xFFFFFFFF))
    path.write_bytes(b"\x89PNG\r\n\x1a\n"
                     + chunk(b"IHDR", struct.pack(">IIBBBBB", w, h, 8, 6, 0, 0, 0))
                     + chunk(b"IDAT", zlib.compress(raw, 9))
                     + chunk(b"IEND", b""))


def shoot(html, out, win_w, win_h):
    subprocess.run(CHROME + [f"--screenshot={out}", f"--window-size={win_w},{win_h}",
                             f"file://{html}"], check=True, capture_output=True)


def render(svg: pathlib.Path, out: pathlib.Path, size: int):
    with tempfile.TemporaryDirectory() as tmp:
        work = pathlib.Path(tmp)
        (work / "logo.svg").write_bytes(svg.read_bytes())
        html = work / "r.html"
        html.write_text('<style>html,body{margin:0;padding:0}'
                        f'img{{display:block;width:{size}px;height:{size}px}}</style>'
                        '<img src="logo.svg">')
        # Calibrate: how much vertical space does this Chrome keep for itself?
        probe = work / "probe.png"
        shoot(html, probe, size, size)
        _, _, rows = decode(probe)
        opaque = [y for y, r in enumerate(rows) if any(r[x * 4 + 3] > 8 for x in range(size))]
        reserved = size - (opaque[-1] + 1)
        shot = work / "shot.png"
        shoot(html, shot, size, size + reserved)
        w, h, rows = decode(shot)
        if len(rows) < size or w < size:
            raise SystemExit(f"render came back {w}x{h}, too small to crop to {size}")
        encode(out, size, size, rows)
        return reserved


if __name__ == "__main__":
    svg, out, size = pathlib.Path(sys.argv[1]), pathlib.Path(sys.argv[2]), int(sys.argv[3])
    print(f"wrote {out} ({size}x{size}); chrome reserved {render(svg, out, size)}px")
