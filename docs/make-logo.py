"""Emit the Hydra emblem: five serpent heads on one body inside a plasmid ring."""
import math, pathlib, sys

BASE = (256, 376)          # where every neck grows from
RING_R, RING_W = 186, 22
ANGLES = [-46, -23, 0, 23, 46]   # wider and the outer heads collide with the ring
LIMB, NECK_W, HEAD_W, HEAD_H = 164, 32, 58, 82
BOW = 14


def head(cx, cy, angle, w, h, eye_side, eye_fill):
    hw, snout, jaw = w / 2, h * 0.58, h * 0.42
    path = (f"M 0 {-snout:.0f} "
            f"C {hw*0.56:.0f} {-snout*0.88:.0f} {hw:.0f} {-snout*0.32:.0f} {hw:.0f} {jaw*0.20:.0f} "
            f"C {hw:.0f} {jaw*0.82:.0f} {hw*0.62:.0f} {jaw:.0f} 0 {jaw:.0f} "
            f"C {-hw*0.62:.0f} {jaw:.0f} {-hw:.0f} {jaw*0.82:.0f} {-hw:.0f} {jaw*0.20:.0f} "
            f"C {-hw:.0f} {-snout*0.32:.0f} {-hw*0.56:.0f} {-snout*0.88:.0f} 0 {-snout:.0f} Z")
    t = f"translate({cx:.0f} {cy:.0f}) rotate({angle})"
    return (f'    <g transform="{t}"><path d="{path}"/></g>\n',
            f'    <g transform="{t}"><ellipse cx="{eye_side*hw*0.28:.0f}" cy="{-snout*0.04:.0f}" '
            f'rx="{w*0.10:.0f}" ry="{w*0.13:.0f}" fill="{eye_fill}"/></g>\n')


def limb(angle_deg, length, bow):
    """A neck radiating from BASE, bowed outward, every control point above it.

    Keeping the curve entirely above the base is what stops the outer necks
    arcing down into a smile - the failure mode of every earlier draft.
    """
    rad = math.radians(angle_deg)
    tipx = BASE[0] + math.sin(rad) * length
    tipy = BASE[1] - math.cos(rad) * length
    midx = BASE[0] + math.sin(rad) * (length * 0.55 + bow)
    midy = BASE[1] - math.cos(rad) * length * 0.55
    return f"M {BASE[0]} {BASE[1]-10} Q {midx:.0f} {midy:.0f} {tipx:.0f} {tipy:.0f}", tipx, tipy


def rod(cx, cy, w, h, angle):
    return (f'    <rect x="{-w/2:.0f}" y="{-h/2:.0f}" width="{w}" height="{h}" rx="{w/2:.0f}" '
            f'transform="translate({cx} {cy}) rotate({angle})"/>\n')


def emblem(bg, mark, eye):
    necks = heads = eyes = ""
    for a in ANGLES:
        d, tx, ty = limb(a, LIMB, BOW if a else 0)
        necks += (f'    <path d="{d}" stroke-width="{NECK_W}" fill="none" '
                  f'stroke="currentColor" stroke-linecap="round"/>\n')
        h, e = head(tx, ty, a, HEAD_W, HEAD_H, 1 if a >= 0 else -1, eye)
        heads += h
        eyes += e

    grad = ('  <defs>\n'
            '    <linearGradient id="tile" x1="0" y1="0" x2="1" y2="1">\n'
            '      <stop offset="0" stop-color="#0d9488"/>\n'
            '      <stop offset="0.5" stop-color="#0ea5e9"/>\n'
            '      <stop offset="1" stop-color="#6366f1"/>\n'
            '    </linearGradient>\n  </defs>\n')
    tile = (f'  <rect width="512" height="512" rx="116" fill="{bg}"/>\n' if bg
            else grad + '  <rect width="512" height="512" rx="116" fill="url(#tile)"/>\n')

    body = (
        "  <!-- The ring is a plasmid: closed circular DNA, which is how most\n"
        "       acquired resistance actually travels between bacteria. -->\n"
        f'    <circle cx="256" cy="256" r="{RING_R}" fill="none" stroke="currentColor" '
        f'stroke-width="{RING_W}"/>\n'
        "  <!-- Necks first, then one body they all grow from, then the heads on\n"
        "       top, so the whole creature reads as a single silhouette. -->\n"
        + necks
        + f'    <ellipse cx="{BASE[0]}" cy="{BASE[1]}" rx="52" ry="24"/>\n'
        + heads
        + "  <!-- Two bacilli, set well inside the ring so they read as cells\n"
          "       rather than as ticks on it. -->\n"
        + rod(166, 352, 22, 52, 26) + rod(346, 352, 22, 52, -26)
        + "  <!-- One eye per head, cut back to the badge colour. Below ~32px they\n"
          "       drop out and the mark degrades to a clean silhouette. -->\n"
        + eyes)

    return ('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 512 512" role="img" '
            'aria-label="Hydra">\n  <title>Hydra</title>\n'
            '  <desc>A five-headed serpent on one body, inside a plasmid ring, '
            'flanked by bacterial rods</desc>\n'
            + tile + f'  <g fill="{mark}" color="{mark}">\n{body}  </g>\n</svg>\n')


if __name__ == "__main__":
    out = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else "logo.svg")
    variant = sys.argv[2] if len(sys.argv) > 2 else "dark"
    if variant == "dark":
        svg = emblem("#0b1220", "#22d3ee", "#0b1220")
    else:
        svg = emblem(None, "#ffffff", "#0ea5e9")
    out.write_text(svg)
    print("wrote", out)
