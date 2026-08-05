"""Generate the Hydra emblem.

A six-headed serpent growing out of a rod-shaped bacterium, inside a ring that
is a plasmid - closed circular DNA, which is how most acquired resistance
actually moves between bacteria.

Two things make it read as a creature rather than a diagram. The necks are drawn
as filled outlines whose width shrinks toward the head: a constant-width stroke
reads as a shaft, a tapering one reads as a neck. And each neck follows an
S-curve, with its head rotated onto the curve's end tangent, so the heads point
where their necks are going.
"""
import math, pathlib, sys

BASE = (256, 386)
RING_R, RING_W = 186, 22

# angle from vertical, length, sway (sign sets which way the S starts),
# width at the body, width where the head joins
# Mirrored pairs: flipping the angle flips the sway, so the two halves undulate
# outward together instead of drifting the same way.
NECKS = [(-60, 150, 58, 38, 16),
         (-33, 174, -52, 38, 16),
         (-10, 190, -44, 40, 17),
         (10, 190, 44, 40, 17),
         (33, 174, 52, 38, 16),
         (60, 150, -58, 38, 16)]


def _pt(angle_deg, along, across):
    rad = math.radians(angle_deg)
    dx, dy = math.sin(rad), -math.cos(rad)
    px, py = math.cos(rad), math.sin(rad)
    return BASE[0] + dx * along + px * across, BASE[1] + dy * along + py * across


def _bezier(p0, c1, c2, p3, t):
    u = 1 - t
    x = (u**3 * p0[0] + 3 * u * u * t * c1[0] + 3 * u * t * t * c2[0] + t**3 * p3[0])
    y = (u**3 * p0[1] + 3 * u * u * t * c1[1] + 3 * u * t * t * c2[1] + t**3 * p3[1])
    dx = (3 * u * u * (c1[0] - p0[0]) + 6 * u * t * (c2[0] - c1[0])
          + 3 * t * t * (p3[0] - c2[0]))
    dy = (3 * u * u * (c1[1] - p0[1]) + 6 * u * t * (c2[1] - c1[1])
          + 3 * t * t * (p3[1] - c2[1]))
    return (x, y), (dx, dy)


def neck(angle_deg, length, sway, w0, w1, steps=22):
    """A tapered, S-curved neck drawn as one filled outline.

    A stroked path has the same width end to end, which reads as a shaft. Walking
    the curve and offsetting by a width that shrinks toward the head gives the
    taper a snake actually has. The sway is lateral only - progress along the
    axis is monotonic - so no neck can arc back down and turn the mark into a
    face, which sank several earlier drafts.
    """
    p0 = (BASE[0], BASE[1] - 6)
    c1 = _pt(angle_deg, length * 0.32, sway)
    c2 = _pt(angle_deg, length * 0.68, -sway * 0.9)
    p3 = _pt(angle_deg, length, sway * 0.12)

    left, right = [], []
    for i in range(steps + 1):
        t = i / steps
        (x, y), (dx, dy) = _bezier(p0, c1, c2, p3, t)
        norm = math.hypot(dx, dy) or 1.0
        nx, ny = -dy / norm, dx / norm
        half = (w0 + (w1 - w0) * t) / 2
        left.append((x + nx * half, y + ny * half))
        right.append((x - nx * half, y - ny * half))

    pts = left + list(reversed(right))
    d = "M " + " L ".join(f"{x:.0f} {y:.0f}" for x, y in pts) + " Z"
    (tip, tangent) = _bezier(p0, c1, c2, p3, 1.0)
    heading = math.degrees(math.atan2(tangent[0], -tangent[1]))
    return d, tip, heading


def head(cx, cy, angle, w, h, eye_side, eye_fill):
    """Snake head in profile: narrow snout, jaw flaring back on one side."""
    hw, snout, jaw = w / 2, h * 0.64, h * 0.36
    path = (f"M 0 {-snout:.0f} "
            f"C {hw*0.74:.0f} {-snout*0.78:.0f} {hw*1.16:.0f} {-snout*0.16:.0f} "
            f"{hw*1.08:.0f} {jaw*0.45:.0f} "
            f"C {hw*1.00:.0f} {jaw*0.95:.0f} {hw*0.55:.0f} {jaw:.0f} 0 {jaw:.0f} "
            f"C {-hw*0.55:.0f} {jaw:.0f} {-hw*0.94:.0f} {jaw*0.80:.0f} {-hw*0.88:.0f} {jaw*0.32:.0f} "
            f"C {-hw*0.82:.0f} {-snout*0.20:.0f} {-hw*0.64:.0f} {-snout*0.78:.0f} 0 {-snout:.0f} Z")
    t = f"translate({cx:.0f} {cy:.0f}) rotate({angle:.0f})"
    return (f'    <g transform="{t}"><path d="{path}"/></g>\n',
            f'    <g transform="{t}"><ellipse cx="{eye_side*hw*0.44:.0f}" cy="{-snout*0.14:.0f}" '
            f'rx="{w*0.115:.0f}" ry="{w*0.155:.0f}" fill="{eye_fill}"/></g>\n')


def rod(cx, cy, w, h, angle):
    return (f'    <rect x="{-w/2:.0f}" y="{-h/2:.0f}" width="{w}" height="{h}" rx="{w/2:.0f}" '
            f'transform="translate({cx} {cy}) rotate({angle})"/>\n')


def emblem(bg, mark, eye):
    necks = heads = eyes = ""
    for angle, length, sway, w0, w1 in NECKS:
        d, tip, heading = neck(angle, length, sway, w0, w1)
        necks += f'    <path d="{d}"/>\n'
        h, e = head(tip[0], tip[1], heading, 50, 78, 1 if angle >= 0 else -1, eye)
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

    body = (f'    <circle cx="256" cy="256" r="{RING_R}" fill="none" stroke="currentColor" '
            f'stroke-width="{RING_W}"/>\n'
            + necks
            + f'    <rect x="{BASE[0]-78}" y="{BASE[1]-27}" width="156" '
              f'height="54" rx="27"/>\n'
            + heads
            + eyes)

    return ('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 512 512" role="img" '
            'aria-label="Hydra">\n  <title>Hydra</title>\n'
            '  <desc>A six-headed serpent growing from a rod-shaped bacterium, '
            'inside a plasmid ring</desc>\n'
            + tile + f'  <g fill="{mark}" color="{mark}">\n{body}  </g>\n</svg>\n')


if __name__ == "__main__":
    out = pathlib.Path(sys.argv[1])
    variant = sys.argv[2] if len(sys.argv) > 2 else "dark"
    out.write_text(emblem("#0b1220", "#22d3ee", "#0b1220") if variant == "dark"
                   else emblem(None, "#ffffff", "#0ea5e9"))
    print("wrote", out)
