"""The start-up banner: a six-headed serpent, the wordmark, and who wrote it.

The drawing is the terminal counterpart of ``docs/logo.svg`` - six heads
fanning out of a rod-shaped bacterium - so the tool looks the same in a
terminal as it does in the README.

The banner is decoration, so it always goes to stderr and never to stdout:
``hydra screen --stdout`` pipes a TSV table into other programs, and a banner
in that stream would corrupt it. ``--quiet``, ``HYDRA_NO_BANNER`` and
``--no-banner`` each switch it off; ``NO_COLOR`` and a non-terminal stream drop
the colour but keep the drawing.
"""

from __future__ import annotations

import os
import sys

from . import __author__, __organisation__, __version__

AUTHOR = __author__
ORGANISATION = __organisation__
TAGLINE = "many heads, one genome"
#: One-line attribution for logs, reports and run manifests.
CREDIT = f"{AUTHOR}, {ORGANISATION}"

_LOGO = (
    "     o     o     ",
    "  o   \\   /   o  ",
    "o  \\   \\ /   /  o",
    " \\__\\___|___/__/ ",
    "     [=====]     ",
)

_LETTERS = {
    "H": ("#   #", "#   #", "#####", "#   #", "#   #"),
    "Y": ("#   #", " # # ", "  #  ", "  #  ", "  #  "),
    "D": ("#### ", "#   #", "#   #", "#   #", "#### "),
    "R": ("#### ", "#   #", "#### ", "#  # ", "#   #"),
    "A": (" ### ", "#   #", "#####", "#   #", "#   #"),
}
_WORDMARK = tuple("  ".join(_LETTERS[c][row] for c in "HYDRA") for row in range(5))

_INDENT = "  "
_GAP = "   "
#: Column the two text lines start at, so they hang under the wordmark.
_TEXT_COLUMN = " " * (len(_INDENT) + len(_LOGO[0]) + len(_GAP))

_BOLD = "\033[1m"
_DIM = "\033[2m"
_CYAN = "\033[36m"
_RESET = "\033[0m"


def _supports_colour(stream) -> bool:
    if os.environ.get("NO_COLOR"):
        return False
    if os.environ.get("TERM", "") == "dumb":
        return False
    try:
        return bool(stream.isatty())
    except (AttributeError, ValueError):
        return False


def _separator(stream) -> str:
    """``·`` where the terminal can encode it, ``-`` on legacy code pages."""
    try:
        "·".encode(stream.encoding or "ascii")
    except (LookupError, UnicodeEncodeError, AttributeError):
        return " - "
    return " · "


def render(colour: bool = False, separator: str = " · ") -> str:
    """The banner as a string, with or without ANSI colour."""
    logo, mark = _LOGO, _WORDMARK
    if colour:
        logo = tuple(f"{_CYAN}{line}{_RESET}" for line in logo)
        mark = tuple(f"{_BOLD}{line}{_RESET}" for line in mark)

    lines = [f"{_INDENT}{art}{_GAP}{text}".rstrip()
             for art, text in zip(logo, mark)]

    version = f"v{__version__}"
    credit = f"{AUTHOR}{separator}{ORGANISATION}"
    if colour:
        version = f"{_DIM}{TAGLINE}{separator}{version}{_RESET}"
        credit = f"{_BOLD}{AUTHOR}{_RESET}{_DIM}{separator}{_RESET}{_CYAN}{ORGANISATION}{_RESET}"
    else:
        version = f"{TAGLINE}{separator}{version}"

    lines += ["", _TEXT_COLUMN + version, _TEXT_COLUMN + credit]
    return "\n".join(lines)


def render_for(stream) -> str:
    """The banner, coloured and punctuated for what ``stream`` can display."""
    return render(_supports_colour(stream), _separator(stream))


def show(stream=None, enabled: bool = True) -> None:
    """Print the banner to stderr unless it has been switched off."""
    if not enabled or os.environ.get("HYDRA_NO_BANNER"):
        return
    stream = sys.stderr if stream is None else stream
    print(render_for(stream), end="\n\n", file=stream)
