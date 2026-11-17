"""Shared types for organism-specific typing modules."""
from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Reference:
    key: str            # e.g. "Harmsen2003"
    authors: str        # "Harmsen D, Claus H, Witte W, et al."
    year: int
    title: str
    journal: str        # journal, volume:pages
    doi: str = ""
    pmid: str = ""
    supports: str = ""  # which call or rule in the module rests on it


@dataclass(frozen=True)
class OrganismDb:
    name: str            # folder name under db_dir/<KEY>/
    source: str          # URL it is fetched from
    version: str         # release, tag, commit or retrieval date
    licence: str         # licence name, or the terms quoted; "none stated" if so
    redistribution: str  # "permitted" | "restricted" | "unstated"
    reference: str       # Reference.key
    notes: str = ""
