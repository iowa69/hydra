"""Database registry, installation and lookup."""

from .registry import DATABASES, DbSpec, resolve_names
from .manager import DatabaseStore

__all__ = ["DATABASES", "DbSpec", "DatabaseStore", "resolve_names"]
