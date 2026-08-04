"""Fetch reference databases from their upstream sources.

Each fetcher downloads what upstream actually publishes and stages it into the
layout the importers in :mod:`hydra_amr.db.manager` already expect, so
``hydra db download NAME`` and ``hydra db import --source DIR`` converge on one
conversion path rather than two that can drift apart.

Only databases with a stable, directly fetchable source live here. The rest are
distributed as landing pages or repositories with no versioned download URL, and
guessing at one would break silently the first time upstream reorganised; for
those, ``hydra db download`` still prints where to get them by hand.
"""

from __future__ import annotations

import gzip
import json
import re
import shutil
import subprocess
import tarfile
import time
import urllib.error
import urllib.request
from pathlib import Path
from typing import Callable

from ..seqio import read_fasta, write_fasta
from ..utils import LOG, HydraError

#: Upstream is asked to identify us, so operators can see what the traffic is.
USER_AGENT = "hydra-amr (+https://github.com/iowa69/hydra)"
TIMEOUT = 120

NCBI_LATEST = ("https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/"
               "AMRFinderPlus/database/latest")

#: Files the protein reference and the mutation catalogues are built from. The
#: per-organism AMR_DNA-* files are discovered from the directory listing.
NCBI_PROT_FILES = ("AMRProt.fa", "fam.tsv", "AMRProt-mutation.tsv",
                   "AMRProt-suppress.tsv", "AMRProt-susceptible.fa",
                   "AMRProt-susceptible.tsv", "taxgroup.tsv")


def _get(url: str, dest: Path, progress: Callable[[str], None] | None = None,
         attempts: int = 4) -> Path:
    """Download *url* to *dest*, retrying a few times before giving up.

    A PubMLST install is over a thousand requests, so a transient failure is
    likely rather than exceptional, and one dropped locus invalidates a whole
    scheme. Only transport-level failures are retried; a 404 is answered once.
    """
    last: Exception | None = None
    for attempt in range(attempts):
        try:
            return _get_once(url, dest, progress)
        except HydraError as exc:
            last = exc
            throttled = "429" in str(exc)
            if "HTTP Error 4" in str(exc) and not throttled:
                raise
            if attempt + 1 < attempts:
                # Being throttled means backing off properly, not trying again
                # immediately with the same pressure that caused it.
                time.sleep((5.0 if throttled else 1.5) * (attempt + 1))
    raise last if last else HydraError(f"could not download {url}")


def _get_once(url: str, dest: Path, progress: Callable[[str], None] | None = None) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    try:
        with urllib.request.urlopen(request, timeout=TIMEOUT) as response:
            total = int(response.headers.get("Content-Length") or 0)
            done = 0
            with open(dest, "wb") as handle:
                while True:
                    block = response.read(1 << 16)
                    if not block:
                        break
                    handle.write(block)
                    done += len(block)
                    if progress and total:
                        progress(f"{dest.name}: {done * 100 // total}%")
    except urllib.error.HTTPError as exc:
        if exc.code in (403, 406):
            # Some hosts (VFDB is one) allowlist a handful of exact User-Agent
            # strings and refuse everything else. Rather than impersonate one,
            # hand the request to curl or wget, which really are those agents.
            if _get_via_tool(url, dest):
                return dest
        raise HydraError(f"could not download {url}: {exc}") from exc
    except (urllib.error.URLError, OSError) as exc:
        raise HydraError(f"could not download {url}: {exc}") from exc
    if dest.stat().st_size == 0:
        raise HydraError(f"{url} returned an empty file")
    return dest


def _get_via_tool(url: str, dest: Path) -> bool:
    """Retry a refused download through curl or wget, if either is installed."""
    for command in (["curl", "-fsSL", "--max-time", str(TIMEOUT), "-o", str(dest), url],
                    ["wget", "-q", "-T", str(TIMEOUT), "-O", str(dest), url]):
        if shutil.which(command[0]) is None:
            continue
        LOG.debug("retrying %s with %s", url, command[0])
        if (subprocess.run(command, capture_output=True).returncode == 0
                and dest.exists() and dest.stat().st_size > 0):
            return True
    return False


def _listing(url: str) -> list[str]:
    """Filenames linked from an FTP-over-HTTP directory index."""
    request = urllib.request.Request(url + "/", headers={"User-Agent": USER_AGENT})
    try:
        with urllib.request.urlopen(request, timeout=TIMEOUT) as response:
            html = response.read().decode("utf-8", "replace")
    except (urllib.error.URLError, OSError) as exc:
        raise HydraError(f"could not list {url}: {exc}") from exc
    return [href for href in re.findall(r'href="([^"]+)"', html)
            if not href.startswith(("/", "?", "http"))]


def _unpack(archive: Path, into: Path) -> Path:
    into.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive) as tar:
        # Refuse members that would escape the extraction directory.
        for member in tar.getmembers():
            target = (into / member.name).resolve()
            if not str(target).startswith(str(into.resolve())):
                raise HydraError(f"{archive.name} contains an unsafe path: {member.name}")
        # Python 3.12+ wants an explicit filter and 3.14 makes 'data' the default;
        # ask for it where it exists so behaviour does not change under us.
        try:
            tar.extractall(into, filter="data")
        except TypeError:
            tar.extractall(into)
    return into


def _write_sequences(records: list[tuple[str, str]], work: Path) -> Path:
    """Write staged records where :func:`_import_nucl` looks for them."""
    if not records:
        raise HydraError("upstream returned no sequences")
    work.mkdir(parents=True, exist_ok=True)
    write_fasta(work / "sequences", records)
    return work


# --------------------------------------------------------------- NCBI protein
def stage_protein(work: Path, progress=None) -> Path:
    """Mirror the AMRFinderPlus data directory, which imports as-is."""
    files = list(NCBI_PROT_FILES)
    listing = _listing(NCBI_LATEST)
    files += [name for name in listing if name.startswith("AMR_DNA-")]
    work.mkdir(parents=True, exist_ok=True)
    for name in files:
        if name not in listing and name not in NCBI_PROT_FILES:
            continue
        try:
            _get(f"{NCBI_LATEST}/{name}", work / name, progress)
        except HydraError:
            # Optional companions come and go between releases; the importer
            # reports anything genuinely required as missing.
            if name in ("AMRProt.fa", "fam.tsv"):
                raise
            LOG.debug("upstream has no %s in this release", name)
    return work


# ------------------------------------------------------------ NCBI nucleotide
def stage_ncbi(work: Path, progress=None) -> Path:
    """AMR_CDS.fa, whose pipe-delimited headers become the ``~~~`` convention."""
    raw = _get(f"{NCBI_LATEST}/AMR_CDS.fa", work / "AMR_CDS.fa", progress)
    records = [(converted, seq) for header, seq in read_fasta(raw)
               if (converted := ncbi_header(header))]
    return _write_sequences(records, work)


def ncbi_header(header: str) -> str:
    """``protAcc|nuclAcc|n|n|gene|family|product`` -> the ``~~~`` convention."""
    fields = header.split("|")
    if len(fields) < 6:
        return ""
    accession = fields[1] or fields[0]
    gene = fields[4] or fields[5]
    product = fields[6].split(" ", 1)[-1] if len(fields) > 6 else ""
    return f"ncbi~~~{gene}~~~{accession}~~~{product.replace('_', ' ')}"


# ----------------------------------------------------------------------- CARD
def stage_card(work: Path, progress=None) -> Path:
    """CARD's homolog-model nucleotides; headers are ``gb|ACC|+|range|ARO:n|name``."""
    archive = _get("https://card.mcmaster.ca/latest/data", work / "card.tar.bz2", progress)
    unpacked = _unpack(archive, work / "raw")
    fasta = unpacked / "nucleotide_fasta_protein_homolog_model.fasta"
    if not fasta.exists():
        raise HydraError(f"CARD archive has no nucleotide_fasta_protein_homolog_model.fasta "
                         f"(found: {', '.join(p.name for p in unpacked.iterdir())})")
    records = [(converted, seq) for header, seq in read_fasta(fasta)
               if (converted := card_header(header))]
    return _write_sequences(records, work)


def card_header(header: str) -> str:
    """``gb|ACC|+|start-end|ARO:n|name [organism]`` -> the ``~~~`` convention."""
    fields = header.split("|")
    if len(fields) < 6:
        return ""
    name, _, organism = fields[5].partition(" [")
    return f"card~~~{name.strip()}~~~{fields[1]}~~~{organism.rstrip('] ').strip()}"


# ---------------------------------------------------- ResFinder/PlasmidFinder
def _stage_cge(work: Path, url: str, tag: str, progress=None) -> Path:
    """The CGE databases: one .fsa per class, headers ``gene_n_accession``."""
    archive = _get(url, work / f"{tag}.tar.gz", progress)
    unpacked = _unpack(archive, work / "raw")
    records = [(cge_header(header, tag, fsa.stem), seq)
               for fsa in sorted(unpacked.rglob("*.fsa"))
               for header, seq in read_fasta(fsa)]
    return _write_sequences(records, work)


def cge_header(header: str, tag: str, drug_class: str) -> str:
    """``gene_copy_accession`` -> the ``~~~`` convention.

    Gene names carry internal underscores, so the split runs from the right: the
    last field is the accession and the one before it the copy number.
    """
    token = header.split()[0]
    parts = token.rsplit("_", 2)
    gene, accession = (parts[0], parts[2]) if len(parts) == 3 else (token, "")
    return f"{tag}~~~{gene.strip('_')}~~~{accession}~~~{drug_class}"


def stage_resfinder(work: Path, progress=None) -> Path:
    return _stage_cge(work, "https://bitbucket.org/genomicepidemiology/resfinder_db/"
                            "get/master.tar.gz", "resfinder", progress)


def stage_plasmidfinder(work: Path, progress=None) -> Path:
    return _stage_cge(work, "https://bitbucket.org/genomicepidemiology/plasmidfinder_db/"
                            "get/master.tar.gz", "plasmidfinder", progress)


# ----------------------------------------------------------------------- VFDB
VFDB_HEADER = re.compile(r"^(?P<id>\S+?)\((?:gb\|)?(?P<acc>[^)]+)\)\s+\((?P<gene>[^)]+)\)\s+"
                         r"(?P<product>.*?)(?:\s+\[[^\]]*\])*$")


def stage_vfdb(work: Path, progress=None) -> Path:
    """VFDB core set; headers are ``VFGnnn(gb|ACC) (gene) product [class] [organism]``."""
    archive = _get("http://www.mgc.ac.cn/VFs/Down/VFDB_setA_nt.fas.gz",
                   work / "vfdb.fas.gz", progress)
    plain = work / "vfdb.fas"
    with gzip.open(archive, "rb") as src, open(plain, "wb") as dst:
        shutil.copyfileobj(src, dst)
    records = [(vfdb_header(header), seq) for header, seq in read_fasta(plain)]
    return _write_sequences(records, work)


def vfdb_header(header: str) -> str:
    """``VFGnnn(gb|ACC) (gene) product [class] [organism]`` -> the ``~~~`` convention.

    VFDB reuses the accession as the gene symbol when no symbol is known, so a
    record named after its accession is upstream's own labelling, not a parse
    failure.
    """
    match = VFDB_HEADER.match(header)
    if not match:
        return f"vfdb~~~{header.split()[0]}~~~~~~{header}"
    return f"vfdb~~~{match['gene']}~~~{match['acc']}~~~{match['product'].strip()}"


# -------------------------------------------------------------------- PubMLST
#: PubMLST's own index of every 7-locus scheme, with a direct URL per locus and
#: per profile table. This is what the scheme names below are derived from.
PUBMLST_INDEX = "https://pubmlst.org/data/dbases.xml"

#: Scheme names PubMLST's index would produce, mapped to the spelling already in
#: circulation. Only these few differ; the rest of the index agrees. Without this
#: an existing --scheme argument would stop resolving after a download.
PUBMLST_ALIASES = {
    "escherichia": "ecoli_achtman_4",
    "ecoli_2": "ecoli",
    "salmonella_2": "senterica_achtman_2",
    "diphtheria_3": "cdiphtheriae",
}

#: How many locus files to pull at once. PubMLST throttles: six concurrent
#: requests draws HTTP 429 and silently costs whole schemes, three does not.
PUBMLST_WORKERS = 3


def pubmlst_scheme_name(species: str, profiles_url: str) -> str:
    """Name one scheme from its database slug and its position in that database.

    ``Acinetobacter baumannii#2`` in ``pubmlst_abaumannii_seqdef`` becomes
    ``abaumannii_2``, which is the spelling already in use. Where the derived
    name differs from that spelling, :data:`PUBMLST_ALIASES` maps it over.
    """
    match = re.search(r"/db/(?:pubmlst_)?([a-z0-9_]+?)_seqdef/schemes/(\d+)/", profiles_url)
    if not match:
        return ""
    slug, scheme_id = match.group(1), int(match.group(2))
    ordinal = re.search(r"#(\d+)$", species.strip())
    index = int(ordinal.group(1)) if ordinal else scheme_id
    name = slug if index == 1 else f"{slug}_{index}"
    return PUBMLST_ALIASES.get(name, name)


def pubmlst_organism(species: str) -> tuple[str, str]:
    """Split ``Escherichia coli#1`` into a genus and species for the scheme map."""
    cleaned = re.sub(r"#\d+$", "", species.strip())
    cleaned = re.sub(r"\s+(spp?\.|species complex|complex)$", "", cleaned, flags=re.I)
    parts = cleaned.split()
    if not parts:
        return "", ""
    return parts[0], (parts[1] if len(parts) > 1 else "")


def _pubmlst_schemes(index_xml: Path) -> list[dict]:
    """Every 7-locus scheme in the index, with its locus and profile URLs."""
    import xml.etree.ElementTree as ET  # noqa: PLC0415 - only a download needs it

    try:
        root = ET.parse(index_xml).getroot()
    except ET.ParseError as exc:
        raise HydraError(f"PubMLST index is not valid XML: {exc}") from exc
    out = []
    for entry in root.findall("species"):
        species = (entry.text or "").strip()
        profiles = entry.find("mlst/database/profiles/url")
        if profiles is None or not profiles.text:
            continue
        name = pubmlst_scheme_name(species, profiles.text.strip())
        if not name:
            continue
        loci = {}
        for locus in entry.findall("mlst/database/loci/locus"):
            url = locus.find("url")
            if locus.text and url is not None and url.text:
                loci[locus.text.strip()] = url.text.strip()
        if loci:
            out.append({"scheme": name, "species": species,
                        "profiles": profiles.text.strip(), "loci": loci})
    return out


def stage_pubmlst(work: Path, progress=None) -> Path:
    """Download every PubMLST 7-locus scheme into the layout the importer wants."""
    from concurrent.futures import ThreadPoolExecutor  # noqa: PLC0415

    index = _get(PUBMLST_INDEX, work / "dbases.xml", progress)
    schemes = _pubmlst_schemes(index)
    if not schemes:
        raise HydraError("PubMLST index listed no 7-locus schemes")
    root = work / "pubmlst"
    root.mkdir(parents=True, exist_ok=True)

    jobs = []
    for scheme in schemes:
        directory = root / scheme["scheme"]
        directory.mkdir(parents=True, exist_ok=True)
        jobs.append((scheme["profiles"], directory / f"{scheme['scheme']}.txt"))
        for locus, url in scheme["loci"].items():
            # The importer strips the locus prefix off each allele name, so the
            # file has to be named for the locus exactly as the profile table
            # spells it.
            jobs.append((url, directory / f"{locus}.tfa"))

    LOG.info("downloading %d schemes (%d files) from PubMLST", len(schemes), len(jobs))
    failures: list[str] = []
    done = 0

    def fetch(job):
        url, dest = job
        try:
            _get(url, dest)
        except HydraError as exc:
            return f"{dest.parent.name}/{dest.name}: {exc}"
        return None

    with ThreadPoolExecutor(max_workers=PUBMLST_WORKERS) as pool:
        for problem in pool.map(fetch, jobs):
            done += 1
            if problem:
                failures.append(problem)
            elif progress and done % 100 == 0:
                progress(f"pubmlst: {done}/{len(jobs)} files")
            if done % 200 == 0:
                LOG.info("  %d/%d files", done, len(jobs))

    complete = _prune_incomplete_schemes(root, schemes)
    if not complete:
        raise HydraError("no PubMLST scheme downloaded completely; "
                         f"first problem: {failures[0] if failures else 'unknown'}")
    if failures:
        LOG.warning("%d PubMLST file(s) failed; %d of %d schemes are complete",
                    len(failures), len(complete), len(schemes))
        for problem in failures[:3]:
            LOG.warning("  %s", problem)
        dropped = sorted({item["scheme"] for item in schemes} - complete)
        if dropped:
            LOG.warning("  schemes not installed: %s", ", ".join(dropped))

    _write_scheme_map(work / "scheme_species_map.tab", schemes, complete)
    return work


def _prune_incomplete_schemes(root: Path, schemes: list[dict]) -> set[str]:
    """Drop schemes whose profile table or any locus failed to download.

    A scheme missing one locus would silently type every isolate as an
    incomplete profile, so it is better not to install it at all.
    """
    complete = set()
    for scheme in schemes:
        directory = root / scheme["scheme"]
        profile = directory / f"{scheme['scheme']}.txt"
        loci_present = all((directory / f"{locus}.tfa").exists()
                           and (directory / f"{locus}.tfa").stat().st_size > 0
                           for locus in scheme["loci"])
        if profile.exists() and profile.stat().st_size > 0 and loci_present:
            complete.add(scheme["scheme"])
        elif directory.exists():
            shutil.rmtree(directory, ignore_errors=True)
    return complete


def _write_scheme_map(path: Path, schemes: list[dict], complete: set[str]) -> None:
    """The scheme -> genus/species table, taken from PubMLST's own labelling."""
    lines = ["#SCHEME\tGENUS\tSPECIES"]
    for scheme in sorted(schemes, key=lambda s: s["scheme"]):
        if scheme["scheme"] not in complete:
            continue
        genus, species = pubmlst_organism(scheme["species"])
        lines.append(f"{scheme['scheme']}\t{genus}\t{species}")
    path.write_text("\n".join(lines) + "\n")


# ------------------------------------------------------ lineage / species data
#: The typing loci and the species sketches ship inside one Kleborate release.
KLEBORATE_LATEST = "https://api.github.com/repos/klebgenomics/Kleborate/releases/latest"
KLEBORATE_FALLBACK = "https://github.com/klebgenomics/Kleborate/archive/refs/tags/v3.2.4.tar.gz"

#: One process fetching both databases should download the release once.
_ARCHIVE_CACHE: dict[str, Path] = {}


def _kleborate_release_url() -> str:
    """The newest release tarball, or a known-good tag if the API is unavailable."""
    try:
        request = urllib.request.Request(KLEBORATE_LATEST, headers={"User-Agent": USER_AGENT})
        with urllib.request.urlopen(request, timeout=TIMEOUT) as response:
            tag = json.load(response).get("tag_name")
        if tag:
            return f"https://github.com/klebgenomics/Kleborate/archive/refs/tags/{tag}.tar.gz"
    except (urllib.error.URLError, OSError, ValueError) as exc:
        LOG.debug("could not read the latest Kleborate release (%s); using %s",
                  exc, KLEBORATE_FALLBACK)
    return KLEBORATE_FALLBACK


def _kleborate_modules(work: Path, progress=None) -> Path:
    """Unpack a Kleborate release and return its module directory."""
    url = _kleborate_release_url()
    cached = _ARCHIVE_CACHE.get(url)
    if cached is None or not cached.exists():
        cached = _get(url, work / "kleborate.tar.gz", progress)
        _ARCHIVE_CACHE[url] = cached
    unpacked = _unpack(cached, work / "raw")
    modules = sorted(unpacked.glob("*/kleborate/modules"))
    if not modules:
        raise HydraError(f"{url} does not contain kleborate/modules")
    return modules[0]


def stage_lineage(work: Path, progress=None) -> Path:
    return _kleborate_modules(work, progress)


def stage_species(work: Path, progress=None) -> Path:
    data = _kleborate_modules(work, progress) / "enterobacterales__species" / "data"
    if not any(data.glob("*.msh")):
        raise HydraError(f"no Mash sketch found under {data}")
    return data


#: Databases that can be fetched without human intervention.
FETCHERS: dict[str, Callable[..., Path]] = {
    "protein": stage_protein,
    "ncbi": stage_ncbi,
    "card": stage_card,
    "resfinder": stage_resfinder,
    "plasmidfinder": stage_plasmidfinder,
    "vfdb": stage_vfdb,
    "pubmlst": stage_pubmlst,
    "lineage": stage_lineage,
    "species": stage_species,
}


def can_fetch(name: str) -> bool:
    return name in FETCHERS


def stage(name: str, work: Path, progress=None) -> Path:
    """Download *name* into *work* and return the directory to import from."""
    if name not in FETCHERS:
        raise HydraError(f"'{name}' has no automatic download; "
                         f"run 'hydra db download' for where to get it by hand")
    return FETCHERS[name](work, progress)
