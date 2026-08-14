#!/usr/bin/env python3
"""Read-only, deterministic compiler for an aggregate TAG acyl-moiety ledger.

It reads concrete xPOOL_AC_EM reactants plus curated TAG multisets, builds a
fresh in-memory COBRA model, validates it, and emits canonical JSON only. It
never edits or writes an SBML model.
"""

from __future__ import annotations

import contextlib
import csv
import hashlib
import inspect
import io
import json
import math
import os
import re
import tempfile
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import cobra
from cobra import Metabolite, Model, Reaction
from cobra.flux_analysis import pfba
from cobra.util.solver import linear_reaction_coefficients

POOL_REACTION_ID = "xPOOL_AC_EM"
EXPECTED_CHAIN_COUNT = 7
EXPECTED_TAG_COUNT = 35
EXPECTED_TERMINAL_ROUTES = 81
EXPECTED_REACTION_COUNT = 95
EXPECTED_GENERATED_METABOLITE_COUNT = 44
SCAFFOLD_FORMULA = "C3H8O3"
SCAFFOLD_CHARGE = 0

_SPEC_ROOT_KEYS = frozenset(
    {
        "schema_version",
        "title",
        "purpose",
        "activation_ready",
        "execution_mode",
        "public_compile_contract",
        "inputs",
        "source_chemical_convention",
        "candidate_gpr_curation",
        "required_blockers",
        "chains",
        "curation_selection",
        "architecture",
        "required_manifest",
        "validation_gates",
        "capability_boundary",
    }
)
_SPEC_TITLE = "Read-only aggregate acyl-moiety ledger contract"
_SPEC_PURPOSE = "Plan an in-memory, non-activatable TAG acyl-moiety ledger without changing model.xml."
_SPEC_EXECUTION_MODE = {
    "source_model_immutable": True,
    "build_scope": "in_memory_only",
    "allowed_artifact": "canonical_json_manifest_only",
    "xml_output_forbidden": True,
    "model_activation_forbidden": True,
    "immutable_input_alias_output_forbidden": True,
}
_SPEC_PUBLIC_COMPILE_CONTRACT = {
    "entrypoint": "compile_dry_run",
    "required_parameters": ["model_path", "curation_csv"],
    "optional_parameters": ["spec_path", "verify_determinism"],
    "forbidden_parameters": ["activation_opt_in"],
}
_SPEC_INPUTS = {
    "source_model": "model.xml",
    "curation": "data/lipid_combo_curation.csv",
    "required_sha256": [
        "model.xml",
        "data/lipid_combo_curation.csv",
        "data/lipid_moiety_ledger_spec.json",
    ],
}
_SPEC_CURATION_SELECTION = {
    "layer": "tri",
    "verdict": "keep",
    "composition_key": "unordered member_chains multiset",
    "expected_tag_outputs": EXPECTED_TAG_COUNT,
}
_SPEC_ARCHITECTURE = {
    "first_acylation_reactions": EXPECTED_CHAIN_COUNT,
    "second_acylation_reactions": EXPECTED_CHAIN_COUNT,
    "terminal_tag_routes": EXPECTED_TERMINAL_ROUTES,
    "biochemical_reactions_exact": EXPECTED_REACTION_COUNT,
    "generated_metabolites_exact": EXPECTED_GENERATED_METABOLITE_COUNT,
    "generated_metabolites_maximum": EXPECTED_GENERATED_METABOLITE_COUNT,
}
_SPEC_REQUIRED_MANIFEST = {
    "canonical_serialization": "UTF-8 JSON, sorted keys, compact separators, trailing newline",
    "must_include": [
        "activation_ready",
        "input_sha256",
        "architecture",
        "source_chemical_convention",
        "chain_ledger_ids",
        "tag_outputs",
        "abundance",
        "candidate_gpr",
        "validation",
        "blockers",
    ],
    "validation_required_fields": [
        "mass_balance",
        "occupancy",
        "chain_identity",
        "ids",
        "budget",
        "fba_probe",
        "abundance",
        "source_chemical_convention",
        "source_model_baseline",
        "source_model_unchanged",
    ],
    "validation_flat_fields": [
        "generated_biochemical_reactions_balanced",
        "generated_metabolites_within_budget",
        "chain_identity_conserved",
        "occupancy_conserved",
        "source_coa_acyl_formula_charge_verified",
        "per_tag_required_chain_cut",
        "all_supply_cut",
        "fba_status",
        "pfba_status",
    ],
}
_SPEC_VALIDATION_GATES = {
    "generated_biochemical_reactions_element_and_charge_balanced": True,
    "charge_balance_interpretation": (
        "Generated reactions always use the declared biochemical-pH anionic canonical "
        "tuples. Source formula/charge and canonical-identity verification is reported "
        "separately and opens a blocker on mismatch; the source model is never repaired."
    ),
    "source_coa_acyl_formula_charge_verified": (
        "The same-compartment CoA and all seven xPOOL_AC_EM acyl-CoAs exactly "
        "match their declared formula, integer charge, chemical identity, and "
        "canonical mappings."
    ),
    "chain_identity_conserved": {
        "chains": EXPECTED_CHAIN_COUNT,
        "statement": "For each chain i, w_i^T S = 0 over generated biochemical reactions.",
    },
    "occupancy_conserved": {
        "statement": (
            "For every generated biochemical reaction, sum_i L_i - mono - 2*di = 0."
        )
    },
    "per_tag_required_chain_cut": (
        "Closing any required chain supply makes that TAG demand exactly zero."
    ),
    "all_supply_cut": (
        "Closing all seven chain supplies makes every TAG demand and every "
        "free-ledger demand exactly zero."
    ),
    "all_supply_cut_internal_flux": (
        "With all chain supplies closed, all 95 internal dry-run reaction maxima "
        "are exactly zero."
    ),
    "fba": "optimal",
    "pfba": (
        "optimal with every selected TAG demand flux equal to 1.0 within numerical tolerance"
    ),
    "determinism": "Two independent plans have byte-identical canonical JSON.",
    "cross_process_determinism": (
        "CLI manifests are byte-identical under distinct PYTHONHASHSEED values."
    ),
    "source_model_unchanged": True,
    "source_model_baseline_scope": "independent_source_baseline_not_integrated_regression",
    "activation_ready": False,
    "xml_output_rejected": True,
    "immutable_input_alias_output_rejected": True,
    "atomic_json_output": True,
}
_SPEC_CAPABILITY_BOUNDARY = {
    "retained": (
        "Aggregate histogram of the seven acyl-chain identities selected only from "
        "tri/keep/member_chains."
    ),
    "curation_inputs_retained": [
        "tri",
        "keep",
        "member_chains",
        "prob",
        "cumulative_coverage",
        "xPOOL abundance weights",
    ],
    "curation_inputs_not_retained": [],
    "toy_supply_upper_bound": 3,
    "not_retained": [
        "DAG pairing identity",
        "sn-1/sn-2/sn-3 position",
        "cross-compartment provenance",
    ],
    "unsupported_uses": [
        "site-specific enzymes",
        "selective lipolysis",
        "isotope tracing",
        "quantitative TAG composition inference without an explicit soft-prior probe",
        "quantitative TAG yield inference",
    ],
}
_SPEC_SOURCE_CHEMICAL_CONVENTION = {
    "name": "biochemical_pH_anionic",
    "interpretation": (
        "The same-compartment CoA and all seven xPOOL_AC_EM acyl-CoAs must use "
        "the explicit biochemical-pH anionic formula/charge tuples below. The "
        "compiler only reads and verifies those tuples and canonical identities; "
        "it never repairs the source model."
    ),
    "annotation_policy": {
        "accepted_value_shapes": ["string", "list_of_strings"],
        "required_canonical_fields": ["bigg", "kegg", "mnxm", "chebi"],
        "wrong_same_namespace_forbidden": True,
    },
    "same_compartment_coa": {
        "chemical_identity": "coenzyme A",
        "formula": "C21H32N7O16P3S",
        "charge": -4,
        "canonical_mapping": {
            "bigg": "coa",
            "kegg": "C00010",
            "mnxm": "MNXM12",
            "chebi": "CHEBI:57287",
        },
    },
    "acyl_coas": {
        "lauroyl": {
            "chemical_identity": "lauroyl-CoA",
            "formula": "C33H54N7O17P3S",
            "charge": -4,
            "canonical_mapping": {
                "bigg": "ddcacoa",
                "kegg": "C01832",
                "mnxm": "MNXM363",
                "chebi": "CHEBI:57375",
            },
        },
        "myristoyl": {
            "chemical_identity": "myristoyl-CoA",
            "formula": "C35H58N7O17P3S",
            "charge": -4,
            "canonical_mapping": {
                "bigg": "tdcoa",
                "kegg": "C02593",
                "mnxm": "MNXM224",
                "chebi": "CHEBI:57385",
            },
        },
        "palmitoyl": {
            "chemical_identity": "palmitoyl-CoA",
            "formula": "C37H62N7O17P3S",
            "charge": -4,
            "canonical_mapping": {
                "bigg": "pmtcoa",
                "kegg": "C00154",
                "mnxm": "MNXM88",
                "chebi": "CHEBI:57379",
            },
        },
        "palmitoleoyl": {
            "chemical_identity": "palmitoleoyl-CoA",
            "formula": "C37H60N7O17P3S",
            "charge": -4,
            "canonical_mapping": {
                "bigg": "hdcoa",
                "kegg": "C21072",
                "mnxm": "MNXM781",
                "chebi": "CHEBI:61540",
            },
        },
        "stearoyl": {
            "chemical_identity": "stearoyl-CoA",
            "formula": "C39H66N7O17P3S",
            "charge": -4,
            "canonical_mapping": {
                "bigg": "stcoa",
                "kegg": "C00412",
                "mnxm": "MNXM272",
                "chebi": "CHEBI:57394",
            },
        },
        "oleoyl": {
            "chemical_identity": "oleoyl-CoA",
            "formula": "C39H64N7O17P3S",
            "charge": -4,
            "canonical_mapping": {
                "bigg": "ocdce9coa",
                "kegg": "C00510",
                "mnxm": "MNXM686",
                "chebi": "CHEBI:57387",
            },
        },
        "linoleoyl": {
            "chemical_identity": "linoleoyl-CoA",
            "formula": "C39H62N7O17P3S",
            "charge": -4,
            "canonical_mapping": {
                "bigg": "lnlccoa",
                "kegg": "C02050",
                "mnxm": "MNXM638",
                "chebi": "CHEBI:57383",
            },
        },
    },
}
_SPEC_BLOCKERS = [
    {
        "id": "read_only_manifest_only_prototype",
        "severity": "blocking",
        "required_manifest_status": "open",
        "effect": (
            "activation_ready must remain false because this compiler only emits "
            "a read-only JSON manifest."
        ),
    },
    {
        "id": "candidate_gpr_compartment_evidence_unresolved",
        "severity": "blocking",
        "required_manifest_status": "open",
        "effect": (
            "activation_ready must remain false: candidate GPRs are supported at "
            "the reaction-family level, but xPOOL_AC_EM is in the endosomal-membrane "
            "context and its source templates are disabled."
        ),
    },
    {
        "id": "aggregate_ledger_lacks_molecular_species_provenance",
        "severity": "blocking",
        "required_manifest_status": "open",
        "effect": (
            "activation_ready must remain false because the aggregate ledger does "
            "not retain DAG pairing, sn position, or cross-compartment provenance."
        ),
    },
    {
        "id": "chemistry_source_not_normalized",
        "severity": "blocking",
        "required_manifest_status": "open_when_source_verification_fails",
        "effect": (
            "activation_ready must remain false whenever the immutable source model "
            "does not exactly match the declared biochemical-pH anionic tuples and "
            "canonical identities."
        ),
    },
]

# These rules deliberately describe candidate reaction-family assignments only.
# They remain manifest metadata on the canonical generated model.  A separate
# evaluation copy can materialize them for a clearly labelled memote experiment.
_CANDIDATE_GPRS: dict[str, dict[str, Any]] = {
    "first_acylation": {
        "gpr": "YALI1C00230g",
        "genes": [{
            "id": "YALI1C00230g", "yali0": "YALI0C00209g", "symbol": "SCT1",
            "function": "glycerol-3-phosphate/dihydroxyacetone phosphate acyltransferase",
            "evidence_status": "curated annotation",
        }],
        "status": "candidate_only",
        "source_urls": ["https://pmc.ncbi.nlm.nih.gov/articles/PMC9505966/"],
        "caveat": "No direct seven-acyl-chain substrate spectrum or endosomal-membrane localization evidence.",
    },
    "second_acylation": {
        "gpr": "YALI1E22736g or YALI1F25951g",
        "genes": [
            {
                "id": "YALI1E22736g", "yali0": "YALI0E18964g", "symbol": "SLC1",
                "function": "1-acyl-sn-glycerol-3-phosphate acyltransferase",
                "evidence_status": "curated annotation",
            },
            {
                "id": "YALI1F25951g", "yali0": "YALI0F19514g", "symbol": "ALE1",
                "function": "lysophospholipid acyltransferase with LPAAT/LPCAT annotation",
                "evidence_status": "curated annotation",
            },
        ],
        "status": "candidate_only",
        "source_urls": ["https://pmc.ncbi.nlm.nih.gov/articles/PMC9505966/"],
        "caveat": "The OR rule denotes model alternatives, not experimentally proven chain interchangeability or endosomal-membrane localization.",
    },
    "terminal_tag_route": {
        "gpr": "YALI1E38810g",
        "genes": [{
            "id": "YALI1E38810g", "yali0": "YALI0E32769g", "symbol": "DGA1",
            "function": "DGAT2-family acyl-CoA:diacylglycerol acyltransferase",
            "evidence_status": "experimentally verified",
        }],
        "status": "candidate_only",
        "source_urls": ["https://pmc.ncbi.nlm.nih.gov/articles/PMC3161177/"],
        "caveat": "DGA1 supports lipid-particle/ER TAG synthesis, not the current endosomal-membrane source context.",
    },
}

_FORMULA_TOKEN = re.compile(r"([A-Z][a-z]?)(\d*)")
_SAFE_CHAIN_ID = re.compile(r"^[a-z][a-z0-9_]*$")
_ELEMENTS = frozenset(
    "H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co "
    "Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb "
    "Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re "
    "Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es "
    "Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og".split()
)


class LedgerError(ValueError):
    """A source input violates the closed read-only ledger contract."""


@dataclass(frozen=True)
class Chain:
    """Concrete acyl-CoA source and its conserved acyl-residue ledger identity."""

    id: str
    curation_member: str
    acyl_coa_id: str
    acyl_coa_name: str
    acyl_coa_formula: str
    acyl_coa_charge: int
    coa_id: str
    coa_name: str
    coa_formula: str
    coa_charge: int
    compartment: str
    residue_formula: str
    residue_charge: int


@dataclass(frozen=True)
class LipidClass:
    """A generated backbone, acyl ledger, or definite TAG metabolite."""

    id: str
    label: str
    formula: str
    charge: int
    compartment: str
    generated: bool
    kind: str


@dataclass(frozen=True)
class TransitionSpec:
    """A deterministic biochemical-reaction specification."""

    id: str
    label: str
    kind: str
    stoichiometry: tuple[tuple[str, int], ...]
    required_chains: tuple[str, ...] = ()
    candidate_gpr: str | None = None


@dataclass
class DryRunResult:
    """In-memory model, compiler structures, and canonical JSON manifest."""

    model: Model
    chains: tuple[Chain, ...]
    lipid_classes: tuple[LipidClass, ...]
    transitions: tuple[TransitionSpec, ...]
    manifest: dict[str, Any]
    canonical_json: bytes

    @property
    def canonical_manifest(self) -> bytes:
        return self.canonical_json


def canonical_json_bytes(value: Mapping[str, Any]) -> bytes:
    """Canonical UTF-8 JSON: sorted keys, compact separators, trailing newline."""

    return (
        json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False)
        + "\n"
    ).encode("utf-8")


def sha256_file(path: str | Path) -> str:
    """Return a file SHA-256 or fail closed."""

    source = Path(path)
    if not source.is_file():
        raise LedgerError(f"required input is not a file: {source}")
    digest = hashlib.sha256()
    with source.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_input_fingerprints(
    input_paths: Mapping[str, str | Path], expected_sha256: Mapping[str, str]
) -> None:
    """Require every immutable input to retain its expected fingerprint."""

    if set(input_paths) != set(expected_sha256):
        raise LedgerError("input path and SHA-256 key sets differ")
    for name, path in sorted(input_paths.items()):
        actual = sha256_file(path)
        if actual != expected_sha256[name]:
            raise LedgerError(f"input fingerprint changed: {name}")


def validate_manifest_output_path(
    output_path: str | Path, input_paths: Mapping[str, str | Path]
) -> Path:
    """Reject an output that aliases an immutable input by path, symlink, or inode."""

    output = Path(output_path)
    if output.suffix.lower() != ".json":
        raise LedgerError("manifest output must have a .json suffix")
    if not output.parent.is_dir():
        raise LedgerError(f"manifest output directory does not exist: {output.parent}")
    if output.is_symlink():
        raise LedgerError("refusing to replace an existing symlink output")
    if output.exists() and not output.is_file():
        raise LedgerError("manifest output must be a regular file when it exists")

    resolved_output = output.resolve(strict=False)
    for name, input_path in sorted(input_paths.items()):
        source = Path(input_path)
        if not source.is_file():
            raise LedgerError(f"immutable input is not a file: {name}")
        if resolved_output == source.resolve(strict=True):
            raise LedgerError(
                f"manifest output aliases immutable input by path: {name}"
            )
        if output.exists():
            try:
                if os.path.samefile(output, source):
                    raise LedgerError(
                        f"manifest output aliases immutable input by inode: {name}"
                    )
            except OSError as error:
                raise LedgerError(
                    f"could not compare output with immutable input: {name}"
                ) from error
    return output


def atomic_write_manifest(
    output_path: str | Path,
    payload: bytes,
    input_paths: Mapping[str, str | Path],
    expected_sha256: Mapping[str, str],
) -> Path:
    """Atomically write JSON after a post-write immutable-input fingerprint check."""

    output = validate_manifest_output_path(output_path, input_paths)
    validate_input_fingerprints(input_paths, expected_sha256)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{output.name}.", suffix=".tmp", dir=str(output.parent)
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        validate_input_fingerprints(input_paths, expected_sha256)
        validate_manifest_output_path(output, input_paths)
        os.replace(temporary, output)
    except BaseException:
        try:
            temporary.unlink(missing_ok=True)
        finally:
            raise
    return output


def parse_formula(formula: str) -> dict[str, int]:
    """Parse only explicit, unparenthesized molecular formulas."""

    if not isinstance(formula, str) or not formula.strip():
        raise LedgerError(f"formula is missing or not explicit: {formula!r}")
    text = formula.strip()
    counts: Counter[str] = Counter()
    position = 0
    while position < len(text):
        match = _FORMULA_TOKEN.match(text, position)
        if match is None:
            raise LedgerError(f"unsupported formula syntax: {formula!r}")
        element, count_text = match.groups()
        if element not in _ELEMENTS:
            raise LedgerError(f"unknown or pseudo element {element!r}: {formula!r}")
        count = int(count_text) if count_text else 1
        if count <= 0:
            raise LedgerError(f"formula has non-positive count: {formula!r}")
        counts[element] += count
        position = match.end()
    return dict(counts)


def formula_from_elements(elements: Mapping[str, int]) -> str:
    """Render non-negative integer elements in Hill order."""

    normalized: dict[str, int] = {}
    for element, count in elements.items():
        if element not in _ELEMENTS:
            raise LedgerError(f"invalid element symbol: {element!r}")
        if isinstance(count, bool) or not isinstance(count, int) or count < 0:
            raise LedgerError(f"invalid element count: {element}={count!r}")
        if count:
            normalized[element] = count
    if not normalized:
        raise LedgerError("formula cannot be empty")
    output: list[str] = []
    for element in ("C", "H"):
        if element in normalized:
            count = normalized.pop(element)
            output.append(element + (str(count) if count != 1 else ""))
    for element in sorted(normalized):
        count = normalized[element]
        output.append(element + (str(count) if count != 1 else ""))
    return "".join(output)


def subtract_formula(
    minuend: Mapping[str, int], subtrahend: Mapping[str, int]
) -> dict[str, int]:
    """Subtract element counts, rejecting negative or non-integer residues."""

    result: dict[str, int] = {}
    for element in sorted(set(minuend) | set(subtrahend)):
        left, right = minuend.get(element, 0), subtrahend.get(element, 0)
        if (
            isinstance(left, bool)
            or isinstance(right, bool)
            or not isinstance(left, int)
            or not isinstance(right, int)
        ):
            raise LedgerError(f"non-integer element count for {element}")
        value = left - right
        if value < 0:
            raise LedgerError(
                f"acyl-CoA minus CoA has negative {element}: {left} - {right}"
            )
        if value:
            result[element] = value
    if not result:
        raise LedgerError("acyl-CoA minus CoA gave an empty residue formula")
    return result


def add_formulas(formulas: Iterable[Mapping[str, int]]) -> dict[str, int]:
    """Sum explicit integer molecular formulas."""

    total: Counter[str] = Counter()
    for formula in formulas:
        for element, count in formula.items():
            if isinstance(count, bool) or not isinstance(count, int) or count < 0:
                raise LedgerError(f"invalid formula component {element}={count!r}")
            total[element] += count
    if not total:
        raise LedgerError("cannot sum zero formulas")
    return dict(total)


def _integer_charge(value: Any, label: str) -> int:
    if isinstance(value, bool) or value is None:
        raise LedgerError(f"{label} lacks explicit integer charge")
    if isinstance(value, int):
        return value
    if isinstance(value, float) and math.isfinite(value) and value.is_integer():
        return int(value)
    raise LedgerError(f"{label} has non-integral charge: {value!r}")


def _source_name_base(name: str | None) -> str:
    if not isinstance(name, str) or not name.strip():
        raise LedgerError("source metabolite has no usable name")
    return name.split("_", 1)[0].strip()


def _chain_id_from_member(member: str) -> str:
    text = member.strip()
    if not text.endswith("-CoA"):
        raise LedgerError(f"source/curation member is not acyl-CoA: {member!r}")
    identifier = text[: -len("-CoA")].lower()
    if not _SAFE_CHAIN_ID.fullmatch(identifier):
        raise LedgerError(f"unsupported chain identifier: {identifier!r}")
    return identifier


def _find_same_compartment_coa(model: Model, compartment: str) -> Metabolite:
    """Find the one named CoA source without accepting its chemistry on trust."""

    matches: list[Metabolite] = []
    for metabolite in model.metabolites:
        if metabolite.compartment != compartment:
            continue
        try:
            base = _source_name_base(metabolite.name).lower()
        except LedgerError:
            continue
        if base in {"coenzyme a", "coa"}:
            matches.append(metabolite)
    if len(matches) != 1:
        raise LedgerError(
            f"expected exactly one same-compartment CoA in {compartment!r}; "
            f"found {[item.id for item in matches]!r}"
        )
    return matches[0]


_ANNOTATION_KEY_ALIASES = {
    "bigg": frozenset({"bigg", "bigg.metabolite", "biggm", "biggmetabolite"}),
    "kegg": frozenset({"kegg", "kegg.compound", "keggc", "keggcompound"}),
    "mnxm": frozenset({"mnxm", "metanetx", "metanetx.chemical", "metanetxchemical"}),
    "chebi": frozenset({"chebi", "chebi.compound", "chebicompound"}),
}


def _annotation_namespace(key: object) -> str | None:
    if not isinstance(key, str):
        return None
    normalized = key.strip().lower().replace("_", ".")
    for namespace, aliases in _ANNOTATION_KEY_ALIASES.items():
        if normalized in aliases:
            return namespace
    return None


def _annotation_string_values(value: object, label: str) -> list[str]:
    """Accept the repository's scalar/list annotations, never an implicit shape."""

    if isinstance(value, str):
        values = [value]
    elif isinstance(value, list):
        values = value
    else:
        raise LedgerError(
            f"{label} must be a string or list of strings, got {type(value).__name__}"
        )
    if not values or any(not isinstance(item, str) or not item.strip() for item in values):
        raise LedgerError(f"{label} has an empty or non-string annotation value")
    return [item.strip() for item in values]


def _normalize_annotation_value(namespace: str, value: str) -> str:
    """Normalize identifier URLs and repository prefixes to canonical comparison text."""

    text = value.strip()
    if "/" in text:
        text = text.rsplit("/", 1)[-1]
    prefix_patterns = {
        "bigg": ("bigg.metabolite:", "bigg:"),
        "kegg": ("kegg.compound:", "kegg:"),
        "mnxm": ("metanetx.chemical:", "metanetx:", "mnxm:"),
        "chebi": ("chebi:",),
    }
    lowered = text.lower()
    for prefix in prefix_patterns[namespace]:
        if lowered.startswith(prefix):
            text = text[len(prefix) :]
            break
    if namespace == "bigg":
        if text.startswith("M_"):
            text = text[2:]
        return text.lower()
    if namespace == "kegg":
        if text.startswith("M_"):
            text = text[2:]
        return text.upper()
    if namespace == "mnxm":
        return text.upper()
    if namespace == "chebi":
        text = text.upper()
        return text if text.startswith("CHEBI:") else f"CHEBI:{text}"
    raise LedgerError(f"unsupported canonical annotation namespace: {namespace}")


def _canonical_annotation_values(metabolite: Metabolite, namespace: str) -> set[str]:
    annotation = getattr(metabolite, "annotation", None)
    if not isinstance(annotation, Mapping):
        raise LedgerError(f"{metabolite.id} has no usable annotation mapping")
    values: set[str] = set()
    found_key = False
    for key, raw_value in annotation.items():
        if _annotation_namespace(key) != namespace:
            continue
        found_key = True
        for value in _annotation_string_values(
            raw_value, f"{metabolite.id} annotation {key!r}"
        ):
            values.add(_normalize_annotation_value(namespace, value))
    if not found_key:
        raise LedgerError(f"{metabolite.id} lacks {namespace} canonical annotation")
    return values


def _source_pool_bindings(
    model: Model,
    source_chemical_convention: Mapping[str, Any],
    pool_reaction_id: str,
) -> tuple[tuple[str, str, Metabolite, Metabolite], ...]:
    """Bind source IDs by named pool membership before chemistry is inspected."""

    expected_acyls = _required_mapping(
        source_chemical_convention.get("acyl_coas"),
        "source_chemical_convention.acyl_coas",
    )
    try:
        pool = model.reactions.get_by_id(pool_reaction_id)
    except KeyError as error:
        raise LedgerError(f"pool reaction missing: {pool_reaction_id}") from error

    bindings: list[tuple[str, str, Metabolite, Metabolite]] = []
    seen_ids: set[str] = set()
    seen_chains: set[str] = set()
    for metabolite, coefficient in pool.metabolites.items():
        if coefficient >= 0:
            continue
        member = _source_name_base(metabolite.name)
        chain_id = _chain_id_from_member(member)
        if chain_id not in expected_acyls:
            raise LedgerError(f"unexpected xPOOL acyl-CoA chain: {chain_id}")
        if metabolite.id in seen_ids or chain_id in seen_chains:
            raise LedgerError(
                f"duplicate concrete acyl chain: {metabolite.id}/{chain_id}"
            )
        seen_ids.add(metabolite.id)
        seen_chains.add(chain_id)
        bindings.append((chain_id, member, metabolite, metabolite))

    if len(bindings) != EXPECTED_CHAIN_COUNT or set(expected_acyls) != seen_chains:
        raise LedgerError(
            f"{pool_reaction_id} must expose exactly the seven declared concrete acyl "
            f"reactants; found {sorted(seen_chains)!r}"
        )
    compartments = {metabolite.compartment for _, _, metabolite, _ in bindings}
    if len(compartments) != 1:
        raise LedgerError("global ledger cannot merge cross-compartment acyl inputs")
    compartment = next(iter(compartments))
    coa = _find_same_compartment_coa(model, compartment)
    return tuple(
        sorted(
            (
                (chain_id, member, metabolite, coa)
                for chain_id, member, metabolite, _ in bindings
            ),
            key=lambda item: (item[0], item[2].id),
        )
    )


def _expected_tuple_report(
    metabolite: Metabolite, expected: Mapping[str, Any], label: str
) -> dict[str, Any]:
    """Audit one immutable source species without changing it or stopping early."""

    expected_formula = str(expected["formula"])
    expected_charge = _integer_charge(expected["charge"], f"{label} expected charge")
    expected_identity = str(expected["chemical_identity"])
    expected_mapping = _required_mapping(expected["canonical_mapping"], f"{label}.canonical_mapping")
    errors: list[str] = []
    try:
        actual_formula = formula_from_elements(parse_formula(metabolite.formula))
    except LedgerError as error:
        actual_formula = None
        errors.append(f"formula: {error}")
    try:
        actual_charge = _integer_charge(metabolite.charge, f"{label} {metabolite.id}")
    except LedgerError as error:
        actual_charge = None
        errors.append(f"charge: {error}")
    try:
        actual_identity = _source_name_base(metabolite.name)
    except LedgerError as error:
        actual_identity = None
        errors.append(f"chemical_identity: {error}")
    if actual_formula is not None and actual_formula != expected_formula:
        errors.append(f"formula: expected {expected_formula}, found {actual_formula}")
    if actual_charge is not None and actual_charge != expected_charge:
        errors.append(f"charge: expected {expected_charge}, found {actual_charge}")
    if actual_identity is not None and actual_identity != expected_identity:
        errors.append(
            f"chemical_identity: expected {expected_identity!r}, found {actual_identity!r}"
        )

    annotation_report: dict[str, Any] = {}
    for namespace, expected_value in sorted(expected_mapping.items()):
        try:
            actual_values = _canonical_annotation_values(metabolite, str(namespace))
            expected_normalized = _normalize_annotation_value(
                str(namespace), str(expected_value)
            )
            passed = actual_values == {expected_normalized}
            if not passed:
                errors.append(
                    f"{namespace}: expected only {expected_normalized!r}, found "
                    f"{sorted(actual_values)!r}"
                )
            annotation_report[str(namespace)] = {
                "expected": expected_normalized,
                "actual": sorted(actual_values),
                "passed": passed,
            }
        except LedgerError as error:
            errors.append(f"{namespace}: {error}")
            annotation_report[str(namespace)] = {
                "expected": str(expected_value),
                "actual": [],
                "passed": False,
            }
    return {
        "source_model_id": metabolite.id,
        "expected": {
            "chemical_identity": expected_identity,
            "formula": expected_formula,
            "charge": expected_charge,
            "canonical_mapping": dict(sorted(expected_mapping.items())),
        },
        "actual": {
            "chemical_identity": actual_identity,
            "formula": actual_formula,
            "charge": actual_charge,
            "canonical_mapping": annotation_report,
        },
        "passed": not errors,
        "errors": errors,
    }


def _source_chemical_verification(
    bindings: Sequence[tuple[str, str, Metabolite, Metabolite]],
    source_chemical_convention: Mapping[str, Any],
) -> dict[str, Any]:
    """Return a complete auditable source chemistry report, even when it fails."""

    coa_expected = _required_mapping(
        source_chemical_convention.get("same_compartment_coa"),
        "source_chemical_convention.same_compartment_coa",
    )
    acyl_expected = _required_mapping(
        source_chemical_convention.get("acyl_coas"),
        "source_chemical_convention.acyl_coas",
    )
    coa = bindings[0][3]
    coa_report = _expected_tuple_report(coa, coa_expected, "same-compartment CoA")
    acyl_reports = {
        chain_id: _expected_tuple_report(
            metabolite,
            _required_mapping(acyl_expected[chain_id], f"acyl_coas.{chain_id}"),
            f"acyl-CoA {chain_id}",
        )
        for chain_id, _member, metabolite, _coa in bindings
    }
    errors = [
        f"same_compartment_coa: {error}"
        for error in coa_report["errors"]
    ] + [
        f"{chain_id}: {error}"
        for chain_id, report in sorted(acyl_reports.items())
        for error in report["errors"]
    ]
    passed = coa_report["passed"] and all(
        report["passed"] for report in acyl_reports.values()
    )
    return {
        "convention": source_chemical_convention["name"],
        "status": "verified" if passed else "not_normalized",
        "tuple_origin": "source_model" if passed else "explicit_canonical_curation_tuple",
        "passed": passed,
        "same_compartment_coa": coa_report,
        "acyl_coas": dict(sorted(acyl_reports.items())),
        "errors": errors,
    }


def _chains_from_bindings(
    bindings: Sequence[tuple[str, str, Metabolite, Metabolite]],
    source_chemical_convention: Mapping[str, Any],
) -> tuple[Chain, ...]:
    """Build ledger chains from declared tuples after binding immutable source IDs."""

    coa_expected = _required_mapping(
        source_chemical_convention["same_compartment_coa"],
        "source_chemical_convention.same_compartment_coa",
    )
    acyl_expected = _required_mapping(
        source_chemical_convention["acyl_coas"],
        "source_chemical_convention.acyl_coas",
    )
    coa_formula = str(coa_expected["formula"])
    coa_charge = _integer_charge(coa_expected["charge"], "declared CoA charge")
    parsed_coa_formula = parse_formula(coa_formula)
    loaded: list[Chain] = []
    for chain_id, member, metabolite, coa in bindings:
        expected = _required_mapping(acyl_expected[chain_id], f"acyl_coas.{chain_id}")
        acyl_formula = str(expected["formula"])
        acyl_charge = _integer_charge(expected["charge"], f"declared {chain_id} charge")
        residue = subtract_formula(parse_formula(acyl_formula), parsed_coa_formula)
        if set(residue) != {"C", "H", "O"} or residue.get("O") != 1:
            raise LedgerError(
                f"declared {chain_id} acyl-CoA minus CoA is not a one-oxygen C/H/O "
                "acyl residue"
            )
        if acyl_charge - coa_charge != 0:
            raise LedgerError(f"declared {chain_id} residue charge is not zero")
        loaded.append(
            Chain(
                id=chain_id,
                curation_member=member,
                acyl_coa_id=metabolite.id,
                acyl_coa_name=metabolite.name,
                acyl_coa_formula=acyl_formula,
                acyl_coa_charge=acyl_charge,
                coa_id=coa.id,
                coa_name=coa.name,
                coa_formula=coa_formula,
                coa_charge=coa_charge,
                compartment=metabolite.compartment,
                residue_formula=formula_from_elements(residue),
                residue_charge=0,
            )
        )
    return tuple(sorted(loaded, key=lambda item: (item.id, item.acyl_coa_id)))


def load_chains(
    model: Model,
    source_chemical_convention: Mapping[str, Any] | str | None = None,
    pool_reaction_id: str = POOL_REACTION_ID,
) -> tuple[Chain, ...]:
    """Strictly verify all source tuples and identities, otherwise fail closed."""

    # Retain the former ``load_chains(model, pool_reaction_id)`` call shape,
    # while making its chemistry checks strict under the shipped v2 contract.
    if isinstance(source_chemical_convention, str):
        if pool_reaction_id != POOL_REACTION_ID:
            raise LedgerError("pool reaction ID was provided twice")
        pool_reaction_id = source_chemical_convention
        source_chemical_convention = None
    if source_chemical_convention is None:
        source_chemical_convention = _validated_source_chemical_convention(
            load_spec(_resolve_spec_path(None))
        )
    bindings = _source_pool_bindings(
        model, source_chemical_convention, pool_reaction_id
    )
    verification = _source_chemical_verification(bindings, source_chemical_convention)
    if not verification["passed"]:
        raise LedgerError(
            "source CoA/acyl-CoA chemical convention verification failed: "
            + "; ".join(verification["errors"])
        )
    return _chains_from_bindings(bindings, source_chemical_convention)


def canonicalize_tag_multiset(
    members: Sequence[str], chain_by_member: Mapping[str, Chain]
) -> tuple[str, str, str]:
    """Translate an unordered curation row into sorted internal chain IDs."""

    if len(members) != 3:
        raise LedgerError(f"TAG must contain exactly 3 chains: {members!r}")
    resolved: list[str] = []
    for member in members:
        source_key = member.strip().lower()
        if source_key not in chain_by_member:
            raise LedgerError(f"curation chain absent from xPOOL_AC_EM: {member!r}")
        resolved.append(chain_by_member[source_key].id)
    return tuple(sorted(resolved))  # type: ignore[return-value]


def load_curated_tag_records(
    curation_csv: str | Path, chains: Sequence[Chain]
) -> tuple[dict[str, Any], ...]:
    """Read 35 kept TAG multisets with their probability and coverage metadata."""

    path = Path(curation_csv)
    if not path.is_file():
        raise LedgerError(f"curation CSV is not a file: {path}")
    chain_by_member = {item.curation_member.lower(): item for item in chains}
    selected: list[dict[str, Any]] = []
    tri_rows: list[tuple[float, float]] = []
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {
            "layer", "n_chains", "member_chains", "verdict", "prob",
            "cumulative_coverage",
        }
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise LedgerError(f"curation lacks required columns: {sorted(required)!r}")
        for row in reader:
            if row.get("layer") != "tri":
                continue
            try:
                row_probability = float(str(row["prob"]))
                row_cumulative = float(str(row["cumulative_coverage"]))
            except (TypeError, ValueError) as error:
                raise LedgerError(f"tri row has invalid probability metadata: {row!r}") from error
            if not (math.isfinite(row_probability) and row_probability > 0.0):
                raise LedgerError(f"tri row has non-positive probability: {row!r}")
            if not (math.isfinite(row_cumulative) and 0.0 < row_cumulative <= 1.0):
                raise LedgerError(f"tri row has invalid cumulative coverage: {row!r}")
            tri_rows.append((row_probability, row_cumulative))
            if row.get("verdict") != "keep":
                continue
            if row.get("n_chains") != "3":
                raise LedgerError(f"kept tri row lacks n_chains=3: {row!r}")
            selected.append({
                "composition": canonicalize_tag_multiset(
                    (row.get("member_chains") or "").split(";"), chain_by_member
                ),
                "prob": row_probability,
                "cumulative_coverage": row_cumulative,
                "curation_rank": len(tri_rows),
            })
    if len(selected) != EXPECTED_TAG_COUNT:
        raise LedgerError(
            f"expected {EXPECTED_TAG_COUNT} kept TAG rows, got {len(selected)}"
        )
    compositions = [item["composition"] for item in selected]
    if len(set(compositions)) != len(compositions):
        raise LedgerError("duplicate unordered kept TAG multiset")
    cumulative = [float(item["cumulative_coverage"]) for item in selected]
    if any(right < left for left, right in zip(cumulative, cumulative[1:])):
        raise LedgerError("kept tri cumulative coverage is not nondecreasing")
    total = sum(float(item["prob"]) for item in selected)
    if not (math.isfinite(total) and 0.0 < total <= 1.0 + 1e-8):
        raise LedgerError("kept tri probabilities do not form a valid coverage")
    running = 0.0
    for probability, declared in tri_rows:
        running += probability
        if abs(running - declared) > 1.1e-4:
            raise LedgerError("tri cumulative coverage disagrees with probability ranking")
    if any(right > left + 1e-12 for (left, _), (right, _) in zip(tri_rows, tri_rows[1:])):
        raise LedgerError("tri probabilities are not ranked from high to low")
    return tuple(sorted(selected, key=lambda item: item["composition"]))


def load_curated_tag_multisets(
    curation_csv: str | Path, chains: Sequence[Chain]
) -> tuple[tuple[str, str, str], ...]:
    """Compatibility wrapper returning only canonical TAG multisets."""

    return tuple(item["composition"] for item in load_curated_tag_records(curation_csv, chains))


def _source_pool_abundance_weights(source: Model, chains: Sequence[Chain]) -> dict[str, float]:
    """Normalize the seven concrete acyl-CoA coefficients in xPOOL_AC_EM."""

    pool = source.reactions.get_by_id(POOL_REACTION_ID)
    raw: dict[str, float] = {}
    for chain in chains:
        metabolite = source.metabolites.get_by_id(chain.acyl_coa_id)
        coefficient = float(pool.metabolites.get(metabolite, 0.0))
        if coefficient >= 0.0:
            raise LedgerError(f"{POOL_REACTION_ID} lacks a reactant coefficient for {chain.id}")
        raw[chain.id] = -coefficient
    total = sum(raw.values())
    if not (math.isfinite(total) and total > 0.0):
        raise LedgerError(f"{POOL_REACTION_ID} abundance weights are invalid")
    return {chain_id: raw[chain_id] / total for chain_id in sorted(raw)}


def _abundance_report(
    tags: Sequence[Mapping[str, Any]], source_weights: Mapping[str, float]
) -> dict[str, Any]:
    """Report the retained mixture without making it a default flux constraint."""

    for tag in tags:
        counts = Counter(str(chain) for chain in tag["chains"])
        permutations = math.factorial(3) / math.prod(
            math.factorial(count) for count in counts.values()
        )
        expected_probability = permutations * math.prod(
            float(source_weights[chain_id]) ** count
            for chain_id, count in counts.items()
        )
        if abs(float(tag["prob"]) - expected_probability) > 1e-6:
            raise LedgerError(
                f"curated TAG probability disagrees with xPOOL multinomial: {tag['id']}"
            )
    coverage = sum(float(tag["prob"]) for tag in tags)
    normalized = sum(float(tag["normalized_prob"]) for tag in tags)
    marginal = {
        chain_id: sum(
            list(tag["chains"]).count(chain_id) * float(tag["normalized_prob"])
            for tag in tags
        ) / 3.0
        for chain_id in sorted(source_weights)
    }
    deltas = {
        chain_id: marginal[chain_id] - float(source_weights[chain_id])
        for chain_id in sorted(source_weights)
    }
    if not _near(normalized, 1.0, 1e-7):
        raise LedgerError(f"normalized kept TAG probabilities do not sum to one: {normalized}")
    return {
        "source_xpool_weights": dict(sorted(source_weights.items())),
        "kept_probability_coverage": coverage,
        "dropped_probability": 1.0 - coverage,
        "normalized_probability_sum": normalized,
        "normalized_probability_semantics": "conditional_on_retained_35_tag_set",
        "probability_model": "unordered_multiset_sum_over_ordered_chain_draws",
        "implied_chain_marginals": marginal,
        "marginal_minus_source_weight": deltas,
        "solver_use": "opt_in_soft_prior_only",
        "cumulative_coverage_use": "curation_qc_only",
    }


def _required_mapping(value: Any, label: str) -> Mapping[str, Any]:
    if type(value) is not dict:
        raise LedgerError(f"spec field {label} must be an object")
    return value


def _required_list(value: Any, label: str) -> list[Any]:
    if type(value) is not list:
        raise LedgerError(f"spec field {label} must be a list")
    return value


def _strict_spec_equal(actual: Any, expected: Any) -> bool:
    """Compare JSON-like schema values without Python numeric type coercion."""

    if type(actual) is not type(expected):
        return False
    if type(expected) is dict:
        return (
            set(actual) == set(expected)
            and all(_strict_spec_equal(actual[key], expected[key]) for key in expected)
        )
    if type(expected) is list:
        return len(actual) == len(expected) and all(
            _strict_spec_equal(left, right) for left, right in zip(actual, expected)
        )
    return actual == expected


def _exact_spec_mapping(
    value: Any, expected: Mapping[str, Any], label: str
) -> dict[str, Any]:
    if type(value) is not dict:
        raise LedgerError(f"spec field {label} must be an object")
    actual = dict(value)
    if not _strict_spec_equal(actual, dict(expected)):
        raise LedgerError(f"spec {label} drifted from its exact contract")
    return actual


def _nonnegative_spec_int(value: Any, label: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise LedgerError(f"spec {label} must be a non-negative integer")
    return value


def _validated_spec_chains(spec: Mapping[str, Any]) -> dict[str, Mapping[str, Any]]:
    raw_chains = _required_list(spec.get("chains"), "chains")
    if len(raw_chains) != EXPECTED_CHAIN_COUNT:
        raise LedgerError("spec must contain exactly seven chains")
    expected_keys = {
        "id",
        "curation_member",
        "acyl_residue_formula",
        "carbon_atoms",
        "double_bonds",
    }
    result: dict[str, Mapping[str, Any]] = {}
    members: set[str] = set()
    for index, raw_item in enumerate(raw_chains):
        item = _required_mapping(raw_item, f"chains[{index}]")
        if set(item) != expected_keys:
            raise LedgerError(f"spec chains[{index}] has unknown or missing fields")
        identifier = item["id"]
        member = item["curation_member"]
        formula = item["acyl_residue_formula"]
        carbon = _nonnegative_spec_int(
            item["carbon_atoms"], f"chains[{index}].carbon_atoms"
        )
        double_bonds = _nonnegative_spec_int(
            item["double_bonds"], f"chains[{index}].double_bonds"
        )
        if (
            not isinstance(identifier, str)
            or not _SAFE_CHAIN_ID.fullmatch(identifier)
            or not isinstance(member, str)
            or member != f"{identifier}-CoA"
            or not isinstance(formula, str)
        ):
            raise LedgerError(f"spec chains[{index}] has invalid id/member/formula")
        expected_hydrogen = 2 * carbon - 2 * double_bonds - 2
        expected_formula = {"C": carbon, "H": expected_hydrogen, "O": 1}
        if expected_hydrogen <= 0 or parse_formula(formula) != expected_formula:
            raise LedgerError(
                f"spec chains[{index}] formula is inconsistent with carbon/double_bonds"
            )
        if formula != formula_from_elements(expected_formula):
            raise LedgerError(f"spec chains[{index}] formula is not canonical")
        if identifier in result or member in members:
            raise LedgerError("spec chain IDs and curation members must be unique")
        result[identifier] = item
        members.add(member)
    if "linoleoyl" not in result:
        raise LedgerError("spec must identify all seven declared source acyl-CoAs")
    return result


def _validated_source_chemical_convention(
    spec: Mapping[str, Any],
) -> Mapping[str, Any]:
    """Require the full biochemical-pH tuple and identity contract verbatim."""

    convention = _exact_spec_mapping(
        spec.get("source_chemical_convention"),
        _SPEC_SOURCE_CHEMICAL_CONVENTION,
        "source_chemical_convention",
    )
    coa = _required_mapping(
        convention["same_compartment_coa"],
        "source_chemical_convention.same_compartment_coa",
    )
    parse_formula(str(coa["formula"]))
    _integer_charge(coa["charge"], "declared same-compartment CoA")
    acyls = _required_mapping(
        convention["acyl_coas"], "source_chemical_convention.acyl_coas"
    )
    if set(acyls) != {
        "lauroyl",
        "myristoyl",
        "palmitoyl",
        "palmitoleoyl",
        "stearoyl",
        "oleoyl",
        "linoleoyl",
    }:
        raise LedgerError("source chemical convention must declare exactly seven acyl-CoAs")
    for chain_id, entry in acyls.items():
        entry = _required_mapping(entry, f"source_chemical_convention.acyl_coas.{chain_id}")
        parse_formula(str(entry["formula"]))
        _integer_charge(entry["charge"], f"declared {chain_id} charge")
    return convention


def _validated_spec_chain_source_contract(
    spec: Mapping[str, Any],
) -> dict[str, Mapping[str, Any]]:
    """Bind every declared ledger residue directly to its declared acyl-CoA."""

    chains = _validated_spec_chains(spec)
    convention = _validated_source_chemical_convention(spec)
    acyl_coas = _required_mapping(
        convention["acyl_coas"], "source_chemical_convention.acyl_coas"
    )
    if set(chains) != set(acyl_coas):
        raise LedgerError(
            "spec chain IDs must exactly equal source_chemical_convention.acyl_coas"
        )

    coa = _required_mapping(
        convention["same_compartment_coa"],
        "source_chemical_convention.same_compartment_coa",
    )
    coa_formula = parse_formula(str(coa["formula"]))
    coa_charge = _integer_charge(coa["charge"], "declared same-compartment CoA")
    for chain_id, chain in chains.items():
        acyl_coa = _required_mapping(
            acyl_coas[chain_id], f"source_chemical_convention.acyl_coas.{chain_id}"
        )
        residue = subtract_formula(
            parse_formula(str(acyl_coa["formula"])), coa_formula
        )
        if set(residue) != {"C", "H", "O"} or residue.get("O") != 1:
            raise LedgerError(
                f"declared {chain_id} acyl-CoA minus CoA is not a C/H/O acyl residue"
            )
        if _integer_charge(
            acyl_coa["charge"], f"declared {chain_id} charge"
        ) != coa_charge:
            raise LedgerError(f"declared {chain_id} residue charge is not zero")
        expected_residue_formula = formula_from_elements(residue)
        if chain["acyl_residue_formula"] != expected_residue_formula:
            raise LedgerError(
                f"spec {chain_id} acyl residue does not equal declared acyl-CoA minus CoA"
            )
    return chains


def _validate_public_compile_contract(spec: Mapping[str, Any]) -> None:
    """Check both the declared and the actual public compiler call contract."""

    _exact_spec_mapping(
        spec.get("public_compile_contract"),
        _SPEC_PUBLIC_COMPILE_CONTRACT,
        "public_compile_contract",
    )
    function = globals().get("compile_dry_run")
    if function is None:
        raise LedgerError(
            "compile_dry_run is unavailable for public-contract validation"
        )
    parameters = list(inspect.signature(function).parameters.values())
    expected_names = ["model_path", "curation_csv", "spec_path", "verify_determinism"]
    if [parameter.name for parameter in parameters] != expected_names:
        raise LedgerError(
            "compile_dry_run parameter names drifted from the specification"
        )
    for parameter in parameters[:2]:
        if (
            parameter.kind is not inspect.Parameter.POSITIONAL_OR_KEYWORD
            or parameter.default is not inspect.Parameter.empty
        ):
            raise LedgerError(
                "compile_dry_run required parameters drifted from contract"
            )
    optional = parameters[2:]
    if (
        any(
            parameter.kind is not inspect.Parameter.KEYWORD_ONLY
            for parameter in optional
        )
        or optional[0].default is not None
        or optional[1].default is not True
        or "activation_opt_in" in inspect.signature(function).parameters
    ):
        raise LedgerError("compile_dry_run optional or forbidden parameters drifted")


def load_spec(spec_path: str | Path) -> dict[str, Any]:
    """Load every normative field of the ledger specification exactly."""

    path = Path(spec_path)
    if not path.is_file():
        raise LedgerError(f"specification is not a file: {path}")
    try:
        with path.open(encoding="utf-8") as handle:
            spec = json.load(handle)
    except (OSError, json.JSONDecodeError) as error:
        raise LedgerError(f"could not parse ledger specification: {path}") from error
    spec = dict(_required_mapping(spec, "root"))
    if set(spec) != _SPEC_ROOT_KEYS:
        raise LedgerError("spec root has unknown or missing normative fields")
    if (
        type(spec["schema_version"]) is not int
        or spec["schema_version"] != 2
        or spec["title"] != _SPEC_TITLE
        or spec["purpose"] != _SPEC_PURPOSE
        or spec["activation_ready"] is not False
    ):
        raise LedgerError("spec schema, title, purpose, or activation state drifted")

    _exact_spec_mapping(spec["execution_mode"], _SPEC_EXECUTION_MODE, "execution_mode")
    _exact_spec_mapping(
        spec["public_compile_contract"],
        _SPEC_PUBLIC_COMPILE_CONTRACT,
        "public_compile_contract",
    )
    _exact_spec_mapping(spec["inputs"], _SPEC_INPUTS, "inputs")
    _exact_spec_mapping(
        spec["candidate_gpr_curation"], _CANDIDATE_GPRS, "candidate_gpr_curation"
    )
    _exact_spec_mapping(
        spec["curation_selection"], _SPEC_CURATION_SELECTION, "curation_selection"
    )
    _exact_spec_mapping(spec["architecture"], _SPEC_ARCHITECTURE, "architecture")
    _exact_spec_mapping(
        spec["required_manifest"], _SPEC_REQUIRED_MANIFEST, "required_manifest"
    )
    _exact_spec_mapping(
        spec["validation_gates"], _SPEC_VALIDATION_GATES, "validation_gates"
    )
    _exact_spec_mapping(
        spec["capability_boundary"], _SPEC_CAPABILITY_BOUNDARY, "capability_boundary"
    )

    blockers = _required_list(spec["required_blockers"], "required_blockers")
    if not _strict_spec_equal(blockers, _SPEC_BLOCKERS):
        raise LedgerError("spec required_blockers drifted from current contract")
    _validated_spec_chain_source_contract(spec)
    return spec


def _validate_spec_against_plan(
    spec: Mapping[str, Any],
    chains: Sequence[Chain],
    tag_multisets: Sequence[tuple[str, str, str]],
) -> None:
    """Cross-check declared ledger tuples against the immutable build contract."""

    expected = _validated_spec_chain_source_contract(spec)
    if set(expected) != {chain.id for chain in chains}:
        raise LedgerError("source xPOOL chain IDs differ from the specification")
    source_convention = _validated_source_chemical_convention(spec)
    source_acyls = _required_mapping(
        source_convention["acyl_coas"], "source_chemical_convention.acyl_coas"
    )
    coa = _required_mapping(
        source_convention["same_compartment_coa"],
        "source_chemical_convention.same_compartment_coa",
    )
    for chain in chains:
        item = expected[chain.id]
        source_acyl = _required_mapping(source_acyls[chain.id], f"acyl_coas.{chain.id}")
        residue = parse_formula(chain.residue_formula)
        if (
            item["curation_member"] != chain.curation_member
            or item["acyl_residue_formula"] != chain.residue_formula
            or residue != parse_formula(str(item["acyl_residue_formula"]))
            or chain.coa_formula != coa["formula"]
            or chain.coa_charge != coa["charge"]
            or chain.acyl_coa_formula != source_acyl["formula"]
            or chain.acyl_coa_charge != source_acyl["charge"]
            or chain.residue_charge != 0
        ):
            raise LedgerError(
                f"declared canonical chain differs from spec for {chain.id}"
            )
    if (
        sum(len(set(composition)) for composition in tag_multisets)
        != EXPECTED_TERMINAL_ROUTES
    ):
        raise LedgerError("curated TAG route budget differs from specification")


def _met(
    identifier: str, name: str, formula: str, charge: int, compartment: str
) -> Metabolite:
    parse_formula(formula)
    _integer_charge(charge, identifier)
    return Metabolite(
        identifier, name=name, formula=formula, charge=charge, compartment=compartment
    )


def _add_reaction(
    model: Model,
    transitions: list[TransitionSpec],
    identifier: str,
    label: str,
    kind: str,
    stoichiometry: Mapping[Metabolite, int],
    required_chains: Sequence[str] = (),
) -> None:
    if model.reactions.has_id(identifier):
        raise LedgerError(f"duplicate generated reaction ID: {identifier}")
    normalized: dict[Metabolite, int] = {}
    for metabolite, coefficient in stoichiometry.items():
        if (
            isinstance(coefficient, bool)
            or not isinstance(coefficient, int)
            or not coefficient
        ):
            raise LedgerError(f"invalid stoichiometry in {identifier}")
        normalized[metabolite] = normalized.get(metabolite, 0) + coefficient
    normalized = {
        item: coefficient for item, coefficient in normalized.items() if coefficient
    }
    reaction = Reaction(identifier, name=label)
    reaction.bounds = (0.0, 1000.0)
    reaction.annotation["lipid_moiety_ledger"] = "read_only_dry_run"
    reaction.annotation["ledger_kind"] = kind
    candidate_gpr = _CANDIDATE_GPRS.get(kind)
    if candidate_gpr is not None:
        reaction.annotation["candidate_gpr_status"] = str(candidate_gpr["status"])
    reaction.add_metabolites(normalized)
    model.add_reactions([reaction])
    transitions.append(
        TransitionSpec(
            id=identifier,
            label=label,
            kind=kind,
            stoichiometry=tuple(
                sorted(
                    (
                        (item.id, coefficient)
                        for item, coefficient in normalized.items()
                    ),
                    key=lambda pair: pair[0],
                )
            ),
            required_chains=tuple(sorted(required_chains)),
            candidate_gpr=(None if candidate_gpr is None else str(candidate_gpr["gpr"])),
        )
    )


def _tag_id(composition: Sequence[str]) -> str:
    return "LM_TAG__" + "__".join(composition)


def _route_id(composition: Sequence[str], last_chain: str) -> str:
    return "LM_TAG_ROUTE__" + "__".join(composition) + "__LAST__" + last_chain


def _build_ledger(
    chains: Sequence[Chain], curated_tags: Sequence[Mapping[str, Any]]
) -> tuple[
    Model,
    tuple[LipidClass, ...],
    tuple[TransitionSpec, ...],
    dict[str, str],
    list[dict[str, Any]],
    dict[str, Any],
]:
    tag_multisets = tuple(
        tuple(str(chain) for chain in item["composition"])
        for item in curated_tags
    )
    if len(chains) != 7 or len(tag_multisets) != 35:
        raise LedgerError("ledger requires seven chains and 35 TAG multisets")
    chain_by_id = {item.id: item for item in chains}
    if len(chain_by_id) != 7:
        raise LedgerError("chain IDs are not unique")
    compartment = chains[0].compartment
    coa_formula, coa_charge = chains[0].coa_formula, chains[0].coa_charge
    if any(
        item.compartment != compartment
        or item.coa_formula != coa_formula
        or item.coa_charge != coa_charge
        for item in chains
    ):
        raise LedgerError("all chains must share compartment and CoA")

    model = Model("lipid_moiety_ledger_dry_run")
    try:
        model.solver = "glpk"
    except Exception as error:
        raise LedgerError("GLPK is required for dry-run FBA probes") from error

    transitions: list[TransitionSpec] = []
    lipid_classes: list[LipidClass] = []
    scaffold = _met(
        "LM_G_SCAFFOLD",
        "glycerol scaffold (dry-run source)",
        SCAFFOLD_FORMULA,
        SCAFFOLD_CHARGE,
        compartment,
    )
    coa = _met(
        "LM_COA_SOURCE",
        "coenzyme A (source proxy)",
        coa_formula,
        coa_charge,
        compartment,
    )
    acyls = {
        item.id: _met(
            f"LM_ACYL_SOURCE__{item.id}",
            f"{item.curation_member} (source proxy)",
            item.acyl_coa_formula,
            item.acyl_coa_charge,
            compartment,
        )
        for item in chains
    }
    mono = _met(
        "LM_MONO_BACKBONE",
        "one-acyl glycerol backbone",
        SCAFFOLD_FORMULA,
        SCAFFOLD_CHARGE,
        compartment,
    )
    di = _met(
        "LM_DI_BACKBONE",
        "two-acyl glycerol backbone",
        SCAFFOLD_FORMULA,
        SCAFFOLD_CHARGE,
        compartment,
    )
    lipid_classes.extend(
        [
            LipidClass(
                mono.id, mono.name, mono.formula, 0, compartment, True, "mono_backbone"
            ),
            LipidClass(di.id, di.name, di.formula, 0, compartment, True, "di_backbone"),
        ]
    )

    ledgers: dict[str, Metabolite] = {}
    ledger_ids: dict[str, str] = {}
    for item in chains:
        ledger = _met(
            f"LM_L__{item.id}",
            f"{item.id} acyl-residue ledger",
            item.residue_formula,
            item.residue_charge,
            compartment,
        )
        ledgers[item.id] = ledger
        ledger_ids[item.id] = ledger.id
        lipid_classes.append(
            LipidClass(
                ledger.id,
                ledger.name,
                ledger.formula,
                int(ledger.charge),
                compartment,
                True,
                "acyl_residue_ledger",
            )
        )

    probability_by_composition = {
        tuple(str(chain) for chain in item["composition"]): {
            "prob": float(item["prob"]),
            "cumulative_coverage": float(item["cumulative_coverage"]),
            "curation_rank": int(item["curation_rank"]),
        }
        for item in curated_tags
    }
    kept_probability = sum(item["prob"] for item in probability_by_composition.values())
    tags: dict[tuple[str, str, str], Metabolite] = {}
    tag_records: list[dict[str, Any]] = []
    for composition in tag_multisets:
        tag_formula = formula_from_elements(
            add_formulas(
                [parse_formula(SCAFFOLD_FORMULA)]
                + [
                    parse_formula(chain_by_id[item].residue_formula)
                    for item in composition
                ]
            )
        )
        tag_charge = sum(chain_by_id[item].residue_charge for item in composition)
        tag = _met(
            _tag_id(composition),
            "TAG(" + "/".join(composition) + ") aggregate output",
            tag_formula,
            tag_charge,
            compartment,
        )
        tags[composition] = tag
        lipid_classes.append(
            LipidClass(
                tag.id,
                tag.name,
                tag.formula,
                int(tag.charge),
                compartment,
                True,
                "tag_output",
            )
        )
        probability = probability_by_composition[composition]
        tag_records.append(
            {
                "id": tag.id,
                "chains": list(composition),
                "formula": tag.formula,
                "charge": int(tag.charge),
                "prob": probability["prob"],
                "normalized_prob": probability["prob"] / kept_probability,
                "normalized_prob_semantics": "conditional_on_retained_35_tag_set",
                "cumulative_coverage": probability["cumulative_coverage"],
                "curation_rank": probability["curation_rank"],
                "terminal_routes": [],
            }
        )

    model.add_metabolites(
        [scaffold, coa, *acyls.values(), mono, di, *ledgers.values(), *tags.values()]
    )

    for item in chains:
        _add_reaction(
            model,
            transitions,
            f"LM_FIRST_ACYLATION__{item.id}",
            f"first acylation with {item.id}",
            "first_acylation",
            {
                scaffold: -1,
                acyls[item.id]: -1,
                mono: 1,
                coa: 1,
                ledgers[item.id]: 1,
            },
            (item.id,),
        )
    for item in chains:
        _add_reaction(
            model,
            transitions,
            f"LM_SECOND_ACYLATION__{item.id}",
            f"second acylation with {item.id}",
            "second_acylation",
            {
                mono: -1,
                acyls[item.id]: -1,
                di: 1,
                coa: 1,
                ledgers[item.id]: 1,
            },
            (item.id,),
        )

    record_by_id = {record["id"]: record for record in tag_records}
    for composition in tag_multisets:
        for last_chain in sorted(set(composition)):
            prior = list(composition)
            prior.remove(last_chain)
            stoich: dict[Metabolite, int] = {
                di: -1,
                acyls[last_chain]: -1,
                tags[composition]: 1,
                coa: 1,
            }
            for prior_chain in prior:
                stoich[ledgers[prior_chain]] = stoich.get(ledgers[prior_chain], 0) - 1
            identifier = _route_id(composition, last_chain)
            _add_reaction(
                model,
                transitions,
                identifier,
                "final TAG route for " + "/".join(composition) + f"; last={last_chain}",
                "terminal_tag_route",
                stoich,
                composition,
            )
            record_by_id[tags[composition].id]["terminal_routes"].append(identifier)

    lipid_classes.sort(key=lambda item: item.id)
    transitions.sort(key=lambda item: item.id)
    tag_records.sort(key=lambda item: item["id"])
    for record in tag_records:
        record["terminal_routes"].sort()

    source_bindings = {
        "scaffold": {
            "proxy_id": scaffold.id,
            "source_model_id": None,
            "formula": scaffold.formula,
            "charge": int(scaffold.charge),
        },
        "coa": {
            "proxy_id": coa.id,
            "source_model_id": chains[0].coa_id,
            "formula": coa.formula,
            "charge": int(coa.charge),
        },
        "acyl_inputs": [
            {
                "chain": item.id,
                "proxy_id": acyls[item.id].id,
                "source_model_id": item.acyl_coa_id,
                "formula": acyls[item.id].formula,
                "charge": int(acyls[item.id].charge),
            }
            for item in chains
        ],
    }
    return (
        model,
        tuple(lipid_classes),
        tuple(transitions),
        dict(sorted(ledger_ids.items())),
        tag_records,
        source_bindings,
    )


def build_candidate_gpr_evaluation_model(model: Model) -> Model:
    """Return a labelled copy with unapproved family-level GPR candidates applied.

    This is for structural evaluation only.  It must never be used as the
    canonical ledger model or interpreted as an activation-ready annotation.
    """

    candidate = model.copy()
    for reaction in candidate.reactions:
        if reaction.annotation.get("lipid_moiety_ledger") != "read_only_dry_run":
            continue
        family = str(reaction.annotation.get("ledger_kind", ""))
        assignment = _CANDIDATE_GPRS.get(family)
        if assignment is None:
            raise LedgerError(f"missing candidate GPR family for {reaction.id}")
        reaction.gene_reaction_rule = str(assignment["gpr"])
        reaction.annotation["candidate_gpr_evaluation"] = "true"
    return candidate


def _validate_mass_balance(
    model: Model, transitions: Sequence[TransitionSpec]
) -> dict[str, Any]:
    errors: dict[str, dict[str, float]] = {}
    for transition in transitions:
        residual = model.reactions.get_by_id(transition.id).check_mass_balance()
        if residual:
            errors[transition.id] = {
                key: float(value) for key, value in sorted(residual.items())
            }
    report = {
        "passed": not errors,
        "reaction_count": len(transitions),
        "unbalanced_reactions": errors,
    }
    if errors:
        raise LedgerError(f"mass/charge balance failure: {errors!r}")
    return report


def _validate_occupancy(
    transitions: Sequence[TransitionSpec], ledger_ids: Mapping[str, str]
) -> dict[str, Any]:
    weights = {
        **{identifier: 1 for identifier in ledger_ids.values()},
        "LM_MONO_BACKBONE": -1,
        "LM_DI_BACKBONE": -2,
    }
    residuals: dict[str, int] = {}
    for transition in transitions:
        residual = sum(
            weights.get(identifier, 0) * coefficient
            for identifier, coefficient in transition.stoichiometry
        )
        if residual:
            residuals[transition.id] = residual
    report = {
        "passed": not residuals,
        "vector": "sum(L_i) - mono_backbone - 2*di_backbone",
        "nonzero_residuals": residuals,
    }
    if residuals:
        raise LedgerError(f"occupancy invariant failure: {residuals!r}")
    return report


def _validate_chain_identity(
    transitions: Sequence[TransitionSpec],
    chains: Sequence[Chain],
    ledger_ids: Mapping[str, str],
    tag_records: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    """Seven vectors include source acyl, ledger, and TAG multiplicities."""

    vectors: dict[str, Any] = {}
    all_residuals: dict[str, dict[str, int]] = {}
    for chain in chains:
        weights: dict[str, int] = {
            f"LM_ACYL_SOURCE__{chain.id}": 1,
            ledger_ids[chain.id]: 1,
        }
        tag_multiplicities: dict[str, int] = {}
        for tag in tag_records:
            count = list(tag["chains"]).count(chain.id)
            if count:
                weights[str(tag["id"])] = count
                tag_multiplicities[str(tag["id"])] = count
        residuals: dict[str, int] = {}
        for transition in transitions:
            residual = sum(
                weights.get(identifier, 0) * coefficient
                for identifier, coefficient in transition.stoichiometry
            )
            if residual:
                residuals[transition.id] = residual
        all_residuals[chain.id] = residuals
        vectors[chain.id] = {
            "acyl_input": f"LM_ACYL_SOURCE__{chain.id}",
            "ledger": ledger_ids[chain.id],
            "tag_multiplicities": dict(sorted(tag_multiplicities.items())),
        }
    report = {
        "passed": not any(all_residuals.values()),
        "vectors": dict(sorted(vectors.items())),
        "nonzero_residuals": dict(sorted(all_residuals.items())),
    }
    if not report["passed"]:
        raise LedgerError(f"chain identity invariant failure: {all_residuals!r}")
    return report


def _validate_ids(
    model: Model, transitions: Sequence[TransitionSpec]
) -> dict[str, Any]:
    reaction_duplicates = sorted(
        key
        for key, count in Counter(item.id for item in model.reactions).items()
        if count > 1
    )
    metabolite_duplicates = sorted(
        key
        for key, count in Counter(item.id for item in model.metabolites).items()
        if count > 1
    )
    transition_duplicates = sorted(
        key
        for key, count in Counter(item.id for item in transitions).items()
        if count > 1
    )
    report = {
        "passed": not (
            reaction_duplicates or metabolite_duplicates or transition_duplicates
        ),
        "duplicate_reaction_ids": reaction_duplicates,
        "duplicate_metabolite_ids": metabolite_duplicates,
        "duplicate_transition_ids": transition_duplicates,
    }
    if not report["passed"]:
        raise LedgerError(f"duplicate generated IDs: {report!r}")
    return report


def _validate_budget(
    classes: Sequence[LipidClass],
    transitions: Sequence[TransitionSpec],
    tags: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    first = sum(item.kind == "first_acylation" for item in transitions)
    second = sum(item.kind == "second_acylation" for item in transitions)
    terminal = sum(item.kind == "terminal_tag_route" for item in transitions)
    generated = sum(item.generated for item in classes)
    route_count = sum(len(item["terminal_routes"]) for item in tags)
    report = {
        "passed": (
            first == EXPECTED_CHAIN_COUNT
            and second == EXPECTED_CHAIN_COUNT
            and terminal == EXPECTED_TERMINAL_ROUTES
            and len(transitions) == EXPECTED_REACTION_COUNT
            and len(tags) == EXPECTED_TAG_COUNT
            and route_count == EXPECTED_TERMINAL_ROUTES
            and generated == EXPECTED_GENERATED_METABOLITE_COUNT
        ),
        "first_acylation_reactions": first,
        "second_acylation_reactions": second,
        "terminal_tag_routes": terminal,
        "biochemical_reactions": len(transitions),
        "tag_outputs": len(tags),
        "generated_metabolites": generated,
    }
    if not report["passed"]:
        raise LedgerError(f"architecture budget failure: {report!r}")
    return report


def _boundary(
    identifier: str, metabolite: Metabolite, coefficient: int, upper: float = 1000.0
) -> Reaction:
    reaction = Reaction(identifier, name=identifier)
    reaction.bounds = (0.0, upper)
    reaction.annotation["lipid_moiety_ledger_probe"] = True
    reaction.add_metabolites({metabolite: coefficient})
    return reaction


def _solve(model: Model) -> tuple[str, float]:
    solution = model.optimize()
    return str(solution.status), float(solution.objective_value)


def _near(value: float, expected: float, tolerance: float = 1e-8) -> bool:
    return math.isfinite(value) and abs(value - expected) <= tolerance


def _run_fba_probes(
    model: Model,
    chains: Sequence[Chain],
    ledger_ids: Mapping[str, str],
    tags: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    """FBA proof: required chains cannot substitute; no TAG/ledger is free."""

    probe = model.copy()
    temporary: list[Reaction] = [
        _boundary(
            "LM_PROBE_SRC_G", probe.metabolites.get_by_id("LM_G_SCAFFOLD"), 1, 1.0
        ),
        _boundary("LM_PROBE_DM_COA", probe.metabolites.get_by_id("LM_COA_SOURCE"), -1),
    ]
    acyl_sources: dict[str, Reaction] = {}
    for chain in chains:
        source = _boundary(
            f"LM_PROBE_SRC_ACYL__{chain.id}",
            probe.metabolites.get_by_id(f"LM_ACYL_SOURCE__{chain.id}"),
            1,
            3.0,
        )
        acyl_sources[chain.id] = source
        temporary.append(source)
    tag_demands: dict[str, Reaction] = {}
    for tag in tags:
        demand = _boundary(
            f"LM_PROBE_DM__{tag['id']}",
            probe.metabolites.get_by_id(str(tag["id"])),
            -1,
        )
        tag_demands[str(tag["id"])] = demand
        temporary.append(demand)
    ledger_demands: dict[str, Reaction] = {}
    for chain_id, ledger_id in sorted(ledger_ids.items()):
        demand = _boundary(
            f"LM_PROBE_DM_LEDGER__{chain_id}",
            probe.metabolites.get_by_id(ledger_id),
            -1,
        )
        ledger_demands[chain_id] = demand
        temporary.append(demand)
    probe.add_reactions(temporary)

    positive_controls: dict[str, dict[str, Any]] = {}
    chain_cuts: dict[str, dict[str, Any]] = {}
    pfba_controls: dict[str, dict[str, Any]] = {}
    pfba_statuses: dict[str, str] = {}
    for tag in tags:
        tag_id = str(tag["id"])
        demand = tag_demands[tag_id]
        probe.objective = demand
        status, value = _solve(probe)
        positive_controls[tag_id] = {
            "status": status,
            "objective": value,
            "passed": status == "optimal" and _near(value, 1.0),
        }
        try:
            pfba_solution = pfba(probe)
        except Exception as error:
            raise LedgerError(f"pFBA probe failed for {tag_id}: {error}") from error
        pfba_status = str(pfba_solution.status)
        pfba_flux = float(pfba_solution.fluxes[demand.id])
        pfba_statuses[tag_id] = pfba_status
        pfba_controls[tag_id] = {
            "status": pfba_status,
            "demand_flux": pfba_flux,
            "passed": pfba_status == "optimal" and _near(pfba_flux, 1.0),
        }
        for chain_id in sorted(set(tag["chains"])):
            source = acyl_sources[chain_id]
            prior_upper = source.upper_bound
            source.upper_bound = 0.0
            cut_status, cut_value = _solve(probe)
            source.upper_bound = prior_upper
            key = f"{tag_id}:{chain_id}"
            chain_cuts[key] = {
                "status": cut_status,
                "objective": cut_value,
                "passed": cut_status == "optimal" and _near(cut_value, 0.0),
            }

    original_upper = {key: source.upper_bound for key, source in acyl_sources.items()}
    for source in acyl_sources.values():
        source.upper_bound = 0.0
    all_tag_values: dict[str, float] = {}
    all_tag_statuses: dict[str, str] = {}
    for tag_id, demand in sorted(tag_demands.items()):
        probe.objective = demand
        status, value = _solve(probe)
        all_tag_statuses[tag_id] = status
        all_tag_values[tag_id] = value
    free_ledger_values: dict[str, float] = {}
    free_ledger_statuses: dict[str, str] = {}
    for chain_id, demand in sorted(ledger_demands.items()):
        probe.objective = demand
        status, value = _solve(probe)
        free_ledger_statuses[chain_id] = status
        free_ledger_values[chain_id] = value

    internal_limits: dict[str, dict[str, Any]] = {}
    for reaction_id in sorted(
        reaction.id
        for reaction in probe.reactions
        if reaction.annotation.get("lipid_moiety_ledger") == "read_only_dry_run"
    ):
        reaction = probe.reactions.get_by_id(reaction_id)
        probe.objective = reaction
        probe.objective_direction = "max"
        forward_status, forward_value = _solve(probe)
        entry: dict[str, Any] = {
            "forward_status": forward_status,
            "forward_maximum": forward_value,
            "reverse_status": None,
            "reverse_maximum": 0.0,
        }
        if reaction.lower_bound < 0.0:
            probe.objective_direction = "min"
            reverse_status, reverse_value = _solve(probe)
            entry["reverse_status"] = reverse_status
            entry["reverse_maximum"] = -reverse_value
        probe.objective_direction = "max"
        entry["passed"] = (
            entry["forward_status"] == "optimal"
            and _near(float(entry["forward_maximum"]), 0.0)
            and (
                entry["reverse_status"] is None
                or (
                    entry["reverse_status"] == "optimal"
                    and _near(float(entry["reverse_maximum"]), 0.0)
                )
            )
        )
        internal_limits[reaction_id] = entry
    for chain_id, upper in original_upper.items():
        acyl_sources[chain_id].upper_bound = upper

    passed = (
        all(item["passed"] for item in positive_controls.values())
        and all(item["passed"] for item in pfba_controls.values())
        and all(item["passed"] for item in chain_cuts.values())
        and all(
            all_tag_statuses[tag_id] == "optimal" and _near(value, 0.0)
            for tag_id, value in all_tag_values.items()
        )
        and all(
            free_ledger_statuses[chain_id] == "optimal" and _near(value, 0.0)
            for chain_id, value in free_ledger_values.items()
        )
        and len(internal_limits) == EXPECTED_REACTION_COUNT
        and all(item["passed"] for item in internal_limits.values())
    )
    report = {
        "passed": passed,
        "positive_controls": dict(sorted(positive_controls.items())),
        "pfba_statuses": dict(sorted(pfba_statuses.items())),
        "pfba_controls": dict(sorted(pfba_controls.items())),
        "required_chain_cuts": dict(sorted(chain_cuts.items())),
        "all_acyl_closed_tag_demands": dict(sorted(all_tag_values.items())),
        "all_acyl_closed_tag_statuses": dict(sorted(all_tag_statuses.items())),
        "all_acyl_closed_free_ledger_demands": dict(sorted(free_ledger_values.items())),
        "all_acyl_closed_free_ledger_statuses": dict(
            sorted(free_ledger_statuses.items())
        ),
        "all_acyl_closed_internal_reaction_maxima": dict(
            sorted(internal_limits.items())
        ),
        "all_acyl_closed_internal_flux_blocked": {
            "passed": len(internal_limits) == EXPECTED_REACTION_COUNT
            and all(item["passed"] for item in internal_limits.values()),
            "reaction_count": len(internal_limits),
        },
    }
    if not passed:
        raise LedgerError(f"FBA probe failure: {report!r}")
    return report


def run_composition_soft_prior(
    model: Model, tags: Sequence[Mapping[str, Any]], *, retention: float = 0.99
) -> dict[str, Any]:
    """Fit optional TAG flux composition without constraining default FBA.

    The probe first maximizes the aggregate TAG demand.  It then retains the
    requested fraction of that optimum and minimizes L1 distance to retained,
    normalized curation probabilities.  It deliberately does not reuse source
    abundance weights as supply bounds.
    """

    if not (0.0 < retention <= 1.0):
        raise LedgerError("soft-prior retention must be in (0, 1]")
    probability_sum = sum(float(tag["normalized_prob"]) for tag in tags)
    if not _near(probability_sum, 1.0, 1e-7):
        raise LedgerError("soft-prior TAG probabilities must sum to one")
    probe = model.copy()
    temporary: list[Reaction] = [
        _boundary("LM_SOFT_SRC_G", probe.metabolites.get_by_id("LM_G_SCAFFOLD"), 1, 1.0),
        _boundary("LM_SOFT_DM_COA", probe.metabolites.get_by_id("LM_COA_SOURCE"), -1),
    ]
    for identifier in sorted(
        met.id for met in probe.metabolites if met.id.startswith("LM_ACYL_SOURCE__")
    ):
        temporary.append(_boundary(f"LM_SOFT_SRC__{identifier}", probe.metabolites.get_by_id(identifier), 1, 3.0))
    demands: dict[str, Reaction] = {}
    for tag in tags:
        tag_id = str(tag["id"])
        demand = _boundary(
            f"LM_SOFT_DM__{tag_id}", probe.metabolites.get_by_id(tag_id), -1
        )
        temporary.append(demand)
        demands[tag_id] = demand
    probe.add_reactions(temporary)
    total_expression = sum(demand.flux_expression for demand in demands.values())
    probe.objective = total_expression
    probe.objective_direction = "max"
    primary = probe.optimize()
    if primary.status != "optimal":
        raise LedgerError(f"soft-prior primary TAG optimization failed: {primary.status}")
    primary_maximum = sum(float(primary.fluxes[demand.id]) for demand in demands.values())
    if not (math.isfinite(primary_maximum) and primary_maximum > 0.0):
        raise LedgerError("soft-prior primary TAG optimum is not positive and finite")
    problem = probe.problem
    keep_total = problem.Constraint(
        total_expression,
        lb=retention * primary_maximum,
        ub=retention * primary_maximum,
        name="LM_SOFT_FIX_TOTAL_TAG",
    )
    additions: list[Any] = [keep_total]
    deviations: dict[str, tuple[Any, Any]] = {}
    for index, tag in enumerate(sorted(tags, key=lambda item: str(item["id"]))):
        tag_id = str(tag["id"])
        plus = problem.Variable(f"LM_SOFT_DEV_POS__{index}", lb=0.0)
        minus = problem.Variable(f"LM_SOFT_DEV_NEG__{index}", lb=0.0)
        target = float(tag["normalized_prob"])
        equality = problem.Constraint(
            plus - minus - demands[tag_id].flux_expression + target * total_expression,
            lb=0.0,
            ub=0.0,
            name=f"LM_SOFT_MATCH__{index}",
        )
        additions.extend([plus, minus, equality])
        deviations[tag_id] = (plus, minus)
    probe.add_cons_vars(additions)
    probe.objective = problem.Objective(
        sum(plus + minus for plus, minus in deviations.values()), direction="min"
    )
    fitted = probe.optimize()
    if fitted.status != "optimal":
        raise LedgerError(f"soft-prior composition fit failed: {fitted.status}")
    total_flux = sum(float(fitted.fluxes[demand.id]) for demand in demands.values())
    per_tag = {
        tag_id: {
            "realized_flux": float(fitted.fluxes[demands[tag_id].id]),
            "target_flux": float(next(tag for tag in tags if str(tag["id"]) == tag_id)["normalized_prob"]) * total_flux,
            "absolute_deviation": abs(
                float(fitted.fluxes[demands[tag_id].id])
                - float(next(tag for tag in tags if str(tag["id"]) == tag_id)["normalized_prob"]) * total_flux
            ),
        }
        for tag_id in sorted(demands)
    }
    return {
        "mode": "opt_in_soft_prior",
        "retention": retention,
        "primary_status": str(primary.status),
        "primary_tag_maximum": primary_maximum,
        "fit_status": str(fitted.status),
        "realized_total_tag_flux": total_flux,
        "l1_deviation": max(0.0, float(fitted.objective_value)),
        "per_tag": per_tag,
    }


def _update_activation_status(manifest: dict[str, Any]) -> None:
    gates = manifest["gates"]
    required = (
        "cobra_mass_balance",
        "occupancy_invariant",
        "chain_identity",
        "no_duplicate_ids",
        "architecture_budget",
        "fba_no_substitution_or_free_lipid",
        "source_model_fingerprint_unchanged",
        "determinism",
        "source_coa_acyl_formula_charge_verified",
    )
    gates["all_required_gates"] = all(bool(gates[name]) for name in required)
    gates["activation_blockers_open"] = any(
        blocker.get("status") == "open" for blocker in manifest["blockers"]
    )
    manifest["activation_ready"] = False


def _reaction_manifest(
    transitions: Sequence[TransitionSpec],
) -> list[dict[str, Any]]:
    return [
        {
            "id": item.id,
            "label": item.label,
            "kind": item.kind,
            "candidate_gpr": (
                None
                if item.candidate_gpr is None
                else dict(_CANDIDATE_GPRS[item.kind])
            ),
            "required_chains": list(item.required_chains),
            "stoichiometry": [
                {"metabolite": metabolite, "coefficient": coefficient}
                for metabolite, coefficient in item.stoichiometry
            ],
        }
        for item in sorted(transitions, key=lambda item: item.id)
    ]


def _manifest(
    *,
    input_sha256: Mapping[str, str],
    chains: Sequence[Chain],
    classes: Sequence[LipidClass],
    transitions: Sequence[TransitionSpec],
    ledger_ids: Mapping[str, str],
    tags: Sequence[Mapping[str, Any]],
    source_proxies: Mapping[str, Any],
    validation: Mapping[str, Any],
) -> dict[str, Any]:
    flat_validation = dict(validation)
    source_chemistry = _required_mapping(
        flat_validation["source_chemical_convention"],
        "validation.source_chemical_convention",
    )
    source_verified = bool(source_chemistry["passed"])
    blockers: list[dict[str, Any]] = [
        {
            "id": blocker["id"],
            "status": "open",
            "severity": blocker["severity"],
            "reason": blocker["effect"],
        }
        for blocker in _SPEC_BLOCKERS
        if blocker["id"] != "chemistry_source_not_normalized"
    ]
    if not source_verified:
        source_blocker = next(
            blocker
            for blocker in _SPEC_BLOCKERS
            if blocker["id"] == "chemistry_source_not_normalized"
        )
        blockers.append(
            {
                "id": source_blocker["id"],
                "status": "open",
                "severity": source_blocker["severity"],
                "reason": source_blocker["effect"],
                "source_chemistry_errors": list(source_chemistry["errors"]),
            }
        )
    fba_probe = _required_mapping(flat_validation["fba_probe"], "validation.fba_probe")
    positive_controls = _required_mapping(
        fba_probe["positive_controls"], "validation.fba_probe.positive_controls"
    )
    pfba_controls = _required_mapping(
        fba_probe["pfba_controls"], "validation.fba_probe.pfba_controls"
    )
    flat_validation.update(
        {
            "generated_biochemical_reactions_balanced": bool(
                _required_mapping(
                    flat_validation["mass_balance"], "validation.mass_balance"
                )["passed"]
            ),
            "generated_metabolites_within_budget": bool(
                _required_mapping(flat_validation["budget"], "validation.budget")[
                    "passed"
                ]
            ),
            "chain_identity_conserved": bool(
                _required_mapping(
                    flat_validation["chain_identity"], "validation.chain_identity"
                )["passed"]
            ),
            "occupancy_conserved": bool(
                _required_mapping(flat_validation["occupancy"], "validation.occupancy")[
                    "passed"
                ]
            ),
            "source_coa_acyl_formula_charge_verified": source_verified,
            "per_tag_required_chain_cut": bool(
                all(
                    item.get("passed") is True
                    for item in _required_mapping(
                        fba_probe["required_chain_cuts"],
                        "validation.fba_probe.required_chain_cuts",
                    ).values()
                )
            ),
            "all_supply_cut": bool(
                all(
                    status == "optimal" and _near(float(value), 0.0)
                    for tag_id, value in _required_mapping(
                        fba_probe["all_acyl_closed_tag_demands"],
                        "validation.fba_probe.all_acyl_closed_tag_demands",
                    ).items()
                    for status in [
                        _required_mapping(
                            fba_probe["all_acyl_closed_tag_statuses"],
                            "validation.fba_probe.all_acyl_closed_tag_statuses",
                        )[tag_id]
                    ]
                )
                and all(
                    status == "optimal" and _near(float(value), 0.0)
                    for chain_id, value in _required_mapping(
                        fba_probe["all_acyl_closed_free_ledger_demands"],
                        "validation.fba_probe.all_acyl_closed_free_ledger_demands",
                    ).items()
                    for status in [
                        _required_mapping(
                            fba_probe["all_acyl_closed_free_ledger_statuses"],
                            "validation.fba_probe.all_acyl_closed_free_ledger_statuses",
                        )[chain_id]
                    ]
                )
            ),
            "fba_status": "optimal"
            if all(item.get("passed") is True for item in positive_controls.values())
            else "failed",
            "pfba_status": "optimal"
            if all(item.get("passed") is True for item in pfba_controls.values())
            else "failed",
            "source_model_unchanged": bool(flat_validation["source_model_unchanged"]),
        }
    )
    manifest: dict[str, Any] = {
        "schema_version": 2,
        "mode": "read_only_in_memory_manifest_only",
        "activation_ready": False,
        "input_sha256": dict(sorted(input_sha256.items())),
        "architecture": {
            "style": "global_acyl_moiety_ledger",
            "first_acylation_reactions": EXPECTED_CHAIN_COUNT,
            "second_acylation_reactions": EXPECTED_CHAIN_COUNT,
            "terminal_tag_routes": EXPECTED_TERMINAL_ROUTES,
            "biochemical_reactions_exact": EXPECTED_REACTION_COUNT,
            "generated_metabolites_exact": EXPECTED_GENERATED_METABOLITE_COUNT,
            "generated_metabolites_maximum": EXPECTED_GENERATED_METABOLITE_COUNT,
            "carry_reactions": 0,
            "ticket_metabolites": 0,
        },
        "source_chemical_convention": dict(source_chemistry),
        "counts": {
            "chains": len(chains),
            "tag_outputs": len(tags),
            "terminal_tag_routes": sum(len(tag["terminal_routes"]) for tag in tags),
            "biochemical_reactions": len(transitions),
            "generated_metabolites": sum(item.generated for item in classes),
        },
        "chains": [
            {
                "id": item.id,
                "curation_member": item.curation_member,
                "source_model_acyl_coa_id": item.acyl_coa_id,
                "source_model_coa_id": item.coa_id,
                "compartment": item.compartment,
                "acyl_coa_formula": item.acyl_coa_formula,
                "acyl_coa_charge": item.acyl_coa_charge,
                "coa_formula": item.coa_formula,
                "coa_charge": item.coa_charge,
                "residue_formula": item.residue_formula,
                "residue_charge": item.residue_charge,
                "ledger_tuple_origin": source_chemistry["tuple_origin"],
            }
            for item in sorted(chains, key=lambda item: item.id)
        ],
        "source_proxies": source_proxies,
        "chain_ledger_ids": dict(sorted(ledger_ids.items())),
        "generated_metabolites": [
            {
                "id": item.id,
                "label": item.label,
                "kind": item.kind,
                "formula": item.formula,
                "charge": item.charge,
                "compartment": item.compartment,
            }
            for item in classes
            if item.generated
        ],
        "tag_outputs": [dict(item) for item in tags],
        "abundance": dict(_required_mapping(flat_validation["abundance"], "validation.abundance")),
        "candidate_gpr": {
            "reaction_count": len(transitions),
            "assigned_reaction_count": sum(item.candidate_gpr is not None for item in transitions),
            "families": dict(_CANDIDATE_GPRS),
            "activation_status": "blocked_pending_compartment_and_substrate_evidence",
        },
        "generated_reactions": _reaction_manifest(transitions),
        "validation": flat_validation,
        "gates": {
            "cobra_mass_balance": bool(validation["mass_balance"]["passed"]),
            "occupancy_invariant": bool(validation["occupancy"]["passed"]),
            "chain_identity": bool(validation["chain_identity"]["passed"]),
            "no_duplicate_ids": bool(validation["ids"]["passed"]),
            "architecture_budget": bool(validation["budget"]["passed"]),
            "fba_no_substitution_or_free_lipid": bool(
                validation["fba_probe"]["passed"]
            ),
            "source_model_fingerprint_unchanged": bool(
                validation["source_model_unchanged"]
            ),
            "source_coa_acyl_formula_charge_verified": source_verified,
            "determinism": False,
        },
        "blockers": blockers,
        "capability_loss": {
            "retained": "Aggregate histogram of seven acyl-chain identities.",
            "blocked_scopes": [
                {
                    "id": "unsupported_transport_scope",
                    "status": "blocked",
                    "reason": "No cross-compartment provenance is represented.",
                },
                {
                    "id": "unsupported_remodeling_scope",
                    "status": "blocked",
                    "reason": "No DAG pairing or selective-lipolysis identity is retained.",
                },
                {
                    "id": "unsupported_topology_scope",
                    "status": "blocked",
                    "reason": "No sn-position or site-specific topology is represented.",
                },
            ],
        },
        "execution_boundary": {
            "source_model_immutable": True,
            "sbml_output_forbidden": True,
            "model_activation_forbidden": True,
            "allowed_artifact": "canonical_json_manifest_only",
        },
    }
    _update_activation_status(manifest)
    return manifest


def _read_model(path: Path) -> Model:
    try:
        cobra.Configuration().solver = "glpk"
    except Exception:
        pass
    with (
        contextlib.redirect_stdout(io.StringIO()),
        contextlib.redirect_stderr(io.StringIO()),
    ):
        try:
            return cobra.io.read_sbml_model(str(path))
        except Exception as error:
            raise LedgerError(f"could not read source SBML: {path}") from error


def _resolve_spec_path(spec_path: str | Path | None) -> Path:
    if spec_path is None:
        return (
            Path(__file__).resolve().parents[1]
            / "data"
            / "lipid_moiety_ledger_spec.json"
        )
    return Path(spec_path)


def _validated_source_objective(source: Model) -> dict[str, Any]:
    """Require the declared source baseline to optimize biomass_C and nothing else."""

    biomass_id = "biomass_C"
    if not source.reactions.has_id(biomass_id):
        raise LedgerError(f"source model lacks required objective reaction {biomass_id}")
    direction = str(source.objective.direction)
    if direction != "max":
        raise LedgerError(f"source objective direction must be max, found {direction!r}")
    coefficients: dict[str, float] = {}
    for reaction, coefficient in linear_reaction_coefficients(source).items():
        numeric = float(coefficient)
        if not math.isfinite(numeric):
            raise LedgerError(
                f"source objective contains a non-finite coefficient for {reaction.id}"
            )
        # The source objective is an exact contract, not a numerical comparison:
        # even a tiny secondary term changes the stated optimization problem.
        if numeric != 0.0:
            coefficients[reaction.id] = numeric
    if coefficients != {biomass_id: 1.0}:
        raise LedgerError(
            "source objective must be the unique biomass_C coefficient +1 objective; "
            f"found {coefficients!r}"
        )
    return {
        "objective_reaction": biomass_id,
        "objective_direction": direction,
        "objective_coefficient": 1.0,
        "objective_expression": coefficients,
    }


def _independent_source_baseline(source: Model) -> dict[str, Any]:
    """One source-model FBA/pFBA check, explicitly not an integrated regression."""

    objective_contract = _validated_source_objective(source)
    fba_solution = source.optimize()
    if fba_solution.status != "optimal":
        raise LedgerError(f"source model FBA is not optimal: {fba_solution.status}")
    biomass_id = str(objective_contract["objective_reaction"])
    try:
        pfba_solution = pfba(source)
    except Exception as error:
        raise LedgerError(f"source model pFBA failed: {error}") from error
    fba_biomass = float(fba_solution.fluxes[biomass_id])
    pfba_biomass = float(pfba_solution.fluxes[biomass_id])
    if pfba_solution.status != "optimal" or not _near(fba_biomass, pfba_biomass):
        raise LedgerError(
            f"source FBA/pFBA biomass mismatch: {fba_biomass} vs {pfba_biomass}"
        )
    return {
        "passed": True,
        "scope": "independent_source_baseline_not_integrated_regression",
        **objective_contract,
        "fba_status": str(fba_solution.status),
        "fba_biomass_flux": fba_biomass,
        "pfba_status": str(pfba_solution.status),
        "pfba_biomass_flux": pfba_biomass,
    }


def _compile_once(
    model_path: str | Path,
    curation_csv: str | Path,
    *,
    spec_path: str | Path | None,
) -> DryRunResult:
    source_path, curation_path = Path(model_path), Path(curation_csv)
    if source_path.suffix.lower() not in {".xml", ".sbml"}:
        raise LedgerError("source model must be .xml or .sbml")
    input_sha256 = {
        "model.xml": sha256_file(source_path),
        "data/lipid_combo_curation.csv": sha256_file(curation_path),
    }
    resolved_spec = _resolve_spec_path(spec_path)
    input_sha256["data/lipid_moiety_ledger_spec.json"] = sha256_file(resolved_spec)
    spec = load_spec(resolved_spec)
    _validate_public_compile_contract(spec)
    source_chemical_convention = _validated_source_chemical_convention(spec)
    source = _read_model(source_path)
    source_baseline = _independent_source_baseline(source)
    source_bindings = _source_pool_bindings(
        source, source_chemical_convention, POOL_REACTION_ID
    )
    source_chemistry = _source_chemical_verification(
        source_bindings, source_chemical_convention
    )
    if source_chemistry["passed"]:
        chains = load_chains(source, source_chemical_convention)
    else:
        # This is an explicitly labelled curation tuple, not a source repair.
        # It keeps the read-only 95-reaction feasibility exercise available while
        # the immutable source mismatch remains an activation-blocking finding.
        chains = _chains_from_bindings(source_bindings, source_chemical_convention)
    curated_tags = load_curated_tag_records(curation_path, chains)
    tags = tuple(item["composition"] for item in curated_tags)
    source_weights = _source_pool_abundance_weights(source, chains)
    _validate_spec_against_plan(spec, chains, tags)
    model, classes, transitions, ledger_ids, tag_records, source_proxies = (
        _build_ledger(chains, curated_tags)
    )
    validation: dict[str, Any] = {
        "mass_balance": _validate_mass_balance(model, transitions),
        "occupancy": _validate_occupancy(transitions, ledger_ids),
        "chain_identity": _validate_chain_identity(
            transitions, chains, ledger_ids, tag_records
        ),
        "ids": _validate_ids(model, transitions),
        "budget": _validate_budget(classes, transitions, tag_records),
        "fba_probe": _run_fba_probes(model, chains, ledger_ids, tag_records),
        "abundance": _abundance_report(tag_records, source_weights),
        "source_chemical_convention": source_chemistry,
        "independent_source_baseline": source_baseline,
        "source_model_baseline": source_baseline,
        "source_model_unchanged": False,
        "curation_unchanged": False,
        "spec_unchanged": False,
        "determinism": {"passed": False, "status": "deferred"},
    }
    if input_sha256["model.xml"] != sha256_file(source_path):
        raise LedgerError(
            "source model fingerprint changed during read-only compilation"
        )
    if input_sha256["data/lipid_combo_curation.csv"] != sha256_file(curation_path):
        raise LedgerError("curation fingerprint changed during compilation")
    if input_sha256["data/lipid_moiety_ledger_spec.json"] != sha256_file(resolved_spec):
        raise LedgerError("spec fingerprint changed during compilation")
    validation["source_model_unchanged"] = True
    validation["curation_unchanged"] = True
    validation["spec_unchanged"] = True
    manifest = _manifest(
        input_sha256=input_sha256,
        chains=chains,
        classes=classes,
        transitions=transitions,
        ledger_ids=ledger_ids,
        tags=tag_records,
        source_proxies=source_proxies,
        validation=validation,
    )
    manifest["gates"]["source_model_fba_pfba"] = bool(source_baseline["passed"])
    return DryRunResult(
        model=model,
        chains=chains,
        lipid_classes=classes,
        transitions=transitions,
        manifest=manifest,
        canonical_json=canonical_json_bytes(manifest),
    )


def compile_dry_run(
    model_path: str | Path,
    curation_csv: str | Path,
    *,
    spec_path: str | Path | None = None,
    verify_determinism: bool = True,
) -> DryRunResult:
    """Compile once with solvers, then validate structure twice without solvers."""

    first = _compile_once(model_path, curation_csv, spec_path=spec_path)
    input_paths = {
        "model.xml": Path(model_path),
        "data/lipid_combo_curation.csv": Path(curation_csv),
        "data/lipid_moiety_ledger_spec.json": _resolve_spec_path(spec_path),
    }
    input_sha256 = first.manifest["input_sha256"]
    curated_tags = tuple(
        {
            "composition": tuple(str(chain) for chain in tag["chains"]),
            "prob": float(tag["prob"]),
            "cumulative_coverage": float(tag["cumulative_coverage"]),
            "curation_rank": int(tag["curation_rank"]),
        }
        for tag in first.manifest["tag_outputs"]
    )

    determinism = {
        "passed": bool(verify_determinism),
        "status": (
            "two_fresh_compilations_byte_identical"
            if verify_determinism
            else "not_requested"
        ),
    }

    def final_manifest(
        classes: Sequence[LipidClass],
        transitions: Sequence[TransitionSpec],
        ledger_ids: Mapping[str, str],
        tags: Sequence[Mapping[str, Any]],
        source_proxies: Mapping[str, Any],
    ) -> dict[str, Any]:
        validation = dict(first.manifest["validation"])
        validation["determinism"] = dict(determinism)
        manifest = _manifest(
            input_sha256=input_sha256,
            chains=first.chains,
            classes=classes,
            transitions=transitions,
            ledger_ids=ledger_ids,
            tags=tags,
            source_proxies=source_proxies,
            validation=validation,
        )
        manifest["gates"]["source_model_fba_pfba"] = bool(
            validation["independent_source_baseline"]["passed"]
        )
        manifest["gates"]["determinism"] = bool(determinism["passed"])
        _update_activation_status(manifest)
        return manifest

    first_manifest = final_manifest(
        first.lipid_classes,
        first.transitions,
        first.manifest["chain_ledger_ids"],
        first.manifest["tag_outputs"],
        first.manifest["source_proxies"],
    )
    if verify_determinism:
        (
            _second_model,
            second_classes,
            second_transitions,
            second_ledger_ids,
            second_tags,
            second_source_proxies,
        ) = _build_ledger(
            first.chains, curated_tags
        )  # pure structural build: no solve
        second_manifest = final_manifest(
            second_classes,
            second_transitions,
            second_ledger_ids,
            second_tags,
            second_source_proxies,
        )
        if canonical_json_bytes(first_manifest) != canonical_json_bytes(
            second_manifest
        ):
            raise LedgerError("fresh structural final manifests are not byte-identical")

    validate_input_fingerprints(input_paths, input_sha256)
    first.manifest = first_manifest
    first.canonical_json = canonical_json_bytes(first_manifest)
    return first


def build_dry_run(
    model_path: str | Path,
    curation_csv: str | Path,
    *,
    spec_path: str | Path | None = None,
    verify_determinism: bool = True,
) -> DryRunResult:
    """Alias for the public compiler."""

    return compile_dry_run(
        model_path,
        curation_csv,
        spec_path=spec_path,
        verify_determinism=verify_determinism,
    )


def plan_lipid_moiety_ledger(
    model_path: str | Path,
    curation_csv: str | Path,
    *,
    spec_path: str | Path | None = None,
    verify_determinism: bool = True,
) -> DryRunResult:
    """Planner-oriented alias for the public compiler."""

    return compile_dry_run(
        model_path,
        curation_csv,
        spec_path=spec_path,
        verify_determinism=verify_determinism,
    )


__all__ = [
    "Chain",
    "DryRunResult",
    "LedgerError",
    "LipidClass",
    "TransitionSpec",
    "add_formulas",
    "atomic_write_manifest",
    "build_candidate_gpr_evaluation_model",
    "build_dry_run",
    "canonical_json_bytes",
    "canonicalize_tag_multiset",
    "compile_dry_run",
    "formula_from_elements",
    "load_chains",
    "load_curated_tag_records",
    "load_curated_tag_multisets",
    "load_spec",
    "parse_formula",
    "plan_lipid_moiety_ledger",
    "run_composition_soft_prior",
    "sha256_file",
    "subtract_formula",
    "validate_input_fingerprints",
    "validate_manifest_output_path",
]
