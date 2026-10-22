"""Opt-in thiamine phosphate chemistry: KEGG formulas and IDs for TMP, TPP, TTP and HET-P, and R1342's EC."""

import copy
import json
from xml.sax.saxutils import unescape

from .config import MODEL

CURATION_PATH = MODEL.curation_file("thiamine_phosphate_chemistry.json")
SPECIES_IDS = {"m977[C_cy]", "m973[C_ex]", "m975[C_cy]", "m971[C_ex]", "m1116[C_mi]", "m976[C_cy]", "m643[C_cy]"}
REACTION_IDS = {"R1342"}
# Every reaction on the species; the balance review assumed exactly these equations.
CONTEXT_IDS = {"R1337", "R1339", "R1341", "R1342", "R1343", "R1344", "R408", "R750", "R751", "R752", "R753",
               "R754", "R933", "biomass_C"}
TPP = "m975[C_cy]"
XML_ENTITIES = {"&apos;": "'", "&quot;": '"'}
# The reviewed chemistry (KEGG fully protonated form, ChEBI parent charge), pinned so an edited data file
# cannot change it silently.
REVIEWED = {
    "thiamine monophosphate": ("C12H18N4O4PS", 1, "C01081", "CHEBI:9533", "MNXM662", "thmmp"),
    "thiamine diphosphate": ("C12H19N4O7P2S", 1, "C00068", "CHEBI:9532", "MNXM256", "thmpp"),
    "thiamine triphosphate": ("C12H20N4O10P3S", 1, "C03028", "CHEBI:9534", "MNXM1341", "thmtp"),
    "4-methyl-5-(2-phosphonooxyethyl)thiazole": ("C6H10NO4PS", 0, "C04327", "CHEBI:17857", "MNXM960", "4mpetz"),
}
MOLECULE_OF = {"m977[C_cy]": "thiamine monophosphate", "m973[C_ex]": "thiamine monophosphate",
               "m975[C_cy]": "thiamine diphosphate", "m971[C_ex]": "thiamine diphosphate",
               "m1116[C_mi]": "thiamine diphosphate", "m976[C_cy]": "thiamine triphosphate",
               "m643[C_cy]": "4-methyl-5-(2-phosphonooxyethyl)thiazole"}


def _unescaped(notes):
    # Notes may read back from SBML with &gt; and &quot; where a build holds > and ".
    return {key: unescape(value, XML_ENTITIES) for key, value in notes.items()}


def _annotation(annotation):
    # A single cross-reference reads back from SBML as a string, several as a list.
    return {key: sorted([value] if isinstance(value, str) else value) for key, value in annotation.items()}


def _notes(spec, row):
    return {"thiamine_phosphate_chemistry_curation_id": spec["curation_id"],
            "thiamine_phosphate_chemistry_status": spec["status"],
            "thiamine_phosphate_chemistry_authorization": spec["authorization"],
            "thiamine_phosphate_chemistry_evidence": row["evidence"],
            "thiamine_phosphate_chemistry_evidence_tier": row["evidence_tier"],
            "thiamine_phosphate_chemistry_rationale": row["rationale"],
            "thiamine_phosphate_chemistry_previous_state": row["previous_state"]}


def _check_spec(spec):
    version = spec.get("schema_version")
    if (type(version) is not int or version != 1
            or spec.get("curation_id") != "THIAMINE-P-CHEMISTRY-20261022"
            or spec.get("status") != "user_authorized_formula_charge_annotation_correction"
            or set(spec.get("species", {})) != SPECIES_IDS
            or set(spec.get("reactions", {})) != REACTION_IDS
            or set(spec.get("context_reactions", {})) != CONTEXT_IDS
            or not isinstance(spec.get("partner_species"), dict)):
        raise ValueError("Unexpected thiamine phosphate chemistry scope/status")
    # The reviewed change: exactly the pinned chemistry per species, and only R1342's EC.
    for mid, row in spec["species"].items():
        formula, charge, kegg, chebi, mnx, bigg = REVIEWED[MOLECULE_OF[mid]]
        after = row.get("after", {})
        if (row.get("molecule") != MOLECULE_OF[mid] or after.get("formula") != formula
                or type(after.get("charge")) is not int or after.get("charge") != charge
                or after.get("annotation") != {"sbo": "SBO:0000247", "kegg.compound": kegg, "chebi": chebi,
                                               "metanetx.chemical": mnx, "bigg.metabolite": bigg}):
            raise ValueError("Thiamine phosphate chemistry spec values differ from the reviewed change")
    r1342 = spec["reactions"]["R1342"]
    if r1342.get("annotation_after") != {**r1342.get("annotation_before", {}), "ec-code": "2.7.4.16"}:
        raise ValueError("Thiamine phosphate chemistry spec values differ from the reviewed change")


def _state(current_notes, notes, row, before, after, label):
    """Return True when the change is present, False when absent; raise on anything else."""
    if not (before or after):
        raise ValueError(f"{label} precondition differs; no edits applied")
    current = _unescaped(current_notes)
    if any(key in current and current[key] != value for key, value in notes.items()):
        raise ValueError(f"Conflicting {label} thiamine phosphate chemistry notes; no edits applied")
    if {k: v for k, v in current.items() if k not in notes} != row["notes_before"]:
        raise ValueError(f"{label} notes outside this curation differ; no edits applied")
    if after:
        # The after-state without these notes came from somewhere else; its evidence would be unrecorded.
        if current != {**row["notes_before"], **notes}:
            raise ValueError(f"{label} already carries the change without this curation; no edits applied")
        return True
    if any(key in current for key in notes):
        raise ValueError(f"{label} carries this curation's notes without the change; no edits applied")
    return False


def apply_thiamine_phosphate_chemistry(model, spec=None):
    """Set the reviewed formulas, charges and IDs and R1342's EC; fail closed on any unexpected model state."""
    if spec is None:
        spec = json.loads(CURATION_PATH.read_text())
    _check_spec(spec)
    if "biomass_C" not in model.reactions or TPP not in {m.id for m in model.reactions.biomass_C.metabolites}:
        raise ValueError("The thiamine phosphate chemistry requires the TPP requirement in biomass_C; "
                         "no edits applied")
    for rid, stoichiometry in spec["context_reactions"].items():
        if rid not in model.reactions:
            raise ValueError(f"{rid} missing; no edits applied")
        reaction = model.reactions.get_by_id(rid)
        current = {m.id: c for m, c in reaction.metabolites.items()}
        if rid == "biomass_C":
            current = {TPP: current[TPP]}
        if current != stoichiometry:
            raise ValueError(f"{rid} differs from the reviewed reaction context; no edits applied")
        # The balance review assumed the partners' chemistry, e.g. neutral ATP; a migrated partner needs a new review.
        for metabolite in reaction.metabolites:
            if (rid != "biomass_C" and metabolite.id not in SPECIES_IDS
                    and [metabolite.formula, metabolite.charge] != spec["partner_species"].get(metabolite.id)):
                raise ValueError(f"{metabolite.id} on {rid} differs from the reviewed partner chemistry; "
                                 "no edits applied")
    done = set()
    for mid, row in spec["species"].items():
        if mid not in model.metabolites:
            raise ValueError(f"{mid} missing; no edits applied")
        metabolite = model.metabolites.get_by_id(mid)
        if (metabolite.name != row["name"] or metabolite.compartment != row["compartment"]
                or sorted(r.id for r in metabolite.reactions) != row["reactions_before"]):
            raise ValueError(f"{mid} identity or reaction context differs; no edits applied")
        current = (metabolite.formula, metabolite.charge, _annotation(metabolite.annotation))
        before = current == (row["before"]["formula"], row["before"]["charge"], _annotation(row["before"]["annotation"]))
        after = current == (row["after"]["formula"], row["after"]["charge"], _annotation(row["after"]["annotation"]))
        if _state(metabolite.notes, _notes(spec, row), row, before, after, mid):
            done.add(mid)
    row = spec["reactions"]["R1342"]
    reaction = model.reactions.get_by_id("R1342")
    if reaction.name != row["name"]:
        raise ValueError("R1342 precondition differs; no edits applied")
    annotation = _annotation(reaction.annotation)
    if _state(reaction.notes, _notes(spec, row), row, annotation == _annotation(row["annotation_before"]),
              annotation == _annotation(row["annotation_after"]), "R1342"):
        done.add("R1342")
    if done and done != SPECIES_IDS | REACTION_IDS:
        raise ValueError("Only part of the thiamine phosphate chemistry is present; no edits applied")
    if not done:
        for mid, row in spec["species"].items():
            metabolite = model.metabolites.get_by_id(mid)
            metabolite.formula = row["after"]["formula"]
            metabolite.charge = row["after"]["charge"]
            metabolite.annotation = copy.deepcopy(row["after"]["annotation"])
            metabolite.notes = {**metabolite.notes, **_notes(spec, row)}
        row = spec["reactions"]["R1342"]
        reaction.annotation = copy.deepcopy(row["annotation_after"])
        reaction.notes = {**reaction.notes, **_notes(spec, row)}
    return {"item": "thiamine phosphate chemistry", "curation_id": spec["curation_id"],
            "status": "already_correct" if done else "applied",
            "changed_species": [] if done else sorted(spec["species"]),
            "changed_reactions": [] if done else sorted(REACTION_IDS), "evidence_status": spec["status"]}
