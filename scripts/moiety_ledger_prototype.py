#!/usr/bin/env python3
"""A mass-balanced proof of concept for aggregate acyl-moiety accounting.

This module deliberately builds a tiny *in-memory* COBRA model.  It does not
read or write ``model.xml`` and is not part of the production rebuild.

The prototype represents phosphatidate as one generic glycerol-phosphate
backbone plus two separately conserved acyl-residue ledger metabolites.  That
avoids making every PA(ChainA/ChainB) molecule an internal metabolite while
retaining the total amount of each acyl chain.  It intentionally cannot retain
sn-1/sn-2 pairing: the two positional routes to PA(16:0/18:1) collapse to the
same aggregate output.  The accompanying tests make both properties explicit.
"""

from __future__ import annotations

from dataclasses import dataclass

from cobra import Metabolite, Model, Reaction


@dataclass(frozen=True)
class Chain:
    """One acyl-CoA and the residue transferred when CoA leaves."""

    id: str
    acyl_coa_formula: str
    residue_formula: str


# Two chains are enough to exercise both positional paths without concealing
# the aggregate-versus-molecular-species distinction under a large network.
CHAINS = (
    Chain("C16_0", "C37H66N7O17P3S", "C16H30O"),
    Chain("C18_1", "C39H68N7O17P3S", "C18H32O"),
)

COA_FORMULA = "C21H36N7O16P3S"
G3P_FORMULA = "C3H9O6P"
PA_16_0_18_1_FORMULA = "C37H71O8P"


def _metabolite(identifier: str, name: str, formula: str) -> Metabolite:
    """Create a neutral demonstration metabolite with a definite formula."""

    return Metabolite(identifier, name=name, formula=formula, charge=0, compartment="c")


def _reaction(identifier: str, name: str, stoichiometry: dict[Metabolite, float]) -> Reaction:
    reaction = Reaction(identifier, name=name)
    reaction.bounds = (0.0, 1000.0)
    reaction.annotation["moiety_ledger_internal"] = True
    reaction.add_metabolites(stoichiometry)
    return reaction


def _boundary(
    identifier: str, metabolite: Metabolite, coefficient: float, upper_bound: float = 1.0
) -> Reaction:
    """Create a source or demand used only to exercise the in-memory model."""

    reaction = Reaction(identifier, name=identifier)
    reaction.bounds = (0.0, upper_bound)
    reaction.annotation["moiety_ledger_boundary"] = True
    reaction.add_metabolites({metabolite: coefficient})
    return reaction


def build_pa_moiety_ledger_demo() -> Model:
    """Build a two-chain PA ledger model with no generic chemical formulas.

    ``lpa_backbone`` and ``pa_backbone`` are class backbones.  Their acyl
    content travels in ``*_moiety_C16_0`` and ``*_moiety_C18_1`` metabolites.
    The model can make aggregate PA(16:0/18:1) through either positional route,
    but cannot say which residue occupies sn-1 or sn-2.
    """

    model = Model("pa_moiety_ledger_demo")
    model.solver = "glpk"

    coa = _metabolite("coa", "coenzyme A", COA_FORMULA)
    g3p = _metabolite("g3p", "glycerol 3-phosphate", G3P_FORMULA)
    lpa_backbone = _metabolite("lpa_backbone", "LPA glycerol-phosphate backbone", G3P_FORMULA)
    pa_backbone = _metabolite("pa_backbone", "PA glycerol-phosphate backbone", G3P_FORMULA)
    pa_16_0_18_1 = _metabolite(
        "pa_16_0_18_1",
        "phosphatidate (16:0/18:1 aggregate output)",
        PA_16_0_18_1_FORMULA,
    )

    acyl_coas: dict[str, Metabolite] = {}
    lpa_moieties: dict[str, Metabolite] = {}
    pa_moieties: dict[str, Metabolite] = {}
    for chain in CHAINS:
        acyl_coas[chain.id] = _metabolite(
            f"acyl_coa_{chain.id}", f"acyl-CoA {chain.id}", chain.acyl_coa_formula
        )
        lpa_moieties[chain.id] = _metabolite(
            f"lpa_moiety_{chain.id}", f"LPA acyl residue {chain.id}", chain.residue_formula
        )
        pa_moieties[chain.id] = _metabolite(
            f"pa_moiety_{chain.id}", f"PA acyl residue {chain.id}", chain.residue_formula
        )

    model.add_metabolites(
        [coa, g3p, lpa_backbone, pa_backbone, pa_16_0_18_1]
        + list(acyl_coas.values())
        + list(lpa_moieties.values())
        + list(pa_moieties.values())
    )

    reactions: list[Reaction] = []
    for chain in CHAINS:
        chain_id = chain.id
        # First acylation: the lipid class and its transferred acyl residue are
        # represented separately, but the combined equation is element-balanced.
        reactions.append(
            _reaction(
                f"LPA_{chain_id}",
                f"first acylation with {chain_id}",
                {
                    acyl_coas[chain_id]: -1,
                    g3p: -1,
                    coa: 1,
                    lpa_backbone: 1,
                    lpa_moieties[chain_id]: 1,
                },
            )
        )
        # Second acylation creates the PA backbone and one newly attached
        # residue.  The pre-existing residue is carried by the transfer below.
        reactions.append(
            _reaction(
                f"PA_ADD_{chain_id}",
                f"second acylation with {chain_id}",
                {
                    lpa_backbone: -1,
                    acyl_coas[chain_id]: -1,
                    coa: 1,
                    pa_backbone: 1,
                    pa_moieties[chain_id]: 1,
                },
            )
        )
        reactions.append(
            _reaction(
                f"LPA_TO_PA_{chain_id}",
                f"carry {chain_id} residue from LPA to PA ledger",
                {lpa_moieties[chain_id]: -1, pa_moieties[chain_id]: 1},
            )
        )

    reactions.append(
        _reaction(
            "ASSEMBLE_PA_C16_0_C18_1",
            "assemble aggregate PA output from backbone and two moieties",
            {
                pa_backbone: -1,
                pa_moieties["C16_0"]: -1,
                pa_moieties["C18_1"]: -1,
                pa_16_0_18_1: 1,
            },
        )
    )

    # Bounded sources make a maximum output of exactly one PA molecule.  The
    # demand reaction is the objective; none of these boundary reactions claim
    # to be chemistry in the real GEM.
    reactions.extend(
        [
            _boundary("SRC_G3P", g3p, 1),
            _boundary("SRC_ACYL_C16_0", acyl_coas["C16_0"], 1),
            _boundary("SRC_ACYL_C18_1", acyl_coas["C18_1"], 1),
            # Two acyl-transfer reactions release two CoA per PA output.
            _boundary("DM_COA", coa, -1, upper_bound=2.0),
            _boundary("DM_PA_C16_0_C18_1", pa_16_0_18_1, -1),
        ]
    )
    model.add_reactions(reactions)
    model.objective = "DM_PA_C16_0_C18_1"
    return model


def internal_mass_balance_errors(model: Model) -> dict[str, dict[str, float]]:
    """Return formula/charge residuals for the ledger's biochemical reactions."""

    return {
        reaction.id: reaction.check_mass_balance()
        for reaction in model.reactions
        if reaction.annotation.get("moiety_ledger_internal") and reaction.check_mass_balance()
    }


def solve_demo() -> dict[str, float | int]:
    """Run the proof of concept and return compact, machine-checkable results."""

    model = build_pa_moiety_ledger_demo()
    solution = model.optimize()
    if solution.status != "optimal":
        raise RuntimeError(f"ledger demo did not solve: {solution.status}")

    results: dict[str, float | int] = {
        "objective": float(solution.objective_value),
        "mass_balance_error_count": len(internal_mass_balance_errors(model)),
    }
    for first_chain in CHAINS:
        with model:
            model.reactions.get_by_id(f"LPA_{first_chain.id}").knock_out()
            positional_solution = model.optimize()
            results[f"without_lpa_{first_chain.id}"] = float(positional_solution.objective_value)

    with model:
        model.reactions.get_by_id("SRC_ACYL_C16_0").knock_out()
        closed_supply = model.optimize()
        results["without_acyl_c16_0_supply"] = float(closed_supply.objective_value)
    return results


if __name__ == "__main__":
    for key, value in solve_demo().items():
        print(f"{key}: {value}")
