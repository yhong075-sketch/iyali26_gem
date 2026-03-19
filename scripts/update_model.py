import logging
from pathlib import Path
from cobra.io import read_sbml_model, write_sbml_model

# Setup clean logging
logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

# Dynamically resolve paths so this script works anywhere
REPO_ROOT = Path(__file__).resolve().parent.parent
STARTING_MODEL_PATH = REPO_ROOT / "data" / "iyli21.xml"
OUTPUT_MODEL_PATH = REPO_ROOT / "model.xml"


def main():
    if not STARTING_MODEL_PATH.exists():
        logger.error(f"Could not find starting model at {STARTING_MODEL_PATH}")
        return

    logger.info(f"Loading raw model: {STARTING_MODEL_PATH.name}")
    model = read_sbml_model(str(STARTING_MODEL_PATH))

    # Model update functions live here
    # model = apply_compartment_fixes(model)
    # model = apply_gene_rules(model)
    # model = add_new_reactions(model)

    logger.info(f"Saving updated model to: {OUTPUT_MODEL_PATH.name}")
    write_sbml_model(model, str(OUTPUT_MODEL_PATH))
    logger.info("Model build complete.")


if __name__ == "__main__":
    main()
