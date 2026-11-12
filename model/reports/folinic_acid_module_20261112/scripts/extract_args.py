"""Write the e584151f build argv (without the script path and --output-model) to args_e584151f.txt."""
import json, sys
from pathlib import Path
build = json.loads(Path(sys.argv[1]).read_text())
argv = build["argv"][1:]
argv = argv[:argv.index("--output-model")]
Path(sys.argv[2]).write_text("\n".join(argv) + "\n")
print("head", build["head"], "dirty", str(build["dirty_status"])[:300])
