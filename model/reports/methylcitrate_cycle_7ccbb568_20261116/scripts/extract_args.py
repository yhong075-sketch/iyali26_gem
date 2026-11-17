"""Write the 7ccbb568 build argv (without the script path and --output-model) to args_7ccbb568.txt."""
import json
import sys
from pathlib import Path

build = json.loads(Path(sys.argv[1]).read_text())
argv = build["argv"][1:]
argv = argv[:argv.index("--output-model")]
Path(sys.argv[2]).write_text("\n".join(argv) + "\n")
print("head", build["head"], "args", len(argv))
