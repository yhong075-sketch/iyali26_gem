"""Build the explicit, optional dipeptide chemical metadata candidate."""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts.gem_annotate.dipeptide_chemistry import SPEC_PATH, build_candidate_file
from scripts.gem_annotate.execution import execution_limits


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', required=True)
    parser.add_argument('--source-sha256', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--patch', default=str(SPEC_PATH))
    parser.add_argument('--enabled', action='store_true')
    args = parser.parse_args()
    with execution_limits(no_solve=True, allow_network=False) as counts:
        result = build_candidate_file(args.source, args.source_sha256, args.output, args.enabled, args.patch)
    print(result['candidate'], result['output_sha256'], 'reload and unchanged LP verified', counts)


if __name__ == '__main__':
    main()
