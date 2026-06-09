#!/bin/bash
# Convenience wrapper for the mouse SCN2A/SCN8A ISM pipeline.
#
# Equivalent to:
#   bash submit_all_ism_scn2a.sh --species mouse "$@"

set -eo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec bash "${SCRIPT_DIR}/submit_all_ism_scn2a.sh" --species mouse "$@"
