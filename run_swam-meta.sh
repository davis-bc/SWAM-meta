#!/usr/bin/env bash
# Run SWAM-meta on Slurm with the default resilient-style profile and emit a
# concise post-run status summary.
# Usage: bash run_swam-meta.sh [input_dir] [output_dir] [extra snakemake args...]

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

INPUT="${1:-$(pwd)/input}"
OUTPUT="${2:-$(pwd)/output}"
shift $(( $# >= 1 ? 1 : 0 ))
shift $(( $# >= 1 ? 1 : 0 ))

mkdir -p "$OUTPUT/logs/run_status"
TIMESTAMP="$(date +%Y%m%d-%H%M%S)"
DRIVER_LOG="$OUTPUT/logs/run_status/snakemake.$TIMESTAMP.log"

echo "Input:      $INPUT"
echo "Output:     $OUTPUT"
echo "Profile:    config/slurm/"
echo "Driver log: $DRIVER_LOG"
echo ""

snakemake \
    --profile config/slurm/ \
    --config in_dir="$INPUT" out_dir="$OUTPUT" \
    --local-cores 1 \
    --quiet rules \
    "$@" 2>&1 | tee "$DRIVER_LOG"

SNAKEMAKE_EXIT=${PIPESTATUS[0]}

python3 workflow/scripts/build_run_report.py \
    --input-dir "$INPUT" \
    --output-dir "$OUTPUT" \
    --driver-log "$DRIVER_LOG"

REPORT_EXIT=$?

echo ""
echo "Run report:        $OUTPUT/logs/run_status/run_report.txt"
echo "Rule status table: $OUTPUT/logs/run_status/sample_rule_status.tsv"

if [ "$SNAKEMAKE_EXIT" -ne 0 ]; then
    exit "$SNAKEMAKE_EXIT"
fi

exit "$REPORT_EXIT"
