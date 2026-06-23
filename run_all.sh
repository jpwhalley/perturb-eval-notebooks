#!/usr/bin/env bash
# run_all.sh — Execute perturbation-evaluation notebooks in dependency order.
#
# Usage:
#   ./run_all.sh                  # print this help (does not start heavy compute)
#   ./run_all.sh --figures        # figure/table notebooks from precomputed CSVs (~2 min)
#   ./run_all.sh --demo           # 2-seed pipeline run to verify installation
#   ./run_all.sh --analysis       # recompute eval CSVs from existing per-seed predictions
#   ./run_all.sh --train-full     # full n=100 CPA + matched-n GEARS (~13 hours)
#
# Notebooks not listed in NOTEBOOKS_READY below are skipped (placeholders).
# Add a notebook's stem to NOTEBOOKS_READY once it has been implemented and
# tested.

set -euo pipefail

LOGDIR="logs"
mkdir -p "$LOGDIR"

# ── Notebook lists (dependency order) ───────────────────────────────────────

DATA_NOTEBOOKS=(
    D01_replogle_atlases
    D02_severity_references
    D03_holdout_specifications
)

ANALYSIS_NOTEBOOKS=(
    P01_train_cpa
    P02_train_gears
    P03_evaluate_severity
    P04_loo_and_winsorisation
)

FIGURE_NOTEBOOKS=(
    01_figure1
    02_table1_mechanism_summary
    03_table2_gears_scope
    04_appendix_threshold_sensitivity
)

# Notebooks that have been implemented and are ready to execute.
# Notebooks not in this list are skipped with a pending message; this prevents
# the scaffold from appearing runnable before the implementations are in place.
NOTEBOOKS_READY=(
    # Figure / table notebooks — read precomputed CSVs, render manuscript outputs.
    01_figure1
    02_table1_mechanism_summary
    03_table2_gears_scope
    04_appendix_threshold_sensitivity

    # Analysis notebooks — read per-seed prediction h5ads, produce eval CSVs.
    # Both gracefully fall through with a helpful message when predictions are absent.
    P03_evaluate_severity
    P04_loo_and_winsorisation

    # D-series — data acquisition. Each gracefully falls through with a helpful
    # message when its prerequisites are absent.
    D01_replogle_atlases
    D02_severity_references
    D03_holdout_specifications

    # Training wrappers — orchestration shells. The substantive training loops
    # delegate to src/perturb_eval/train_{cpa,gears}.py (currently documented
    # stubs). Each notebook gracefully exits with a usage message when neither
    # QUICK_DEMO nor FULL_TRAINING is set, and with a prerequisite-check message
    # when the D-series outputs are absent.
    P01_train_cpa
    P02_train_gears
)

# ── Parse arguments ─────────────────────────────────────────────────────────

print_help() {
    sed -n '/^# Usage:/,/^# Notebooks/p' "$0" | sed 's/^# *//' | sed '$d'
}

MODE="${1:-help}"

case "$MODE" in
    --figures)
        NOTEBOOKS=("${FIGURE_NOTEBOOKS[@]}")
        ;;
    --demo)
        NOTEBOOKS=("${DATA_NOTEBOOKS[@]}" "${ANALYSIS_NOTEBOOKS[@]}" "${FIGURE_NOTEBOOKS[@]}")
        export QUICK_DEMO=1
        ;;
    --analysis)
        NOTEBOOKS=("${ANALYSIS_NOTEBOOKS[@]}" "${FIGURE_NOTEBOOKS[@]}")
        ;;
    --train-full)
        NOTEBOOKS=("${DATA_NOTEBOOKS[@]}" "${ANALYSIS_NOTEBOOKS[@]}" "${FIGURE_NOTEBOOKS[@]}")
        export FULL_TRAINING=1
        ;;
    -h|--help|help|"")
        print_help
        exit 0
        ;;
    *)
        echo "Unknown mode: $MODE"
        echo
        print_help
        exit 1
        ;;
esac

# ── Helper: is the notebook marked ready? ───────────────────────────────────

is_ready() {
    local nb="$1"
    # Guard against unbound-variable error under `set -u` when the array is empty.
    if [[ ${#NOTEBOOKS_READY[@]} -eq 0 ]]; then
        return 1
    fi
    for ready in "${NOTEBOOKS_READY[@]}"; do
        if [[ "$ready" == "$nb" ]]; then
            return 0
        fi
    done
    return 1
}

# ── Execute ─────────────────────────────────────────────────────────────────

TOTAL=${#NOTEBOOKS[@]}
EXECUTED=0
SKIPPED=0
START_ALL=$(date +%s)

echo "=== Mode: $MODE (${TOTAL} notebook(s) in dependency chain) ==="
echo

for nb in "${NOTEBOOKS[@]}"; do
    NB_FILE="${nb}.ipynb"
    LOG_FILE="${LOGDIR}/${nb}.log"

    if [[ ! -f "$NB_FILE" ]]; then
        echo "  SKIP  $NB_FILE (file not found)"
        SKIPPED=$((SKIPPED + 1))
        continue
    fi

    if ! is_ready "$nb"; then
        echo "  SKIP  $NB_FILE (marked pending; not yet in NOTEBOOKS_READY)"
        SKIPPED=$((SKIPPED + 1))
        continue
    fi

    EXECUTED=$((EXECUTED + 1))
    echo -n "  [$EXECUTED] $NB_FILE ... "
    START=$(date +%s)

    if uv run jupyter nbconvert \
        --to notebook \
        --execute \
        --inplace \
        --ExecutePreprocessor.timeout=7200 \
        --ExecutePreprocessor.startup_timeout=120 \
        "$NB_FILE" \
        > "$LOG_FILE" 2>&1; then
        ELAPSED=$(( $(date +%s) - START ))
        echo "OK (${ELAPSED}s)"
    else
        ELAPSED=$(( $(date +%s) - START ))
        echo "FAILED (${ELAPSED}s) — see $LOG_FILE"
        echo
        echo "=== STOPPED: $NB_FILE failed ==="
        tail -20 "$LOG_FILE"
        exit 1
    fi
done

TOTAL_TIME=$(( $(date +%s) - START_ALL ))
echo
echo "=== Done: $EXECUTED executed, $SKIPPED skipped, ${TOTAL_TIME}s total ==="
