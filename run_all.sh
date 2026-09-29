#!/usr/bin/env bash
#
# Re-run every experiment of this repository in one go, on the current machine.
#
#   ./run_all.sh                 # instantiate, stamp the environment, run all six + summary
#   ./run_all.sh --quick         # EXP{4,6}_QUICK=1: smoke test, writes *_quick files
#   ./run_all.sh --only 2,6      # re-run selected experiments only
#   ./run_all.sh --help
#
# Everything is written by the Julia scripts themselves into results/ (and is
# overwritten in place); this script only fixes the order, records the machine,
# keeps a log per experiment, and reports which ones succeeded.
#
# The experiments are run strictly one after another: exp2 measures nanosecond-
# to-microsecond timings, which are only meaningful on an otherwise idle machine.

set -euo pipefail

ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
JULIA="${JULIA:-julia}"

QUICK=0
KEEP_GOING=0
INSTANTIATE=1
DRY_RUN=0
ONLY=""

usage() {
  cat <<'EOF'
Usage: ./run_all.sh [options]

  -q, --quick           Set EXP4_QUICK=1 and EXP6_QUICK=1. Experiments 4 and 6
                        then run a tiny case list and write *_quick.json /
                        summary_*_quick.md, leaving their committed results
                        alone. Experiments 1-3 have no quick mode: they run in
                        full and DO overwrite their committed results (including
                        exp2's timings, measured under whatever load the machine
                        is under). Add --only 4,6 to touch nothing else.
  -o, --only LIST       Comma- or space-separated experiment numbers to run,
                        e.g. --only 2,6. make_summary.jl runs only when 1, 2
                        and 3 all ran successfully.
  -k, --keep-going      Run the remaining experiments after a failure instead of
                        stopping at the first one.
      --no-instantiate  Skip Pkg.instantiate() (assume the environment is ready).
  -n, --dry-run         Print what would run, run nothing.
  -h, --help            This message.

Environment:
  JULIA                 Julia binary to use (default: julia).
  JULIA_NUM_THREADS     Passed through; recorded in results/environment.md.

Exit status is nonzero if any experiment failed. Commit results/ only when all
of them are green: a partial re-run leaves the directory half old, half new.
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    -q|--quick)        QUICK=1 ;;
    -k|--keep-going)   KEEP_GOING=1 ;;
    -n|--dry-run)      DRY_RUN=1 ;;
    --no-instantiate)  INSTANTIATE=0 ;;
    -o|--only)         ONLY="${2:-}"; shift ;;
    --only=*)          ONLY="${1#*=}" ;;
    -h|--help)         usage; exit 0 ;;
    *) printf 'run_all.sh: unknown option %s\n\n' "$1" >&2; usage >&2; exit 2 ;;
  esac
  shift
done

# Experiment number -> script. Order is the run order; make_summary.jl consumes
# the tables of 1-3 and so must come last.
EXP_NUMS=(1 2 3 4 5 6)
EXP_FILES=(
  exp1_exactness.jl
  exp2_sparsity.jl
  exp3_conformity.jl
  exp4_convergence.jl
  exp5_interpolation.jl
  exp6_curlcurl.jl
)

selected() {
  [[ -z "$ONLY" ]] && return 0
  local n="$1" want
  for want in ${ONLY//,/ }; do
    [[ "$want" == "$n" ]] && return 0
  done
  return 1
}

if [[ "$DRY_RUN" == 0 ]] && ! command -v "$JULIA" >/dev/null 2>&1; then
  printf 'run_all.sh: %s not found on PATH (set JULIA=/path/to/julia)\n' "$JULIA" >&2
  exit 127
fi

STAMP="$(date +%Y%m%d-%H%M%S)"
LOGDIR="$ROOT/logs/run-$STAMP"
RESULTS="$ROOT/results"

printf '=== bmm-2026-experiments: re-run %s ===\n' "$STAMP"
printf 'repo   : %s\n' "$ROOT"
printf 'julia  : %s (%s)\n' "$JULIA" "$("$JULIA" --version 2>/dev/null || echo 'version unknown')"
printf 'logs   : %s\n' "$LOGDIR"
if [[ "$QUICK" == 1 ]]; then
  printf 'mode   : QUICK (exp4/exp6 smoke test, *_quick outputs)\n'
  if [[ -z "$ONLY" ]] || selected 1 || selected 2 || selected 3; then
    printf '         experiments 1-3 have no quick mode: they run in full and\n'
    printf '         overwrite their committed results (--only 4,6 to avoid that)\n'
  fi
fi
[[ -n "$ONLY" ]] && printf 'only   : %s\n' "$ONLY"
printf '\n'

if [[ "$DRY_RUN" == 1 ]]; then
  [[ "$INSTANTIATE" == 1 ]] && printf 'would run: Pkg.instantiate()\n'
  printf 'would run: %s --project=. env_stamp.jl\n' "$JULIA"
  for i in "${!EXP_NUMS[@]}"; do
    selected "${EXP_NUMS[$i]}" && printf 'would run: %s --project=. %s\n' "$JULIA" "${EXP_FILES[$i]}"
  done
  if selected 1 && selected 2 && selected 3; then
    printf 'would run: %s --project=. make_summary.jl\n' "$JULIA"
  else
    printf 'would skip: make_summary.jl (experiments 1-3 are not all selected)\n'
  fi
  exit 0
fi

mkdir -p "$LOGDIR" "$RESULTS"

if [[ "$QUICK" == 1 ]]; then
  export EXP4_QUICK=1 EXP6_QUICK=1
fi

# ── Environment ───────────────────────────────────────────────────────────────
# Instantiating under a different Julia than the one recorded in Manifest.toml
# rewrites the manifest; that is fine, but the regenerated manifest has to be
# committed with the results for the recorded environment to stay truthful.
MANIFEST_BEFORE="$(git -C "$ROOT" status --porcelain -- Manifest.toml 2>/dev/null || true)"

# Parse every script before anything is overwritten: a syntax error in a late
# script would otherwise surface hours in, on top of half-regenerated results.
# Meta.parseall reports syntax errors as :error nodes inside the returned
# expression rather than by throwing, so the tree has to be walked.
printf -- '--- parse check ---\n'
if ! ( cd "$ROOT" && "$JULIA" --project=. -e '
  function parse_errors(ex, acc = String[])
    ex isa Expr || return acc
    if ex.head === :error || ex.head === :incomplete
      push!(acc, sprint(showerror, ex.args[1]))
    else
      for a in ex.args; parse_errors(a, acc); end
    end
    acc
  end
  bad = 0
  for f in ARGS
    errs = parse_errors(Meta.parseall(read(f, String); filename = f))
    isempty(errs) && continue
    global bad += 1
    println(stderr, "syntax error in ", f, ":\n", first(errs))
  end
  bad == 0 ? println("all scripts parse") : exit(1)' \
    env_stamp.jl exp1_exactness.jl exp2_sparsity.jl exp3_conformity.jl \
    exp4_convergence.jl exp5_interpolation.jl exp6_curlcurl.jl make_summary.jl \
  ) 2>&1 | tee "$LOGDIR/parse.log"; then
  printf '\nParse check failed: nothing was run and results/ is untouched.\n' >&2
  exit 1
fi
printf '\n'

if [[ "$INSTANTIATE" == 1 ]]; then
  printf -- '--- instantiating the environment ---\n'
  "$JULIA" --project="$ROOT" -e 'using Pkg; Pkg.instantiate()' 2>&1 | tee "$LOGDIR/instantiate.log"
  printf '\n'
fi

MANIFEST_AFTER="$(git -C "$ROOT" status --porcelain -- Manifest.toml 2>/dev/null || true)"
MANIFEST_CHANGED=0
if [[ "$MANIFEST_BEFORE" != "$MANIFEST_AFTER" ]]; then
  MANIFEST_CHANGED=1
  printf '!!! Manifest.toml changed while instantiating (most likely a different\n'
  printf '!!! Julia version than the one it records). Review and commit it with\n'
  printf '!!! the new results:  git diff Manifest.toml\n\n'
fi

printf -- '--- recording the machine ---\n'
"$JULIA" --project="$ROOT" "$ROOT/env_stamp.jl" ${ONLY:+"$ONLY"} 2>&1 | tee "$LOGDIR/env_stamp.log"
printf '\n'

# ── Experiments ───────────────────────────────────────────────────────────────
declare -a RAN_NUMS=() RAN_STATUS=() RAN_SECONDS=()
FAILURES=0

run_script() {
  local num="$1" script="$2" log="$LOGDIR/${2%.jl}.log" start elapsed status=0

  printf -- '--- [%s] %s ---\n' "$num" "$script"
  [[ "$script" == exp2_sparsity.jl ]] &&
    printf '    (microbenchmarks: leave the machine idle until this finishes)\n'

  start=$SECONDS
  if ( cd "$ROOT" && "$JULIA" --project=. "$script" ) 2>&1 | tee "$log"; then
    status=0
  else
    status="${PIPESTATUS[0]}"
  fi
  elapsed=$((SECONDS - start))

  RAN_NUMS+=("$num")
  RAN_SECONDS+=("$elapsed")
  if [[ "$status" == 0 ]]; then
    RAN_STATUS+=("ok")
    printf '    done in %ds\n\n' "$elapsed"
  else
    RAN_STATUS+=("FAILED (exit $status)")
    FAILURES=$((FAILURES + 1))
    printf '    FAILED after %ds (exit %s), see %s\n\n' "$elapsed" "$status" "$log"
    [[ "$KEEP_GOING" == 1 ]] || return 1
  fi
  return 0
}

TOTAL_START=$SECONDS
RUN_ABORTED=0

for i in "${!EXP_NUMS[@]}"; do
  selected "${EXP_NUMS[$i]}" || continue
  if ! run_script "${EXP_NUMS[$i]}" "${EXP_FILES[$i]}"; then
    RUN_ABORTED=1
    break
  fi
done

# make_summary.jl concatenates the tables of experiments 1-3; running it on a
# mixture of fresh and stale tables would produce exactly the inconsistency a
# re-run is meant to remove.
summary_inputs_fresh() {
  local n
  for n in 1 2 3; do
    selected "$n" || return 1
    local k
    for k in "${!RAN_NUMS[@]}"; do
      [[ "${RAN_NUMS[$k]}" == "$n" && "${RAN_STATUS[$k]}" == "ok" ]] && continue 2
    done
    return 1
  done
  return 0
}

if [[ "$RUN_ABORTED" == 0 ]] && summary_inputs_fresh; then
  run_script summary make_summary.jl || RUN_ABORTED=1
else
  printf -- '--- skipping make_summary.jl (experiments 1-3 were not all re-run) ---\n\n'
fi

TOTAL=$((SECONDS - TOTAL_START))

# ── Report ────────────────────────────────────────────────────────────────────
printf '=== summary ===\n'
printf '%-10s %-8s %s\n' "target" "time" "status"
for i in "${!RAN_NUMS[@]}"; do
  printf '%-10s %-8s %s\n' "${RAN_NUMS[$i]}" "${RAN_SECONDS[$i]}s" "${RAN_STATUS[$i]}"
done
printf 'total: %dm%ds\n' $((TOTAL / 60)) $((TOTAL % 60))

{
  printf '\n## Run log (%s)\n\n' "$STAMP"
  printf 'Produced by `run_all.sh`'
  [[ "$QUICK" == 1 ]] && printf ' with `--quick` (experiments 4 and 6 wrote `*_quick` files)'
  [[ -n "$ONLY" ]] && printf ', restricted to experiments `%s`' "$ONLY"
  printf '.\n\n| target | wall time | status |\n|---|---|---|\n'
  for i in "${!RAN_NUMS[@]}"; do
    printf '| %s | %ss | %s |\n' "${RAN_NUMS[$i]}" "${RAN_SECONDS[$i]}" "${RAN_STATUS[$i]}"
  done
  printf '\nTotal wall time: %dm%ds. Logs: `logs/run-%s/`.\n' \
    $((TOTAL / 60)) $((TOTAL % 60)) "$STAMP"
  [[ "$MANIFEST_CHANGED" == 1 ]] &&
    printf '\n`Manifest.toml` was rewritten by `Pkg.instantiate()` during this run.\n'
} >> "$RESULTS/environment.md"

if [[ "$FAILURES" -gt 0 ]]; then
  printf '\n%d target(s) failed: results/ now mixes new and old data. Fix and re-run\n' "$FAILURES"
  printf 'before committing (git checkout -- results/ restores the committed state).\n'
  exit 1
fi

if [[ "$QUICK" == 1 ]]; then
  printf '\nQuick mode: experiments 4 and 6 wrote *_quick files only.\n'
  overwrote_committed=0
  for i in "${!RAN_NUMS[@]}"; do
    case "${RAN_NUMS[$i]}" in 1|2|3|summary) overwrote_committed=1 ;; esac
  done
  if [[ "$overwrote_committed" == 1 ]]; then
    printf 'The experiments 1-3 that ran have no quick mode: their results were\n'
    printf 'regenerated in place.\n'
  else
    printf 'No committed table was overwritten (environment.md aside).\n'
  fi
  printf -- '--quick is a smoke test; re-run without it for publishable data.\n'
else
  printf '\nAll targets green. Review and commit results/ together with Manifest.toml.\n'
fi

# Explicit: the exit status of a green run must not depend on the last printf.
exit 0
