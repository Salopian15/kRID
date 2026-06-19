#!/bin/bash
# ============================================================================
# aggregate.sh  (v6-aware)
#
# Aggregate per-query summary.tsv files into a cohort-level view.
#
# Usage:
#   ./aggregate.sh <project_workspace>
#
#
# Env vars (read from environment, defaults match the pipeline):
#   JACCARD_THRESHOLD       0.85
#   CONTAINMENT_THRESHOLD   0.95
# ============================================================================

# Note to future jack check for self v self runs
set -euo pipefail

# Add exclude self, was causing issues with self vs self checks of the EIAR database.
EXCLUDE_SELF=0
POSITIONAL=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --exclude-self) EXCLUDE_SELF=1; shift ;;
        -h|--help) echo "Usage: $0 [--exclude-self] <project_workspace>"; exit 0 ;;
        --*) echo "Unknown option: $1" >&2; exit 1 ;;
        *) POSITIONAL+=("$1"); shift ;;
    esac
done
set -- "${POSITIONAL[@]}"

PROJ_DIR="${1:-}"
[[ -n "${PROJ_DIR}" && -d "${PROJ_DIR}" ]] \
    || { echo "Usage: $0 [--exclude-self] <project_workspace>" >&2; exit 1; }

WORKDIR="${PROJ_DIR}/query_runs"
PANEL_UNION_DIR="${PROJ_DIR}/panel_union"
EXCLUDED="${PANEL_UNION_DIR}/panel_excluded.tsv"
SIZES="${PANEL_UNION_DIR}/accession_sizes.tsv"
SIZES_VAR="${PANEL_UNION_DIR}/accession_sizes_var.tsv"
U_FILE="${PANEL_UNION_DIR}/panel_union_size.txt"

JACCARD_THRESHOLD="${JACCARD_THRESHOLD:-0.85}"
CONTAINMENT_THRESHOLD="${CONTAINMENT_THRESHOLD:-0.95}"

COHORT_TOP="${PROJ_DIR}/cohort_top_hits.tsv"
COHORT_ALL="${PROJ_DIR}/cohort_all_hits.tsv"
COHORT_FREQ="${PROJ_DIR}/cohort_panel_frequency.tsv"
COHORT_REPORT="${PROJ_DIR}/cohort_report.txt"

COHORT_TOP_VAR="${PROJ_DIR}/cohort_top_hits_var.tsv"
COHORT_ALL_VAR="${PROJ_DIR}/cohort_all_hits_var.tsv"
COHORT_FREQ_VAR="${PROJ_DIR}/cohort_panel_frequency_var.tsv"

# Fallback header used only if no summary.tsv exists to get from.
DEFAULT_HEADER="sample	rank	panel	na	nb	n_intersect	n_union	jaccard	containment	shared_state_rate	containment_q	na_var	nb_var	n_intersect_var	jaccard_var	containment_var	flag	flag_var"

# Resolve column indices by header name .

COL_sample=""; COL_rank=""; COL_panel=""; COL_jaccard=""; COL_containment=""
COL_ssr=""; COL_containment_q=""; COL_na_var=""; COL_nb_var=""; COL_ni_var=""
COL_jaccard_var=""; COL_containment_var=""; COL_flag=""; COL_flag_var=""

resolve_cols() {
    local hdr="$1"
    local i=0 name
    while IFS= read -r name; do
        i=$((i+1))
        case "$name" in
            sample)            COL_sample=$i ;;
            rank)              COL_rank=$i ;;
            panel)             COL_panel=$i ;;
            jaccard)           COL_jaccard=$i ;;
            containment)       COL_containment=$i ;;
            shared_state_rate) COL_ssr=$i ;;
            containment_q)     COL_containment_q=$i ;;
            na_var)            COL_na_var=$i ;;
            nb_var)            COL_nb_var=$i ;;
            n_intersect_var)   COL_ni_var=$i ;;
            jaccard_var)       COL_jaccard_var=$i ;;
            containment_var)   COL_containment_var=$i ;;
            flag)              COL_flag=$i ;;
            flag_var)          COL_flag_var=$i ;;
        esac
    done < <(printf '%s\n' "$hdr" | tr '\t' '\n')
}

# Get the header from the first available summary.tsv.
SUMMARY_HEADER="${DEFAULT_HEADER}"
first_summary="$(find "${WORKDIR}" -mindepth 2 -maxdepth 2 -name summary.tsv 2>/dev/null | head -1 || true)"
if [[ -n "${first_summary}" && -e "${first_summary}" ]]; then
    SUMMARY_HEADER="$(head -1 "${first_summary}")"
fi
resolve_cols "${SUMMARY_HEADER}"

# Sanity: flag and jaccard columns are mandatory for the report to make sense (learn from past jacks mistakes).
[[ -n "${COL_flag}" ]]    || { echo "ERROR: could not locate 'flag' column in summary header." >&2; exit 1; }
[[ -n "${COL_jaccard}" ]] || { echo "ERROR: could not locate 'jaccard' column in summary header." >&2; exit 1; }
[[ -n "${COL_panel}" ]]   || COL_panel=3
[[ -n "${COL_sample}" ]]  || COL_sample=1
[[ -n "${COL_containment}" ]] || COL_containment=9
[[ -n "${COL_ssr}" ]]     || COL_ssr=10

# Check if variable mask is being used here. This is the part that (ideally) uses the more meaningful segment of the genome.
HAS_VAR=0
if [[ -n "${COL_flag_var}" && -n "${COL_jaccard_var}" ]]; then
    if find "${WORKDIR}" -mindepth 2 -maxdepth 2 -name summary_var.tsv 2>/dev/null | grep -q .; then
        HAS_VAR=1
    fi
fi

# Collate one per-query summary into the cohort all-hits + top-hit files.

collate_one() {
    local src="$1" call="$2" top="$3"
    if (( EXCLUDE_SELF )); then
        awk -F'\t' -v sc="${COL_sample}" -v pc="${COL_panel}" \
            'NR>1 && $sc != $pc' "${src}" >> "${call}"
        awk -F'\t' -v sc="${COL_sample}" -v pc="${COL_panel}" \
            'NR>1 && $sc != $pc {print; exit}' "${src}" >> "${top}"
    else
        tail -n +2 "${src}" >> "${call}"
        awk -F'\t' 'NR == 2' "${src}" >> "${top}"
    fi
}

# Collate per-query summary.tsv files (full ranking, non masked).
echo "${SUMMARY_HEADER}" > "${COHORT_TOP}"
echo "${SUMMARY_HEADER}" > "${COHORT_ALL}"

n_total=0
n_failed=0

shopt -s nullglob
for sample_dir in "${WORKDIR}"/*/; do
    summary_tsv="${sample_dir}summary.tsv"

    if [[ ! -e "${summary_tsv}" ]]; then
        # No machine-readable summary — was this sample marked FAILED in the text summary?
        if [[ -e "${sample_dir}summary.txt" ]] \
                && grep -q "^FAILED:" "${sample_dir}summary.txt" 2>/dev/null; then
            n_failed=$((n_failed + 1))
        fi
        continue
    fi

    collate_one "${summary_tsv}" "${COHORT_ALL}" "${COHORT_TOP}"
    n_total=$((n_total + 1))
done
shopt -u nullglob

# Collate per-query summary_var.tsv files (variable-fraction ranking).

n_var_total=0
if (( HAS_VAR )); then
    echo "${SUMMARY_HEADER}" > "${COHORT_TOP_VAR}"
    echo "${SUMMARY_HEADER}" > "${COHORT_ALL_VAR}"
    shopt -s nullglob
    for sample_dir in "${WORKDIR}"/*/; do
        sv="${sample_dir}summary_var.tsv"
        [[ -e "${sv}" ]] || continue
        collate_one "${sv}" "${COHORT_ALL_VAR}" "${COHORT_TOP_VAR}"
        n_var_total=$((n_var_total + 1))
    done
    shopt -u nullglob
fi

# Panel hit frequency (legacy, and var if present).

echo -e "panel\tn_top_hits" > "${COHORT_FREQ}"
if (( n_total > 0 )); then
    tail -n +2 "${COHORT_TOP}" | cut -f"${COL_panel}" \
        | sort | uniq -c | sort -k1,1 -nr \
        | awk '{print $2 "\t" $1}' >> "${COHORT_FREQ}"
fi

if (( HAS_VAR )); then
    echo -e "panel\tn_top_hits" > "${COHORT_FREQ_VAR}"
    if (( n_var_total > 0 )); then
        tail -n +2 "${COHORT_TOP_VAR}" | cut -f"${COL_panel}" \
            | sort | uniq -c | sort -k1,1 -nr \
            | awk '{print $2 "\t" $1}' >> "${COHORT_FREQ_VAR}"
    fi
fi

# Distribution helper: min/q25/median/q75/max/mean of a numeric column.

dist_of() {  # dist_of <tsv_with_header> <col_index> <n_rows>
    local file="$1" col="$2" n="$3"
    if (( n <= 0 )); then echo "- - - - - -"; return; fi
    local tmp; tmp=$(mktemp)
    tail -n +2 "${file}" | cut -f"${col}" | sort -g > "${tmp}"
    awk -v N="${n}" '
        BEGIN {
            i25 = int(N*0.25); if (i25 < 1) i25 = 1
            i50 = int(N*0.5);  if (i50 < 1) i50 = 1
            i75 = int(N*0.75); if (i75 < 1) i75 = 1
        }
        NR == 1   { mn = $1 }
        NR == i25 { q25 = $1 }
        NR == i50 { q50 = $1 }
        NR == i75 { q75 = $1 }
        { sum += $1; if ($1+0 > mx+0) mx = $1 }
        END {
            if (NR > 0) printf "%.4f %.4f %.4f %.4f %.4f %.4f\n", mn, q25, q50, q75, mx, sum/NR
            else        print "- - - - - -"
        }
    ' "${tmp}"
    rm -f "${tmp}"
}

read jmin jq25 jmed jq75 jmax jmean < <(dist_of "${COHORT_TOP}" "${COL_jaccard}" "${n_total}")

jvmin="-"; jvq25="-"; jvmed="-"; jvq75="-"; jvmax="-"; jvmean="-"
if (( HAS_VAR )); then
    read jvmin jvq25 jvmed jvq75 jvmax jvmean \
        < <(dist_of "${COHORT_TOP_VAR}" "${COL_jaccard_var}" "${n_var_total}")
fi

# ---------------------------------------------------------------------------
# 4. Flag tallies (legacy flag, and flag_var if present).
# ---------------------------------------------------------------------------
n_red=$(awk  -F'\t' -v c="${COL_flag}" 'NR>1 && $c=="REDUNDANT"' "${COHORT_TOP}" | wc -l)
n_sub=$(awk  -F'\t' -v c="${COL_flag}" 'NR>1 && $c=="SUBSET"'    "${COHORT_TOP}" | wc -l)
n_none=$(awk -F'\t' -v c="${COL_flag}" 'NR>1 && $c!="REDUNDANT" && $c!="SUBSET"' "${COHORT_TOP}" | wc -l)

nv_red=0; nv_sub=0; nv_none=0
if (( HAS_VAR )); then
    nv_red=$(awk  -F'\t' -v c="${COL_flag_var}" 'NR>1 && $c=="REDUNDANT"' "${COHORT_TOP_VAR}" | wc -l)
    nv_sub=$(awk  -F'\t' -v c="${COL_flag_var}" 'NR>1 && $c=="SUBSET"'    "${COHORT_TOP_VAR}" | wc -l)
    nv_none=$(awk -F'\t' -v c="${COL_flag_var}" 'NR>1 && $c!="REDUNDANT" && $c!="SUBSET"' "${COHORT_TOP_VAR}" | wc -l)
fi

# ---------------------------------------------------------------------------
# 5. Panel sanity stats.
# ---------------------------------------------------------------------------
n_panel=0
n_excl=0
[[ -e "${SIZES}"    ]] && n_panel=$(($(wc -l < "${SIZES}")    - 1))
[[ -e "${EXCLUDED}" ]] && n_excl=$(($(wc -l < "${EXCLUDED}") - 1))
(( n_panel < 0 )) && n_panel=0
(( n_excl  < 0 )) && n_excl=0

U="?"
[[ -e "${U_FILE}" ]] && U=$(cat "${U_FILE}")

n_var_set="?"
[[ -e "${SIZES_VAR}" ]] && n_var_set="present"

pct() { awk -v n="$1" -v d="$2" 'BEGIN{
    if (d>0) printf "%5.1f%%", 100*n/d; else printf "  -  " }'; }

# Make a more readable summary of the results

{
    echo "================================================================"
    echo "  TEFF COHORT REDUNDANCY REPORT"
    echo "================================================================"
    echo "  workspace      : ${PROJ_DIR}"
    echo "  generated      : $(date)"
    echo "  panel size     : ${n_panel} retained (${n_excl} excluded by sanity check)"
    echo "  panel |U|      : ${U}"
    echo "  query count    : ${n_total} processed, ${n_failed} failed/missing"
    if (( EXCLUDE_SELF )); then
        echo "  mode           : --exclude-self (all-vs-all; self-matches dropped, top hit = best non-self)"
    fi
    echo "  thresholds     : jaccard >= ${JACCARD_THRESHOLD}, containment >= ${CONTAINMENT_THRESHOLD}"
    if (( HAS_VAR )); then
        echo "  variable-fraction scoring : ACTIVE (${n_var_total}/${n_total} queries scored)"
    else
        echo "  variable-fraction scoring : not present in this workspace"
    fi
    echo "================================================================"
    echo ""
    echo "### LEGACY  (whole-set Jaccard — same as v5; coverage-sensitive) ###"
    echo ""
    echo "FLAG SUMMARY  (top-1 match per query)"
    printf "  %-12s  %4d  (%s)\n" "REDUNDANT" "${n_red}"  "$(pct ${n_red}  ${n_total})"
    printf "  %-12s  %4d  (%s)\n" "SUBSET"    "${n_sub}"  "$(pct ${n_sub}  ${n_total})"
    printf "  %-12s  %4d  (%s)\n" "no flag"   "${n_none}" "$(pct ${n_none} ${n_total})"
    echo ""
    echo "TOP-MATCH JACCARD DISTRIBUTION  (top-1 per query)"
    printf "  %-8s = %s\n" "min"    "${jmin}"
    printf "  %-8s = %s\n" "q25"    "${jq25}"
    printf "  %-8s = %s\n" "median" "${jmed}"
    printf "  %-8s = %s\n" "q75"    "${jq75}"
    printf "  %-8s = %s\n" "max"    "${jmax}"
    printf "  %-8s = %s\n" "mean"   "${jmean}"
    echo ""
    echo "CONFIRMED REDUNDANT PAIRS  (jaccard >= ${JACCARD_THRESHOLD} AND containment >= ${CONTAINMENT_THRESHOLD})"
    if (( n_red > 0 )); then
        printf "  %-20s  %-20s  %8s  %8s  %8s\n" "query" "panel" "jaccard" "contain." "ssr"
        echo "  --------------------------------------------------------------------------"
        awk -F'\t' -v c="${COL_flag}" '$c == "REDUNDANT"' "${COHORT_TOP}" \
            | sort -t$'\t' -k"${COL_jaccard}","${COL_jaccard}" -gr \
            | awk -F'\t' -v sm="${COL_sample}" -v pn="${COL_panel}" -v jc="${COL_jaccard}" \
                  -v ct="${COL_containment}" -v ss="${COL_ssr}" \
                  '{ printf "  %-20s  %-20s  %8.4f  %8.4f  %8.4f\n", $sm, $pn, $jc, $ct, $ss }'
    else
        echo "  (none)"
    fi
    echo ""
    echo "SUBSET CANDIDATES  (containment >= ${CONTAINMENT_THRESHOLD} but jaccard < ${JACCARD_THRESHOLD})"
    if (( n_sub > 0 )); then
        printf "  %-20s  %-20s  %8s  %8s  %8s\n" "query" "panel" "jaccard" "contain." "ssr"
        echo "  --------------------------------------------------------------------------"
        awk -F'\t' -v c="${COL_flag}" '$c == "SUBSET"' "${COHORT_TOP}" \
            | sort -t$'\t' -k"${COL_jaccard}","${COL_jaccard}" -gr \
            | awk -F'\t' -v sm="${COL_sample}" -v pn="${COL_panel}" -v jc="${COL_jaccard}" \
                  -v ct="${COL_containment}" -v ss="${COL_ssr}" \
                  '{ printf "  %-20s  %-20s  %8.4f  %8.4f  %8.4f\n", $sm, $pn, $jc, $ct, $ss }'
    else
        echo "  (none)"
    fi
    echo ""
    echo "PANEL HIT FREQUENCY  (legacy top-1 most often, top 20)"
    printf "  %-20s  %s\n" "panel" "n_top_hits"
    echo "  ----------------------------------------"
    if (( n_total > 0 )); then
        tail -n +2 "${COHORT_FREQ}" | head -20 \
            | awk -F'\t' '{ printf "  %-20s  %d\n", $1, $2 }'
    else
        echo "  (no queries processed)"
    fi
    echo ""

    # Variable mask run
    
    if (( HAS_VAR )); then
        echo "================================================================"
        echo "### VARIABLE-FRACTION  (core-masked; coverage-robust signal) ###"
        echo "    ranking + flags below use jaccard_var / containment_var,"
        echo "    i.e. only discriminating (non-core) k-mers. This is the view"
        echo "    intended for the low-coverage cross-genebank reconciliation."
        echo "================================================================"
        echo ""
        echo "FLAG SUMMARY  (variable-fraction top-1 match per query)"
        printf "  %-12s  %4d  (%s)\n" "REDUNDANT" "${nv_red}"  "$(pct ${nv_red}  ${n_var_total})"
        printf "  %-12s  %4d  (%s)\n" "SUBSET"    "${nv_sub}"  "$(pct ${nv_sub}  ${n_var_total})"
        printf "  %-12s  %4d  (%s)\n" "no flag"   "${nv_none}" "$(pct ${nv_none} ${n_var_total})"
        echo ""
        echo "TOP-MATCH JACCARD_VAR DISTRIBUTION  (variable-fraction top-1 per query)"
        printf "  %-8s = %s\n" "min"    "${jvmin}"
        printf "  %-8s = %s\n" "q25"    "${jvq25}"
        printf "  %-8s = %s\n" "median" "${jvmed}"
        printf "  %-8s = %s\n" "q75"    "${jvq75}"
        printf "  %-8s = %s\n" "max"    "${jvmax}"
        printf "  %-8s = %s\n" "mean"   "${jvmean}"
        echo "  (if this spread is wide where the legacy one above was flat,"
        echo "   the core-masking has recovered real discriminating signal.)"
        echo ""
        echo "CONFIRMED REDUNDANT PAIRS  (jaccard_var >= ${JACCARD_THRESHOLD} AND containment_var >= ${CONTAINMENT_THRESHOLD})"
        if (( nv_red > 0 )); then
            printf "  %-20s  %-20s  %8s  %8s  %8s\n" "query" "panel" "jac_var" "con_var" "con_q"
            echo "  --------------------------------------------------------------------------"
            awk -F'\t' -v c="${COL_flag_var}" '$c == "REDUNDANT"' "${COHORT_TOP_VAR}" \
                | sort -t$'\t' -k"${COL_jaccard_var}","${COL_jaccard_var}" -gr \
                | awk -F'\t' -v sm="${COL_sample}" -v pn="${COL_panel}" -v jv="${COL_jaccard_var}" \
                      -v cv="${COL_containment_var}" -v cq="${COL_containment_q}" \
                      '{ printf "  %-20s  %-20s  %8.4f  %8.4f  %8.4f\n", $sm, $pn, $jv, $cv, (cq?$cq:0) }'
        else
            echo "  (none)"
        fi
        echo ""
        echo "SUBSET CANDIDATES  (containment_var >= ${CONTAINMENT_THRESHOLD} but jaccard_var < ${JACCARD_THRESHOLD})"
        if (( nv_sub > 0 )); then
            printf "  %-20s  %-20s  %8s  %8s  %8s\n" "query" "panel" "jac_var" "con_var" "con_q"
            echo "  --------------------------------------------------------------------------"
            awk -F'\t' -v c="${COL_flag_var}" '$c == "SUBSET"' "${COHORT_TOP_VAR}" \
                | sort -t$'\t' -k"${COL_jaccard_var}","${COL_jaccard_var}" -gr \
                | awk -F'\t' -v sm="${COL_sample}" -v pn="${COL_panel}" -v jv="${COL_jaccard_var}" \
                      -v cv="${COL_containment_var}" -v cq="${COL_containment_q}" \
                      '{ printf "  %-20s  %-20s  %8.4f  %8.4f  %8.4f\n", $sm, $pn, $jv, $cv, (cq?$cq:0) }'
        else
            echo "  (none)"
        fi
        echo ""
        echo "PANEL HIT FREQUENCY  (variable-fraction top-1 most often, top 20)"
        printf "  %-20s  %s\n" "panel" "n_top_hits"
        echo "  ----------------------------------------"
        if (( n_var_total > 0 )); then
            tail -n +2 "${COHORT_FREQ_VAR}" | head -20 \
                | awk -F'\t' '{ printf "  %-20s  %d\n", $1, $2 }'
        else
            echo "  (no queries with variable-fraction scoring)"
        fi
        echo ""
        if (( n_var_total < n_total )); then
            echo "NOTE: ${n_var_total} of ${n_total} queries had variable-fraction scoring."
            echo "      The remainder fell back to legacy-only (mask inactive for that query)."
            echo ""
        fi
    fi

    if (( n_excl > 0 )); then
        echo "EXCLUDED PANEL ACCESSIONS  (failed sanity check)"
        printf "  %-20s  %-14s  %s\n" "accession" "n_kmers" "reason"
        echo "  --------------------------------------------------"
        tail -n +2 "${EXCLUDED}" | awk -F'\t' '{ printf "  %-20s  %-14s  %s\n", $1, $2, $3 }'
        echo ""
    fi
    echo "================================================================"
    echo "  Generated files:"
    echo "    ${COHORT_TOP}"
    echo "    ${COHORT_ALL}"
    echo "    ${COHORT_FREQ}"
    if (( HAS_VAR )); then
        echo "    ${COHORT_TOP_VAR}"
        echo "    ${COHORT_ALL_VAR}"
        echo "    ${COHORT_FREQ_VAR}"
    fi
    echo "    ${COHORT_REPORT}"
    echo "================================================================"
} > "${COHORT_REPORT}"

# Print to stdout too, so the SLURM job log captures it
cat "${COHORT_REPORT}"
