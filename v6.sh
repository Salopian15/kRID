#!/bin/bash
# ============================================================================
# v6.s
#
# Self-contained orchestrator: writes all SLURM scripts (including the
# per-sample query pipeline) into <project_workspace>/scripts/, then
# submits them as a dependency chain.
#
# Usage:
#   ./v6.sh [--ignore-corrupt] <panel_input> <query_input> <project_workspace>
#
# Supports:
#   - FASTQ files (trimmed via fastp, kmc min_count=4)
#   - FASTA files / Assemblies (no fastp, kmc min_count=1)
#   - Precompiled KMC databases (*.kmc_pre / *.kmc_suf)
#
# Note to future Jack, pipeline is still specific to JIC compute nodes, 
# this will need fixing/changing in future, sincerely past jack
# ============================================================================
set -euo pipefail

# ---------------------------------------------------------------------------
# Logging helpers
# ---------------------------------------------------------------------------
log()  { echo "[$(date '+%H:%M:%S')] $*"; }
warn() { echo "[$(date '+%H:%M:%S')] WARN: $*" >&2; }
die()  { echo "[$(date '+%H:%M:%S')] ERROR: $*" >&2; exit 1; }

# ---------------------------------------------------------------------------
# Argument Parsing
# ---------------------------------------------------------------------------
IGNORE_CORRUPT=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --ignore-corrupt)
            IGNORE_CORRUPT=1
            shift
            ;;
        -*)
            die "Unknown option: $1"
            ;;
        *)
            break
            ;;
    esac
done

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 [--ignore-corrupt] <panel_input_dir> <query_input_dir> <project_workspace>"
    echo ""
    echo "  --ignore-corrupt  : Discard corrupted samples and continue pipeline."
    echo "  panel_input_dir   : directory containing FASTQ, FASTA, or precompiled KMC DBs"
    echo "  query_input_dir   : directory containing FASTQ, FASTA, or precompiled KMC DBs"
    echo "  project_workspace : output directory (created if it doesn't exist)"
    echo ""
    echo "Environment knobs:"
    echo "  MIN_PANEL_KMERS         minimum k-mers for a panel accession to be kept (default 1000000)"
    echo "  JACCARD_THRESHOLD       Jaccard cutoff for REDUNDANT flag (default 0.85)"
    echo "  CONTAINMENT_THRESHOLD   Containment cutoff for REDUNDANT/SUBSET flag (default 0.95)"
    echo "  QUERY_CI                FASTQ query k-mer min-count (default 2; v5 used 4)"
    echo "  PANEL_CI                FASTQ panel-build k-mer min-count (default 4)"
    echo "  MASK_CORE               also score on the variable k-mer fraction (default 1)"
    echo "  CORE_FREQ_FRAC          k-mers in >= this frac of accessions are 'core' (default 0.95)"
    echo "  MASK_PAR                parallel kmc_tools workers in the mask build (default 12)"
    exit 1
fi

# ---------------------------------------------------------------------------
# Pre-run: check if on a SLURM head node
# ---------------------------------------------------------------------------
command -v sbatch >/dev/null 2>&1 \
    || die "sbatch not found — this script must be run on a SLURM head node."
command -v squeue >/dev/null 2>&1 \
    || die "squeue not found — SLURM tools are not available on this node."

# ---------------------------------------------------------------------------
# Resolve paths
# ---------------------------------------------------------------------------
PANEL_DIR=$(realpath "$1" 2>/dev/null) \
    || die "Cannot resolve panel_fastq_dir '$1' — does it exist?"
QUERY_DIR=$(realpath "$2" 2>/dev/null) \
    || die "Cannot resolve query_fastq_dir '$2' — does it exist?"

[[ -d "$PANEL_DIR" ]] || die "Panel directory does not exist: $PANEL_DIR"
[[ -d "$QUERY_DIR" ]] || die "Query directory does not exist: $QUERY_DIR"
[[ -r "$PANEL_DIR" ]] || die "Panel directory is not readable: $PANEL_DIR"
[[ -r "$QUERY_DIR" ]] || die "Query directory is not readable: $QUERY_DIR"

PROJ_DIR=$(realpath "$3")

PANEL_DBS="${PROJ_DIR}/panel_dbs"
PANEL_UNION_DIR="${PROJ_DIR}/panel_union"
WORKDIR="${PROJ_DIR}/query_runs"
SCRIPTS_DIR="${PROJ_DIR}/scripts"

mkdir -p "$PANEL_DBS" "$PANEL_UNION_DIR" "$WORKDIR" "$SCRIPTS_DIR" \
    || die "Failed to create project directories under $PROJ_DIR — check permissions."

log "Project workspace : $PROJ_DIR"
log "Panel input       : $PANEL_DIR"
log "Query input       : $QUERY_DIR"
if [ "$IGNORE_CORRUPT" -eq 1 ]; then
    log "Ignore Corrupt    : ENABLED (Corrupt samples will be skipped gracefully)"
fi

# ---------------------------------------------------------------------------
# Helper: validate job ID returned by sbatch --parsable
# ---------------------------------------------------------------------------
validate_job_id() {
    local jid="$1"
    local phase="$2"
    [[ "$jid" =~ ^[0-9]+$ ]] \
        || die "Phase '$phase' submission failed — sbatch returned an unexpected value: '$jid'. Check the SLURM logs."
    log "  Job ID: $jid"
}

# ---------------------------------------------------------------------------
# Helper: Generate Sample Sheet from Directory
# ---------------------------------------------------------------------------
generate_sheet() {
    local in_dir=$1
    local out_sheet=$2
    local label=$3

    echo -e "sample_id\ttype\tr1\tr2" > "$out_sheet"
    local found=0

    declare -A processed
    shopt -s nullglob

    # 1. Look for precompiled KMC files
    for db in "${in_dir}"/*.kmc_pre; do
        local base="${db%.kmc_pre}"
        local sid=$(basename "$base")
        if [[ -z "${processed[$sid]:-}" ]]; then
            echo -e "${sid}\tKMC\t${base}\t" >> "$out_sheet"
            processed[$sid]=1
            (( found++ )) || true
        fi
    done

    # 2. Look for FASTA formats (Genome assemblies)
    for f in "${in_dir}"/*.fa "${in_dir}"/*.fasta "${in_dir}"/*.fna "${in_dir}"/*.fa.gz "${in_dir}"/*.fasta.gz; do
        local sid=$(basename "$f" | sed -E 's/\.(fa|fasta|fna)(\.gz)?$//')
        if [[ -z "${processed[$sid]:-}" ]]; then
            echo -e "${sid}\tFASTA\t${f}\t" >> "$out_sheet"
            processed[$sid]=1
            (( found++ )) || true
        fi
    done

    # 3. Look for FASTQ formats
    for f in "${in_dir}"/*.fastq "${in_dir}"/*.fq "${in_dir}"/*.fastq.gz "${in_dir}"/*.fq.gz; do
        # Exclude R2 to avoid double processing
        local is_r2=0
        if [[ "$f" =~ _R2 || "$f" =~ _2\.f ]]; then is_r2=1; fi

        if [[ $is_r2 -eq 1 ]]; then continue; fi

        local sid=$(basename "$f" | sed -E 's/_(R1|1).*//; s/\.(fastq|fq)(\.gz)?$//')

        if [[ -z "${processed[$sid]:-}" ]]; then
            local r2=""
            if [[ "$f" =~ _R1 ]]; then r2="${f/_R1/_R2}"; fi
            if [[ "$f" =~ _1\.f ]]; then r2="${f/_1./_2.}"; fi

            if [[ -n "$r2" && -e "$r2" ]]; then
                echo -e "${sid}\tFASTQ\t${f}\t${r2}" >> "$out_sheet"
            else
                echo -e "${sid}\tFASTQ\t${f}\t" >> "$out_sheet"
            fi
            processed[$sid]=1
            (( found++ )) || true
        fi
    done
    shopt -u nullglob

    if [[ $found -eq 0 ]]; then
        die "No valid input files found in $label directory: $in_dir"
    fi

    log "Found $found $label sample(s)."
}

PANEL_SHEET="${PROJ_DIR}/panel_samples.tsv"
QUERY_SHEET="${PROJ_DIR}/query_samples.tsv"

log "Scanning input directories..."
generate_sheet "$PANEL_DIR" "$PANEL_SHEET" "panel"
generate_sheet "$QUERY_DIR" "$QUERY_SHEET" "query"

PANEL_COUNT=$(($(wc -l < "$PANEL_SHEET") - 1))
QUERY_COUNT=$(($(wc -l < "$QUERY_SHEET") - 1))

log "Panel: $PANEL_COUNT sample(s)  |  Query: $QUERY_COUNT sample(s)"

# ============================================================================
# Write SLURM scripts, note to future jack this is still specific to JIC nodes, sincerely past Jack
# ============================================================================
log "Writing SLURM scripts to $SCRIPTS_DIR ..."

# ---------------------------------------------------------------------------
# Phase 0: Per-sample query pipeline
# ---------------------------------------------------------------------------
cat << 'EOF_SLURM_HEAD' > "${SCRIPTS_DIR}/teff_query_pipeline.slurm"
#!/bin/bash -e
#SBATCH -p jic-medium,nbi-medium
#SBATCH -t 06:00:00
#SBATCH --mem=128G
#SBATCH -c 16
#SBATCH -J teff_query
#SBATCH -o slurm-%A_%a.out

EOF_SLURM_HEAD

echo "IGNORE_CORRUPT=${IGNORE_CORRUPT}" >> "${SCRIPTS_DIR}/teff_query_pipeline.slurm"

cat << 'QUERY_PIPELINE_EOF' >> "${SCRIPTS_DIR}/teff_query_pipeline.slurm"

_stage="init"
trap 'echo "ERROR [${BASH_SOURCE[0]}] stage=${_stage} line=${LINENO}: command exited with status $?" >&2' ERR

module load fastp kmc || echo "WARN: module load failed, assuming binaries are in PATH"

EXISTING_DBS_DIR="${EXISTING_DBS_DIR:-${HOME}/scratch/teff_kmcdatabases}"
PANEL_UNION_DB="${PANEL_UNION_DB:-${HOME}/scratch/teff/unified_teff_kmers}"
HEADER_FILE="${HEADER_FILE:-${HOME}/scratch/teff/unified_list.txt}"
EXISTING_SIZES="${EXISTING_SIZES:-${HOME}/scratch/teff/results/accession_sizes.tsv}"
PANEL_U_FILE="${PANEL_U_FILE:-${HOME}/scratch/teff/results/panel_union_size.txt}"
WORKDIR="${WORKDIR:-${PWD}/query_runs}"

# v6: variable-fraction (core-masked) scoring inputs. VAR_DB is the panel
# variable k-mer set built in Phase 2; EXISTING_SIZES_VAR holds |accession n var|.
# If MASK_CORE=0 or VAR_DB is absent, all *_var metrics are reported as NA and
# the run behaves exactly like v5.
MASK_CORE="${MASK_CORE:-1}"
VAR_DB="${VAR_DB:-}"
EXISTING_SIZES_VAR="${EXISTING_SIZES_VAR:-}"

KMER_LENGTH=31
FASTP_THREADS=8
FASTP_QUAL=20
FASTP_MIN_LEN=50
KMC_MEMORY_GB=128
KMC_THREADS=16
TOP_N=20

# Thresholds (configurable via environment)
JACCARD_THRESHOLD="${JACCARD_THRESHOLD:-0.85}"
CONTAINMENT_THRESHOLD="${CONTAINMENT_THRESHOLD:-0.95}"
REDUNDANCY_THRESHOLD="${REDUNDANCY_THRESHOLD:-0.97}"   # legacy SSR threshold (diagnostic only)
QUERY_CI="${QUERY_CI:-2}"                              # v6: FASTQ k-mer min-count (v5 used 4)
PAIR_THREADS=16

set -euo pipefail

usage() {
    echo "Usage: sbatch --array=1-N teff_query_pipeline.slurm --sample-sheet samples.tsv"
    echo "  --jaccard-threshold      override JACCARD_THRESHOLD"
    echo "  --containment-threshold  override CONTAINMENT_THRESHOLD"
    echo "  --threshold              legacy alias for --shared-state-threshold"
    exit 1
}

SAMPLE_ID=""; TYPE=""; R1=""; R2=""; SAMPLE_SHEET=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --sample-sheet)           SAMPLE_SHEET="$2"; shift 2 ;;
        --jaccard-threshold)      JACCARD_THRESHOLD="$2"; shift 2 ;;
        --containment-threshold)  CONTAINMENT_THRESHOLD="$2"; shift 2 ;;
        --threshold|--shared-state-threshold)
                                  REDUNDANCY_THRESHOLD="$2"; shift 2 ;;
        -h|--help)                usage ;;
        *) echo "ERROR: Unknown argument: $1" >&2; usage ;;
    esac
done

if [[ -n "${SAMPLE_SHEET}" ]]; then
    line=$(tail -n +2 "${SAMPLE_SHEET}" | sed -n "${SLURM_ARRAY_TASK_ID}p")
    SAMPLE_ID=$(echo "${line}" | cut -f1)
    TYPE=$(echo "${line}" | cut -f2)
    R1=$(echo "${line}" | cut -f3)
    R2=$(echo "${line}" | cut -f4)
fi

[[ -n "${SAMPLE_ID}" ]] || { echo "ERROR: --sample-sheet requires array task execution." >&2; exit 1; }

sample_dir="${WORKDIR}/${SAMPLE_ID}"
mkdir -p "${sample_dir}" || { echo "ERROR: Failed to create working directory: ${sample_dir}" >&2; exit 1; }
cd "${sample_dir}"

echo "============================================================"
echo "  teff query pipeline (fv4)"
echo "============================================================"
echo "  sample_id : ${SAMPLE_ID}"
echo "  type      : ${TYPE}"
echo "  R1/DB     : ${R1}"
echo "  job       : ${SLURM_JOB_ID:-N/A} / task ${SLURM_ARRAY_TASK_ID:-N/A}"
echo "  thresholds: jaccard>=${JACCARD_THRESHOLD}  containment>=${CONTAINMENT_THRESHOLD}"
echo "  query_ci  : ${QUERY_CI}   mask_core: ${MASK_CORE}"
echo "============================================================"

# STAGE 1 — QC + Trim (FASTQ ONLY)
_stage="fastp-trim"
echo; echo ">>> [1/4] Prep / Trim"

trim_R1=""
trim_R2=""

if [[ "$TYPE" == "FASTQ" ]]; then
    mkdir -p trim
    trim_R1="trim/${SAMPLE_ID}_R1.fq.gz"
    [[ -n "${R2}" ]] && trim_R2="trim/${SAMPLE_ID}_R2.fq.gz"

    if [[ ! -e "${trim_R1}" ]]; then
        if [[ -n "${R2}" ]]; then
            if ! fastp -i "${R1}" -I "${R2}" -o "${trim_R1}" -O "${trim_R2}" -j "trim/${SAMPLE_ID}.fastp.json" -h "trim/${SAMPLE_ID}.fastp.html" -q ${FASTP_QUAL} -l ${FASTP_MIN_LEN} -w ${FASTP_THREADS}; then
                if [[ "${IGNORE_CORRUPT}" == "1" ]]; then
                    echo "WARN: fastp failed. Ignoring due to --ignore-corrupt."
                    echo "FAILED: Corrupted input at fastp" > summary.txt
                    exit 0
                else
                    echo "ERROR: fastp failed." >&2; exit 1
                fi
            fi
        else
            if ! fastp -i "${R1}" -o "${trim_R1}" -j "trim/${SAMPLE_ID}.fastp.json" -h "trim/${SAMPLE_ID}.fastp.html" -q ${FASTP_QUAL} -l ${FASTP_MIN_LEN} -w ${FASTP_THREADS}; then
                if [[ "${IGNORE_CORRUPT}" == "1" ]]; then
                    echo "WARN: fastp failed. Ignoring due to --ignore-corrupt."
                    echo "FAILED: Corrupted input at fastp" > summary.txt
                    exit 0
                else
                    echo "ERROR: fastp failed." >&2; exit 1
                fi
            fi
        fi
    fi
elif [[ "$TYPE" == "FASTA" ]]; then
    echo "    [skip] FASTA input, skipping fastp"
    trim_R1="${R1}"
elif [[ "$TYPE" == "KMC" ]]; then
    echo "    [skip] Precompiled KMC input, skipping fastp"
fi

# STAGE 2 — KMC count
_stage="kmc-count"
echo; echo ">>> [2/4] KMC count"
mkdir -p kmc tmp_kmc
query_db="kmc/${SAMPLE_ID}_ci${QUERY_CI}"

if [[ "$TYPE" == "KMC" ]]; then
    echo "    Linking existing KMC database..."
    ln -sf "${R1}.kmc_pre" "${query_db}.kmc_pre"
    ln -sf "${R1}.kmc_suf" "${query_db}.kmc_suf"
else
    if [[ ! -e "${query_db}.kmc_pre" ]]; then
        kmc_input="kmc/input_${SAMPLE_ID}.txt"
        echo "${trim_R1}" > "${kmc_input}"
        [[ -n "${trim_R2}" ]] && echo "${trim_R2}" >> "${kmc_input}"

        if [[ "$TYPE" == "FASTA" ]]; then
            if ! kmc -k${KMER_LENGTH} -ci1 -fm -m${KMC_MEMORY_GB} -t${KMC_THREADS} @"${kmc_input}" "${query_db}" tmp_kmc/; then
                if [[ "${IGNORE_CORRUPT}" == "1" ]]; then
                    echo "WARN: kmc failed on FASTA. Ignoring due to --ignore-corrupt."
                    echo "FAILED: Corrupted input at kmc count" > summary.txt
                    exit 0
                else
                    exit 1
                fi
            fi
        elif [[ "$TYPE" == "FASTQ" ]]; then
            if ! kmc -k${KMER_LENGTH} -ci${QUERY_CI} -fq -m${KMC_MEMORY_GB} -t${KMC_THREADS} @"${kmc_input}" "${query_db}" tmp_kmc/; then
                if [[ "${IGNORE_CORRUPT}" == "1" ]]; then
                    echo "WARN: kmc failed on FASTQ. Ignoring due to --ignore-corrupt."
                    echo "FAILED: Corrupted input at kmc count" > summary.txt
                    exit 0
                else
                    exit 1
                fi
            fi
        fi
        rm -rf tmp_kmc/* "${kmc_input}"
    fi
fi

count_kmers() { kmc_tools transform "$1" dump /dev/stdout 2>/dev/null | wc -l; }

echo "    counting query k-mers..."
N_B=$(count_kmers "${query_db}")
echo "    |query| = ${N_B} k-mers"

[[ $N_B -gt 0 ]] || { echo "ERROR: KMC DB empty." >&2; exit 1; }

# v6: decide whether variable-fraction (core-masked) scoring is available for
# this run, and if so count the query's variable k-mers (nb_var = |query n VAR|).
MASK_ACTIVE=0
N_B_VAR="NA"
if [[ "${MASK_CORE}" == "1" && -n "${VAR_DB}" && -e "${VAR_DB}.kmc_pre" ]]; then
    MASK_ACTIVE=1
    qvar_db="kmc/${SAMPLE_ID}_ci${QUERY_CI}_var"
    if kmc_tools simple "${query_db}" "${VAR_DB}" intersect "${qvar_db}" >/dev/null 2>&1; then
        N_B_VAR=$(count_kmers "${qvar_db}")
        rm -f "${qvar_db}.kmc_pre" "${qvar_db}.kmc_suf"
    else
        echo "    WARN: could not intersect query with VAR_DB; var metrics -> NA"
        MASK_ACTIVE=0
    fi
    echo "    |query n variable| = ${N_B_VAR} k-mers"
else
    echo "    [skip] core-masking inactive (MASK_CORE=${MASK_CORE}, VAR_DB='${VAR_DB}') — legacy metrics only"
fi
export MASK_ACTIVE VAR_DB

# STAGE 3 — establish sizes
_stage="panel-statistics"
echo; echo ">>> [3/4] gather panel statistics"

if [[ -e "${PANEL_U_FILE}" ]]; then
    U=$(cat "${PANEL_U_FILE}")
else
    U=$(count_kmers "${PANEL_UNION_DB}")
    echo "${U}" > "${PANEL_U_FILE}" || true
fi

declare -A SIZE_OF
while IFS=$'\t' read -r sname snk; do
    SIZE_OF["${sname}"]="${snk}"
done < <(tail -n +2 "${EXISTING_SIZES}")

# v6: per-accession variable-fraction sizes (na_var), if the mask was built.
declare -A SIZE_VAR_OF
if [[ "${MASK_ACTIVE}" == "1" && -n "${EXISTING_SIZES_VAR}" && -e "${EXISTING_SIZES_VAR}" ]]; then
    while IFS=$'\t' read -r sname snk; do
        SIZE_VAR_OF["${sname}"]="${snk}"
    done < <(tail -n +2 "${EXISTING_SIZES_VAR}")
fi

# STAGE 4 — pairwise comparison
_stage="pairwise-comparison"
echo; echo ">>> [4/4] pairwise comparison"
mkdir -p inter

compare_one() {
    local ename="$1"
    local edb="${EXISTING_DBS_DIR}/${ename}"
    if [[ ! -e "${edb}.kmc_pre" ]]; then
        echo -e "${ename}\tNA\tNA"
        return
    fi
    local idb="inter/${ename}_v_${SAMPLE_ID}"
    kmc_tools simple "${edb}" "${query_db}" intersect "${idb}" >/dev/null 2>&1 || { echo -e "${ename}\tNA\tNA"; return; }
    local n_inter=$(kmc_tools transform "${idb}" dump /dev/stdout 2>/dev/null | wc -l)
    # v6: variable-fraction intersection = (panel n query) n VAR. Reuses the
    # intersection db we just built, so it's one extra intersect per pair.
    local n_inter_var="NA"
    if [[ "${MASK_ACTIVE}" == "1" ]]; then
        local ivar="inter/${ename}_v_${SAMPLE_ID}_var"
        if kmc_tools simple "${idb}" "${VAR_DB}" intersect "${ivar}" >/dev/null 2>&1; then
            n_inter_var=$(kmc_tools transform "${ivar}" dump /dev/stdout 2>/dev/null | wc -l)
            rm -f "${ivar}.kmc_pre" "${ivar}.kmc_suf"
        fi
    fi
    rm -f "${idb}.kmc_pre" "${idb}.kmc_suf"
    echo -e "${ename}\t${n_inter}\t${n_inter_var}"
}
export -f compare_one
export EXISTING_DBS_DIR query_db SAMPLE_ID MASK_ACTIVE VAR_DB

intersect_raw="intersections.tsv"
cat "${HEADER_FILE}" | grep -v '^$' | xargs -I{} -P ${PAIR_THREADS} bash -c 'compare_one "$@"' _ {} > "${intersect_raw}"
rmdir inter 2>/dev/null || true

_stage="scoring"
echo; echo ">>> writing rates and summary"
rates_tsv="rates.tsv"
rates_var_tsv="rates_var.tsv"
rates_hdr="rank\texisting\tn_a\tn_b\tn_intersect\tn_union\tjaccard\tcontainment\tshared_state_rate\tcontainment_q\tna_var\tnb_var\tn_intersect_var\tjaccard_var\tcontainment_var"

# Compute per-pair metrics. Legacy (whole-set) Jaccard/Containment/SSR are
# unchanged from v5; v6 appends a directional query containment (ni/nb) and,
# when the mask is active, the variable-fraction metrics.
#   jaccard_var      = ni_var / (na_var + nb_var - ni_var)
#   containment_var  = ni_var / nb_var          (query-in-panel, directional)
rates_body=$(mktemp)
awk -F'\t' -v U="${U}" -v NB="${N_B}" -v NBV="${N_B_VAR}" -v mask="${MASK_ACTIVE}" \
    -v sizes_file="${EXISTING_SIZES}" -v sizes_var_file="${EXISTING_SIZES_VAR}" '
    BEGIN {
        while ((getline line < sizes_file) > 0) {
            n = split(line, f, "\t"); if (f[1] == "sample") continue; size[f[1]] = f[2] + 0
        }
        close(sizes_file)
        have_var = (mask == 1 && NBV != "NA" && sizes_var_file != "")
        if (have_var) {
            while ((getline line < sizes_var_file) > 0) {
                n = split(line, g, "\t"); if (g[1] == "sample") continue; size_var[g[1]] = g[2] + 0
            }
            close(sizes_var_file)
            nbv = NBV + 0
        }
    }
    {
        ename = $1; ni = $2; niv = $3
        if (ni == "NA" || !(ename in size)) next
        na = size[ename]
        if (na <= 0 || NB <= 0) next
        nu = na + NB - ni
        if (nu <= 0) next
        jaccard       = ni / nu
        min_ab        = (na < NB) ? na : NB
        containment   = (min_ab > 0) ? ni / min_ab : 0
        containment_q = (NB > 0) ? ni / NB : 0
        ssr           = (U > 0) ? (2*ni + U - na - NB) / U : 0
        legacy = sprintf("%s\t%d\t%d\t%d\t%d\t%.6f\t%.6f\t%.6f\t%.6f", \
                         ename, na, NB, ni, nu, jaccard, containment, ssr, containment_q)
        if (have_var && niv != "NA" && (ename in size_var)) {
            nav = size_var[ename]
            nuv = nav + nbv - niv
            jv  = (nuv > 0) ? niv / nuv : 0
            cv  = (nbv > 0) ? niv / nbv : 0
            printf "%s\t%d\t%d\t%d\t%.6f\t%.6f\n", legacy, nav, nbv, niv, jv, cv
        } else {
            printf "%s\tNA\tNA\tNA\tNA\tNA\n", legacy
        }
    }
' "${intersect_raw}" > "${rates_body}"

# Legacy view: ranked by whole-set Jaccard (field 6) — identical ordering to v5.
{ echo -e "${rates_hdr}"; sort -t$'\t' -k6,6 -gr "${rates_body}" | awk -F'\t' 'BEGIN{OFS="\t"} {print NR, $0}'; } > "${rates_tsv}"

# v6 view: ranked by variable-fraction Jaccard (field 13) so the best masked
# match leads, even if it isn't the legacy top-1.
if [[ "${MASK_ACTIVE}" == "1" ]]; then
    { echo -e "${rates_hdr}"; sort -t$'\t' -k13,13 -gr "${rates_body}" | awk -F'\t' 'BEGIN{OFS="\t"} {print NR, $0}'; } > "${rates_var_tsv}"
fi
rm -f "${rates_body}"

# --- Machine-readable summary (one row per top hit, easy to aggregate) ----
# v6: legacy column positions are preserved ($7 jaccard, $8 containment), so
# anything parsing the v5 summary still works; new columns and a variable-
# fraction flag are appended at the end.
emit_summary() {
    local src="$1" out="$2"
    {
        echo -e "sample\trank\tpanel\tna\tnb\tn_intersect\tn_union\tjaccard\tcontainment\tshared_state_rate\tcontainment_q\tna_var\tnb_var\tn_intersect_var\tjaccard_var\tcontainment_var\tflag\tflag_var"
        awk -F'\t' -v sid="${SAMPLE_ID}" -v top="${TOP_N}" \
            -v jt="${JACCARD_THRESHOLD}" -v ct="${CONTAINMENT_THRESHOLD}" '
            NR == 1 { next }
            NR <= top + 1 {
                jaccard     = $7 + 0
                containment = $8 + 0
                flag = "."
                if (jaccard >= jt && containment >= ct) flag = "REDUNDANT"
                else if (containment >= ct)              flag = "SUBSET"
                flag_var = "."
                if ($14 != "NA") {
                    jv = $14 + 0; cv = $15 + 0
                    if (jv >= jt && cv >= ct) flag_var = "REDUNDANT"
                    else if (cv >= ct)         flag_var = "SUBSET"
                }
                print sid "\t" $0 "\t" flag "\t" flag_var
            }
        ' "${src}"
    } > "${out}"
}

summary_tsv="summary.tsv"
emit_summary "${rates_tsv}" "${summary_tsv}"               # legacy ranking (whole-set Jaccard)
if [[ "${MASK_ACTIVE}" == "1" && -e "${rates_var_tsv}" ]]; then
    emit_summary "${rates_var_tsv}" "summary_var.tsv"      # variable-fraction ranking
fi

# --- Human-readable summary (with z-scores for diagnostic ssr) ------------
# v6: factored so we can emit both the legacy (whole-set) ranking and, when the
# mask is active, a variable-fraction ranking. summary.txt keeps v5 ordering.
emit_human() {
    # Args: src out mode [zcol] [zlabel]
    #   zcol   = field in the rates file to express as a z-score (default 7 =
    #            whole-set jaccard). The z is now built on a PRESENCE-BASED
    #            metric (jaccard / jaccard_var), not the old SSR, so it measures
    #            how far the top hit's jaccard stands above this query's panel-
    #            wide mean — a flat (coverage-artifact) field gives the top hit
    #            z~0, a real duplicate spikes high.
    #   zlabel = column header for that z (default jac_z).
    local src="$1" out="$2" mode="$3" zcol="${4:-7}" zlabel="${5:-jac_z}"
    {
        echo "============================================================"
        echo "  TEFF QUERY RESULT — ${SAMPLE_ID}   [${mode}]"
        echo "============================================================"
        echo "  date            : $(date)"
        echo "  query k-mers    : ${N_B}"
        echo "  query var k-mers: ${N_B_VAR}"
        echo "  panel |U|       : ${U}"
        echo "  jaccard cutoff  : ${JACCARD_THRESHOLD}"
        echo "  contain. cutoff : ${CONTAINMENT_THRESHOLD}"
        echo "  query_ci        : ${QUERY_CI}    mask_core: ${MASK_CORE} (active=${MASK_ACTIVE})"
        echo "  ${zlabel}           : SDs above this query's panel-wide mean (presence-based)"
        echo ""
        printf "  %-5s  %-20s  %8s  %8s  %8s  %8s  %8s  %s\n" \
            "rank" "panel" "jaccard" "contain." "jac_var" "con_var" "${zlabel}" "flag"
        echo "  ------------------------------------------------------------------------------------------"
        awk -F'\t' -v top="${TOP_N}" -v jt="${JACCARD_THRESHOLD}" -v ct="${CONTAINMENT_THRESHOLD}" -v zc="${zcol}" '
            # Pass 1: mean/sd of the chosen presence-based column (zc) across the
            # whole panel for this query, to express the top hits as a z-score.
            NR == FNR {
                if (FNR == 1) next
                if ($zc == "NA") next
                r = $zc + 0; n++; sum += r; sumsq += r*r; next
            }
            FNR == 1 {
                mean = (n>0) ? sum/n : 0
                v    = (n>0) ? (sumsq/n - mean*mean) : 0
                sd   = (v > 0) ? sqrt(v) : 0
                next
            }
            FNR <= top + 1 {
                jaccard     = $7 + 0
                containment = $8 + 0
                zval = $zc
                if (zval == "NA") z_s = "    NA  "
                else              z_s = sprintf("%+8.2f", (sd > 0) ? ((zval + 0) - mean) / sd : 0)
                jv_s = ($14 == "NA") ? "    NA  " : sprintf("%8.4f", $14 + 0)
                cv_s = ($15 == "NA") ? "    NA  " : sprintf("%8.4f", $15 + 0)
                flag = ""
                if ($14 != "NA") {
                    jv = $14 + 0; cv = $15 + 0
                    if (jv >= jt && cv >= ct) flag = "REDUNDANT"
                    else if (cv >= ct)         flag = "SUBSET"
                } else {
                    if (jaccard >= jt && containment >= ct) flag = "REDUNDANT"
                    else if (containment >= ct)              flag = "SUBSET"
                }
                printf "  %-5s  %-20s  %8.4f  %8.4f  %s  %s  %s  %s\n", \
                    $1, $2, jaccard, containment, jv_s, cv_s, z_s, flag
            }
        ' "${src}" "${src}"
        echo "============================================================"
    } > "${out}"
}

summary_txt="summary.txt"
# Legacy report: z on whole-set jaccard (field 7). Variable-fraction report:
# z on jaccard_var (field 14), so each report's z tracks the metric it ranks by.
emit_human "${rates_tsv}" "${summary_txt}" "ranked by whole-set jaccard (legacy)" 7 "jac_z"
if [[ "${MASK_ACTIVE}" == "1" && -e "${rates_var_tsv}" ]]; then
    emit_human "${rates_var_tsv}" "summary_var.txt" "ranked by variable-fraction jaccard (v6)" 14 "jvar_z"
fi

echo "Done."
QUERY_PIPELINE_EOF
chmod +x "${SCRIPTS_DIR}/teff_query_pipeline.slurm"

# ---------------------------------------------------------------------------
# Phase 1: Panel DB Builder
# ---------------------------------------------------------------------------
cat << 'EOF_SLURM_P1' > "${SCRIPTS_DIR}/01_build_panel.slurm"
#!/bin/bash -e
#SBATCH -p jic-medium,nbi-medium
#SBATCH -t 04:00:00
#SBATCH --mem=128G
#SBATCH -c 16
#SBATCH -J teff_panel
#SBATCH -o slurm-%A_%a.out
#SBATCH --requeue
EOF_SLURM_P1

echo "IGNORE_CORRUPT=${IGNORE_CORRUPT}" >> "${SCRIPTS_DIR}/01_build_panel.slurm"

cat << 'EOF_P1' >> "${SCRIPTS_DIR}/01_build_panel.slurm"

module load fastp kmc || echo "WARN: module load failed"
set -euo pipefail

trap 'echo "ERROR [01_build_panel] task=${SLURM_ARRAY_TASK_ID:-N/A} line=${LINENO}: exited with status $?" >&2' ERR

SHEET="$1"
OUTDIR="$2"
KMER_LENGTH=31
KMC_THREADS=16
PANEL_CI="${PANEL_CI:-4}"   # v6: FASTQ panel-build min-count (default 4, as v5)

line=$(tail -n +2 "${SHEET}" | sed -n "${SLURM_ARRAY_TASK_ID}p")
SAMPLE_ID=$(echo "$line" | cut -f1)
TYPE=$(echo "$line" | cut -f2)
R1=$(echo "$line" | cut -f3)
R2=$(echo "$line" | cut -f4)

echo "Building KMC database for panel sample: ${SAMPLE_ID} (Type: ${TYPE})"

if [[ -e "${OUTDIR}/${SAMPLE_ID}.kmc_pre" && -e "${OUTDIR}/${SAMPLE_ID}.kmc_suf" ]]; then
    echo "    [skip] KMC database already exists"
else
    if [[ "$TYPE" == "KMC" ]]; then
        echo "    Precompiled DB found, symlinking..."
        ln -sf "${R1}.kmc_pre" "${OUTDIR}/${SAMPLE_ID}.kmc_pre"
        ln -sf "${R1}.kmc_suf" "${OUTDIR}/${SAMPLE_ID}.kmc_suf"
    else
        mkdir -p "${OUTDIR}/tmp_${SAMPLE_ID}"
        kmc_input="${OUTDIR}/input_${SAMPLE_ID}.txt"
        echo "$R1" > "$kmc_input"
        [[ -n "$R2" && -e "$R2" ]] && echo "$R2" >> "$kmc_input"

        if [[ "$TYPE" == "FASTA" ]]; then
            if ! kmc -k${KMER_LENGTH} -ci1 -fm -m128 -t${KMC_THREADS} @"${kmc_input}" "${OUTDIR}/${SAMPLE_ID}" "${OUTDIR}/tmp_${SAMPLE_ID}/"; then
                if [[ "${IGNORE_CORRUPT}" == "1" ]]; then
                    echo "WARN: kmc failed on FASTA ${SAMPLE_ID}. Ignoring."
                    rm -rf "${OUTDIR}/tmp_${SAMPLE_ID}" "$kmc_input"
                    exit 0
                else
                    echo "ERROR: kmc failed." >&2; exit 1
                fi
            fi
        elif [[ "$TYPE" == "FASTQ" ]]; then
            if ! kmc -k${KMER_LENGTH} -ci${PANEL_CI} -fq -m128 -t${KMC_THREADS} @"${kmc_input}" "${OUTDIR}/${SAMPLE_ID}" "${OUTDIR}/tmp_${SAMPLE_ID}/"; then
                if [[ "${IGNORE_CORRUPT}" == "1" ]]; then
                    echo "WARN: kmc failed on FASTQ ${SAMPLE_ID}. Ignoring."
                    rm -rf "${OUTDIR}/tmp_${SAMPLE_ID}" "$kmc_input"
                    exit 0
                else
                    echo "ERROR: kmc failed." >&2; exit 1
                fi
            fi
        fi

        rm -rf "${OUTDIR}/tmp_${SAMPLE_ID}" "$kmc_input"
    fi
    echo "Finished ${SAMPLE_ID}"
fi
EOF_P1

# ---------------------------------------------------------------------------
# Phase 2: Panel Union & Metadata  (with sanity check)
# ---------------------------------------------------------------------------
cat << 'EOF' > "${SCRIPTS_DIR}/02_union_panel.slurm"
#!/bin/bash -e
#SBATCH -p jic-medium,nbi-medium
#SBATCH -t 12:00:00
#SBATCH --mem=256G
#SBATCH -c 16
#SBATCH -J teff_union
#SBATCH -o slurm-%j.out

module load kmc || echo "WARN: module load failed"
set -euo pipefail

PANEL_DBS="$1"
UNION_DIR="$2"
PANEL_SHEET="${3:-}"
PANEL_UNION_DB="${UNION_DIR}/unified_teff_kmers"
EXISTING_SIZES="${UNION_DIR}/accession_sizes.tsv"
HEADER_FILE="${UNION_DIR}/unified_list.txt"
PANEL_U_FILE="${UNION_DIR}/panel_union_size.txt"
EXCLUDED_FILE="${UNION_DIR}/panel_excluded.tsv"
MISSING_FILE="${UNION_DIR}/panel_missing.tsv"

# Inherited from the orchestrator via --export. When 1, accessions that
# Phase 1 deliberately skipped (corrupt input) are tolerated here; only
# UNEXPECTED gaps (e.g. a NODE_FAIL that --requeue couldn't recover) are fatal.
IGNORE_CORRUPT="${IGNORE_CORRUPT:-0}"

# Configurable sanity floor — accessions with fewer k-mers than this
# are excluded from the panel.  Defaults to 1e6 (any plausible accession
# should have orders of magnitude more than this for k=31).
MIN_PANEL_KMERS="${MIN_PANEL_KMERS:-1000000}"

# v6: variable-fraction mask controls and outputs.
MASK_CORE="${MASK_CORE:-1}"
CORE_FREQ_FRAC="${CORE_FREQ_FRAC:-0.95}"
VAR_DB="${UNION_DIR}/panel_variable"
EXISTING_SIZES_VAR="${UNION_DIR}/accession_sizes_var.tsv"

# ---------------------------------------------------------------------------
# Completeness assert (runs because Phase 2 depends on Phase 1 via afterany).
# Every accession listed in the panel sheet must have a KMC db on disk.
# A gap here means Phase 1 never finished that task — most commonly a
# NODE_FAIL that --requeue couldn't recover. Report it loudly rather than
# unioning a silently-incomplete panel.
# ---------------------------------------------------------------------------
if [[ -n "${PANEL_SHEET}" && -e "${PANEL_SHEET}" ]]; then
    echo "0. Verifying panel completeness against sheet: ${PANEL_SHEET}"
    echo -e "sample_id\trow\treason" > "$MISSING_FILE"
    n_expected=0
    n_missing=0
    row=0
    while IFS=$'\t' read -r sid _type _r1 _r2; do
        row=$((row+1))
        [[ -z "${sid}" ]] && continue
        n_expected=$((n_expected+1))
        if [[ ! -e "${PANEL_DBS}/${sid}.kmc_pre" || ! -e "${PANEL_DBS}/${sid}.kmc_suf" ]]; then
            echo -e "${sid}\t${row}\tmissing_kmc_db" >> "$MISSING_FILE"
            n_missing=$((n_missing+1))
        fi
    done < <(tail -n +2 "${PANEL_SHEET}")

    if (( n_missing > 0 )); then
        missing_rows=$(tail -n +2 "$MISSING_FILE" | cut -f2 | paste -sd, -)
        echo "    ${n_missing} of ${n_expected} expected panel KMC database(s) are MISSING from ${PANEL_DBS}:" >&2
        tail -n +2 "$MISSING_FILE" | awk -F'\t' '{printf "         - %s (sheet row %s)\n", $1, $2}' >&2
        echo "       Logged to: ${MISSING_FILE}" >&2
        echo "" >&2
        echo "       These accessions never built — usually a NODE_FAIL that --requeue could not recover." >&2
        echo "       Rebuild just those rows, then re-run the union:" >&2
        echo "         sbatch --requeue --array=${missing_rows} \\" >&2
        echo "             <scripts_dir>/01_build_panel.slurm ${PANEL_SHEET} ${PANEL_DBS}" >&2
        echo "         sbatch --dependency=afterok:<rebuild_jobid> \\" >&2
        echo "             <scripts_dir>/02_union_panel.slurm ${PANEL_DBS} ${UNION_DIR} ${PANEL_SHEET}" >&2
        echo "       (<scripts_dir> is the 'scripts' folder in your project workspace.)" >&2
        if [[ "${IGNORE_CORRUPT}" == "1" ]]; then
            echo "" >&2
            echo "    --ignore-corrupt is set: proceeding with the ${n_expected} - ${n_missing} accession(s) that DID build." >&2
            echo "    (Use this only if the gaps are genuinely corrupt inputs you intend to drop.)" >&2
        else
            echo "" >&2
            echo "    Refusing to build an incomplete panel union. (Re-run with --ignore-corrupt to" >&2
            echo "    proceed anyway with only the accessions that built.)" >&2
            exit 1
        fi
    else
        echo "    All ${n_expected} expected panel database(s) present."
        rm -f "$MISSING_FILE"
    fi
else
    echo "0. WARN: no panel sheet passed — skipping completeness check; will union" >&2
    echo "        whatever .kmc_pre files are present in ${PANEL_DBS}." >&2
fi

n_dbs=$(ls "${PANEL_DBS}"/*.kmc_pre 2>/dev/null | wc -l)
[[ $n_dbs -gt 0 ]] || { echo "ERROR: No panel KMC DBs generated." >&2; exit 1; }

echo "1. Counting k-mers per panel accession and applying sanity floor (>= ${MIN_PANEL_KMERS} k-mers)..."
echo -e "sample\tn_kmers" > "$EXISTING_SIZES"
echo -e "sample\tn_kmers\treason" > "$EXCLUDED_FILE"
> "$HEADER_FILE"

idx=1
config="${UNION_DIR}/kmc_union.config"
echo "INPUT:" > "$config"
declare -a sets

n_good=0
n_bad=0

for db_pre in "${PANEL_DBS}"/*.kmc_pre; do
    [ -e "$db_pre" ] || continue
    db_base="${db_pre%.kmc_pre}"
    sid=$(basename "$db_base")

    count=$(kmc_tools transform "$db_base" dump /dev/stdout 2>/dev/null | wc -l)

    if (( count < MIN_PANEL_KMERS )); then
        echo "    WARN: ${sid} has only ${count} k-mers (< ${MIN_PANEL_KMERS}) — EXCLUDING"
        echo -e "${sid}\t${count}\tbelow_min_kmers" >> "$EXCLUDED_FILE"
        n_bad=$((n_bad+1))
        continue
    fi

    echo "$sid" >> "$HEADER_FILE"
    echo -e "${sid}\t${count}" >> "$EXISTING_SIZES"

    echo "  set${idx} = ${db_base}" >> "$config"
    sets+=("set${idx}")
    idx=$((idx+1))
    n_good=$((n_good+1))
done

echo "    Panel composition: ${n_good} retained, ${n_bad} excluded"
[[ $n_good -gt 0 ]] || { echo "ERROR: no panel accessions passed sanity check." >&2; exit 1; }

if (( n_bad > 0 )); then
    echo "    Excluded accessions logged to: ${EXCLUDED_FILE}"
fi

echo "OUTPUT:" >> "$config"
expr=$(IFS=+ ; echo "${sets[*]}")
echo "  ${PANEL_UNION_DB} = ${expr}" >> "$config"

echo "2. Building KMC union from ${n_good} panel accession(s)..."
kmc_tools complex "$config" || { echo "ERROR: complex union failed." >&2; exit 1; }

echo "3. Counting union |U|..."
U=$(kmc_tools transform "${PANEL_UNION_DB}" dump /dev/stdout 2>/dev/null | wc -l)
echo "$U" > "$PANEL_U_FILE"
echo "    |U| = ${U} k-mers"

# ---------------------------------------------------------------------------
# v6: variable-fraction mask. Build a per-k-mer panel occurrence count (each
# accession contributes presence=1, summed across accessions), then keep only
# k-mers present in FEWER than CORE_FREQ_FRAC of accessions. Conserved core
# k-mers (in ~all accessions) carry no discriminating signal and dominate
# whole-set Jaccard at low coverage; Phase 3's *_var metrics are computed
# against this variable set. Non-fatal on failure: warns, leaves VAR_DB absent,
# and Phase 3 then reports legacy metrics only.
# ---------------------------------------------------------------------------
if [[ "${MASK_CORE}" == "1" ]]; then
    # MASK_PAR = how many kmc_tools transforms/intersects to run at once. The
    # per-accession set_counts and na_var steps are independent and were the
    # cause of the v6.0 Phase-2 wall-clock blowout when run serially (~2x220
    # full passes over the panel). They are now fanned out with xargs -P.
    # set_counts/intersect are streaming + low-memory, so oversubscribing cores
    # is cheap; default 12 leaves headroom on a 16-core node.
    MASK_PAR="${MASK_PAR:-12}"
    echo "4. Building variable-fraction mask (core_freq_frac=${CORE_FREQ_FRAC}, parallel=${MASK_PAR})..."
    core_min=$(awk -v n="${n_good}" -v f="${CORE_FREQ_FRAC}" 'BEGIN{v=int(n*f); if(v<1)v=1; print v}')
    echo "    core threshold: present in >= ${core_min} of ${n_good} accession(s) is masked out"
    bin_dir="${UNION_DIR}/bin_presence"
    freq_cfg="${UNION_DIR}/kmc_freq.config"
    freq_db="${UNION_DIR}/panel_freq"

    # IMPORTANT: this build runs as a STANDALONE subshell whose exit status we
    # capture explicitly (mask_rc), NOT as an `if (...)` condition. Bash
    # suppresses `set -e` inside a subshell used as an if/while/&& condition, so
    # a mid-build failure (e.g. kmc_tools complex hitting an operand/counter
    # cap) would otherwise be ignored and produce a corrupt mask. Every critical
    # command therefore carries an explicit `|| exit 1`, and we re-check the
    # VAR_DB exists at the end before declaring success.
    (
        set -e
        mkdir -p "${bin_dir}" || exit 1

        # --- 4a. Binarise each accession (count -> presence=1), in parallel ---
        binarise_one() {
            local sid="$1"
            kmc_tools transform "${PANEL_DBS}/${sid}" set_counts 1 "${bin_dir}/${sid}" >/dev/null 2>&1
        }
        export -f binarise_one
        export PANEL_DBS bin_dir
        echo "    [4a] binarising ${n_good} accessions (set_counts 1)..."
        grep -v '^[[:space:]]*$' "${HEADER_FILE}" \
            | xargs -P "${MASK_PAR}" -I{} bash -c 'binarise_one "$@"' _ {} \
            || { echo "    ERROR: a set_counts worker failed" >&2; exit 1; }

        # Assert every expected binarised db exists before summing.
        while read -r sid; do
            [ -z "$sid" ] && continue
            [ -e "${bin_dir}/${sid}.kmc_pre" ] \
                || { echo "    ERROR: binarised db missing for ${sid}" >&2; exit 1; }
        done < "${HEADER_FILE}"

        # --- 4b. Build the complex-sum config in panel order ---
        echo "INPUT:" > "${freq_cfg}"
        i=1
        declare -a bsets
        while read -r sid; do
            [ -z "$sid" ] && continue
            echo "  s${i} = ${bin_dir}/${sid}" >> "${freq_cfg}"
            bsets+=("s${i}")
            i=$((i+1))
        done < "${HEADER_FILE}"
        echo "OUTPUT:" >> "${freq_cfg}"
        sumexpr=$(IFS=+; echo "${bsets[*]}")
        echo "  ${freq_db} = ${sumexpr}" >> "${freq_cfg}"
        echo "OUTPUT_PARAMS:" >> "${freq_cfg}"
        echo "  -ci1 -cs$((n_good + 1))" >> "${freq_cfg}"

        # --- 4c. Sum (single, unparallelizable step) and reduce to var set ---
        echo "    [4c] summing occurrence counts across accessions (kmc_tools complex)..."
        kmc_tools complex "${freq_cfg}" \
            || { echo "    ERROR: kmc_tools complex failed" >&2; exit 1; }
        # variable set = occurrence count in [1, core_min-1]
        kmc_tools transform "${freq_db}" -ci1 -cx$((core_min - 1)) reduce "${VAR_DB}" \
            || { echo "    ERROR: kmc_tools reduce failed" >&2; exit 1; }
        # Positive confirmation the mask actually materialised.
        [ -e "${VAR_DB}.kmc_pre" ] \
            || { echo "    ERROR: VAR_DB not produced" >&2; exit 1; }
    )
    mask_rc=$?

    if [[ ${mask_rc} -eq 0 ]]; then
        nvar=$(kmc_tools transform "${VAR_DB}" dump /dev/stdout 2>/dev/null | wc -l)
        echo "    |variable set| = ${nvar} k-mers"

        # --- 5. Per-accession variable k-mer counts (na_var), in parallel ---
        # Each worker writes its own <sid>.cnt file; results are collated in
        # panel order afterwards so concurrent writes can't interleave.
        echo "5. Counting per-accession variable k-mers (na_var, parallel=${MASK_PAR})..."
        navar_dir="${UNION_DIR}/tmp_navar"
        mkdir -p "${navar_dir}"
        navar_one() {
            local sid="$1"
            local tmpv="${UNION_DIR}/tmp_var_${sid}"
            local cvar=0
            if kmc_tools simple "${PANEL_DBS}/${sid}" "${VAR_DB}" intersect "${tmpv}" >/dev/null 2>&1; then
                cvar=$(kmc_tools transform "${tmpv}" dump /dev/stdout 2>/dev/null | wc -l)
                rm -f "${tmpv}.kmc_pre" "${tmpv}.kmc_suf"
            fi
            printf '%s\t%s\n' "${sid}" "${cvar}" > "${navar_dir}/${sid}.cnt"
        }
        export -f navar_one
        export PANEL_DBS VAR_DB UNION_DIR navar_dir
        grep -v '^[[:space:]]*$' "${HEADER_FILE}" \
            | xargs -P "${MASK_PAR}" -I{} bash -c 'navar_one "$@"' _ {}

        echo -e "sample\tn_kmers_var" > "${EXISTING_SIZES_VAR}"
        while read -r sid; do
            [ -z "$sid" ] && continue
            if [ -e "${navar_dir}/${sid}.cnt" ]; then
                cat "${navar_dir}/${sid}.cnt"
            else
                printf '%s\t0\n' "${sid}"
            fi
        done < "${HEADER_FILE}" >> "${EXISTING_SIZES_VAR}"

        rm -rf "${bin_dir}" "${navar_dir}" "${freq_cfg}" \
               "${freq_db}.kmc_pre" "${freq_db}.kmc_suf"
        echo "    variable mask ready: ${VAR_DB} (na_var -> ${EXISTING_SIZES_VAR})"
    else
        echo "    WARN: variable-mask build failed — Phase 3 will use legacy metrics only." >&2
        rm -rf "${bin_dir}" "${UNION_DIR}/tmp_navar" 2>/dev/null || true
        rm -f "${freq_cfg}" "${freq_db}.kmc_pre" "${freq_db}.kmc_suf" 2>/dev/null || true
        rm -f "${VAR_DB}.kmc_pre" "${VAR_DB}.kmc_suf" 2>/dev/null || true
    fi
else
    echo "4. MASK_CORE=0 — skipping variable-fraction mask (legacy metrics only)."
fi

echo "Done Phase 2."
EOF

# ---------------------------------------------------------------------------
# Phase 3: Query Submission Wrapper
# ---------------------------------------------------------------------------
cat << 'EOF' > "${SCRIPTS_DIR}/03_run_queries.slurm"
#!/bin/bash -e
#SBATCH -p jic-medium,nbi-medium
#SBATCH -t 06:00:00
#SBATCH --mem=128G
#SBATCH -c 16
#SBATCH -J teff_query
#SBATCH -o slurm-%A_%a.out

module load fastp kmc || echo "WARN: module load failed"
set -euo pipefail

export EXISTING_DBS_DIR="$1"
export PANEL_UNION_DB="$2/unified_teff_kmers"
export HEADER_FILE="$2/unified_list.txt"
export EXISTING_SIZES="$2/accession_sizes.tsv"
export PANEL_U_FILE="$2/panel_union_size.txt"
export VAR_DB="$2/panel_variable"
export EXISTING_SIZES_VAR="$2/accession_sizes_var.tsv"
export WORKDIR="$3"
QUERY_SHEET="$4"
PIPELINE_SCRIPT="$5"

bash "${PIPELINE_SCRIPT}" --sample-sheet "${QUERY_SHEET}"
EOF

# ============================================================================
# Job Submission Orchestration
# ============================================================================
log "Submitting SLURM pipeline..."

J1=$(sbatch --parsable --array=1-${PANEL_COUNT} \
    --export=ALL,PANEL_CI="${PANEL_CI:-4}" \
    "${SCRIPTS_DIR}/01_build_panel.slurm" "$PANEL_SHEET" "$PANEL_DBS")
validate_job_id "$J1" "Phase 1 (panel build)"

# Phase 2 depends on Phase 1 via afterany (NOT afterok): it must run even if
# some Phase 1 tasks failed, so its completeness assert can report exactly
# which accessions are missing instead of leaving a silent stuck dependency.
# Inherits MIN_PANEL_KMERS, IGNORE_CORRUPT, and (v6) the core-mask knobs
# MASK_CORE / CORE_FREQ_FRAC, which drive the variable-fraction build.
J2=$(sbatch --parsable --dependency=afterany:${J1} \
    --export=ALL,MIN_PANEL_KMERS="${MIN_PANEL_KMERS:-1000000}",IGNORE_CORRUPT="${IGNORE_CORRUPT}",MASK_CORE="${MASK_CORE:-1}",CORE_FREQ_FRAC="${CORE_FREQ_FRAC:-0.95}",MASK_PAR="${MASK_PAR:-12}" \
    "${SCRIPTS_DIR}/02_union_panel.slurm" "$PANEL_DBS" "$PANEL_UNION_DIR" "$PANEL_SHEET")
validate_job_id "$J2" "Phase 2 (union)"

# Phase 3 depends on Phase 2 succeeding (afterok). If Phase 2 aborts on an
# incomplete panel, Phase 3 simply never becomes runnable.
J3=$(sbatch --parsable --dependency=afterok:${J2} --array=1-${QUERY_COUNT} \
    --export=ALL,JACCARD_THRESHOLD="${JACCARD_THRESHOLD:-0.85}",CONTAINMENT_THRESHOLD="${CONTAINMENT_THRESHOLD:-0.95}",REDUNDANCY_THRESHOLD="${REDUNDANCY_THRESHOLD:-0.97}",QUERY_CI="${QUERY_CI:-2}",MASK_CORE="${MASK_CORE:-1}" \
    "${SCRIPTS_DIR}/03_run_queries.slurm" \
    "$PANEL_DBS" "$PANEL_UNION_DIR" "$WORKDIR" "$QUERY_SHEET" \
    "${SCRIPTS_DIR}/teff_query_pipeline.slurm")
validate_job_id "$J3" "Phase 3 (queries)"

# ---------------------------------------------------------------------------
# Phase 4: Cohort aggregation (depends on Phase 3, runs even if some queries failed)
# ---------------------------------------------------------------------------
SELF_DIR="$(dirname "$(realpath "$0")")"
AGGREGATE_SCRIPT="${SELF_DIR}/aggregate.sh"

if [[ ! -e "${AGGREGATE_SCRIPT}" ]]; then
    warn "aggregate.sh not found next to fv4.sh (${AGGREGATE_SCRIPT})"
    warn "Cohort report will NOT be auto-generated. Run aggregate.sh manually after Phase 3 completes."
    J4=""
else
    # Copy aggregate.sh into the scripts dir so it travels with the workspace
    cp "${AGGREGATE_SCRIPT}" "${SCRIPTS_DIR}/aggregate.sh"
    chmod +x "${SCRIPTS_DIR}/aggregate.sh"

    J4=$(sbatch --parsable \
        --dependency=afterany:${J3} \
        --export=ALL,JACCARD_THRESHOLD="${JACCARD_THRESHOLD:-0.85}",CONTAINMENT_THRESHOLD="${CONTAINMENT_THRESHOLD:-0.95}" \
        -p jic-short,nbi-short -t 00:15:00 --mem=4G -c 1 \
        -J teff_aggregate \
        -o "${PROJ_DIR}/slurm-aggregate-%j.out" \
        --wrap="bash ${SCRIPTS_DIR}/aggregate.sh ${PROJ_DIR}")
    validate_job_id "$J4" "Phase 4 (aggregate)"
fi

echo ""
log "All jobs submitted successfully."
echo ""
printf "  %-30s %s\n" "Phase 1 (panel build):"    "job ${J1}  [${PANEL_COUNT} tasks, --requeue]"
printf "  %-30s %s\n" "Phase 2 (union + sanity):" "job ${J2}  [afterany ${J1}, asserts completeness]"
printf "  %-30s %s\n" "Phase 3 (queries):"        "job ${J3}  [${QUERY_COUNT} tasks, afterok ${J2}]"
if [[ -n "${J4}" ]]; then
    printf "  %-30s %s\n" "Phase 4 (cohort report):"  "job ${J4}  [afterany ${J3}]"
fi
echo ""
log "Monitor progress      : squeue -u \$USER"
log "If Phase 2 reports MISSING accessions (NODE_FAIL gaps), rebuild just those rows:"
log "  sbatch --requeue --array=<rows> ${SCRIPTS_DIR}/01_build_panel.slurm ${PANEL_SHEET} ${PANEL_DBS}"
log "  then re-run: sbatch --dependency=afterok:<jobid> ${SCRIPTS_DIR}/02_union_panel.slurm ${PANEL_DBS} ${PANEL_UNION_DIR} ${PANEL_SHEET}"
log "Phase 1 logs          : ${SCRIPTS_DIR}/../slurm-${J1}_*.out"
log "Phase 2 sanity report : ${PANEL_UNION_DIR}/panel_excluded.tsv  (low-kmer exclusions)"
log "Phase 2 missing report: ${PANEL_UNION_DIR}/panel_missing.tsv   (only if gaps found)"
log "Phase 3 logs          : ${WORKDIR}/<sample_id>/slurm-${J3}_*.out"
log "Per-query summaries   : ${WORKDIR}/<sample_id>/summary.{txt,tsv}"
if [[ -n "${J4}" ]]; then
    log "Cohort report         : ${PROJ_DIR}/cohort_report.txt  (after Phase 4 completes)"
    log "Cohort top hits       : ${PROJ_DIR}/cohort_top_hits.tsv"
fi
