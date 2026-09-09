#!/bin/bash
# Is the cross-node exchange host-staged, and can the environment fix it?
#
#     run_ucx.sh <exe> <results_dir>
#
# The GPU traces (job 5139028) show that in rect_jacobi at 8 ranks / 4 per node
# the two ranks holding the single cross-node link move ~694 MB through HOST
# memory per steady window -- ~57 MB/step, which is that config's entire exchange
# volume -- while the six intra-node ranks move 25 MB and use cuda_ipc P2P. The
# same config spread 2 ranks/node stages NOTHING on any rank. So the cross-node
# path is not GPU-direct at 4 ranks/node, and that, not link count, is what
# separates 738 us from 150 us.
#
# This sweep is ENVIRONMENT ONLY -- no solver code, no rebuild, so it needs no
# bit-exactness gate. Each variant runs the anomaly's exact configuration and
# placement; only UCX's view of GPU memory changes.
set -uo pipefail

EXE="${1:?usage: run_ucx.sh <exe> <results_dir>}"
RES="${2:?usage: run_ucx.sh <exe> <results_dir>}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CFG="$HERE/../configs"
NSTEPS="${NSTEPS:-200}"

CONFIG="${CONFIG:-rect_jacobi}"
RANKS="${RANKS:-8}"; NODES="${NODES:-2}"

# name|env assignments
#   base       reproduces every previous run in this campaign
#   nocache    UCX_MEMTYPE_CACHE=n has been set by every submit script here since
#              the campaign began. It disables UCX's memory-type cache, so every
#              buffer's type is re-resolved per operation; it is a workaround
#              flag, and nobody has ever measured what it costs. Unsetting it is
#              the first thing to try.
#   gdr        force GPUDirect RDMA on for the IB transport.
#   both       unset the cache AND force GDR.
command -v mpirun >/dev/null || { echo "ERROR: mpirun not on PATH" >&2; exit 1; }
[ -x "$EXE" ] || { echo "ERROR: solver not executable: $EXE" >&2; exit 1; }
mkdir -p "$RES"
per_node=$(( RANKS / NODES ))
echo "=== ucx sweep start $(date '+%F %T')  $CONFIG ${RANKS}x${NODES} nsteps=$NSTEPS"

# name:env assignments[:GPUS permutation]
#
#   base      reproduces every previous run in this campaign
#   cache     UCX_MEMTYPE_CACHE=n has been set by every submit script here since
#             the campaign began -- a workaround flag nobody has measured.
#   gdr       force GPUDirect RDMA on for the IB transport.
#   gpumap    THE TOPOLOGY FIX, and no code changes. nvidia-smi topo -m on these
#             nodes puts all three HCAs on NUMA 0 with GPU0/GPU1 (NODE), while
#             GPU2/GPU3 reach every NIC only across the inter-socket link (SYS),
#             where GPUDirect RDMA does not apply. comm.f90 assigns
#             device = local_rank mod num_devices, so at 4 ranks/node the Morton
#             chain's cross-node ranks -- LOCAL 3 on the low node and LOCAL 0 on
#             the high one -- land on GPU3 and GPU0. GPU3 is the bad one, and it
#             is exactly the rank the traces caught staging 694 MB through host
#             memory. The permutation 0,2,3,1 puts local ranks 0 AND 3 (both
#             cross-node ends) on GPU0/GPU1 and leaves the purely intra-node
#             ranks 1,2 on GPU2/GPU3, where NVLink serves them equally well.
declare -a VAR_NAME VAR_ENV VAR_GPUS
add() { VAR_NAME+=("$1"); VAR_ENV+=("$2"); VAR_GPUS+=("${3:-}"); }
# VARIANTS entries are  name:ENV+ENV:gpu/gpu/gpu/gpu
# SLASHES, not commas, in the GPU list: `sbatch --export` uses the comma as its
# OWN separator, so a variant carrying "0,2,3,1" arrives chopped into GPUS=0 plus
# three junk entries -- which silently pins every local rank to one GPU and reads
# as a 4x slowdown rather than as the malformed input it is. Cost two runs.
if [ -n "${VARIANTS:-}" ]; then
    for v in $VARIANTS; do
        add "${v%%:*}" "$(echo "${v#*:}" | cut -d: -f1 | tr '+' ' ')" \
            "$(echo "$v" | cut -s -d: -f3 | tr '/' ',')"
    done
else
    add base   "UCX_MEMTYPE_CACHE=n"
    add cache  "UCX_MEMTYPE_CACHE=y"
    add gdr    "UCX_MEMTYPE_CACHE=n UCX_IB_GPU_DIRECT_RDMA=y"
    add gpumap "UCX_MEMTYPE_CACHE=n" "0,2,3,1"
    add gpumap_cache "UCX_MEMTYPE_CACHE=y" "0,2,3,1"
fi

# Per-rank GPU pinning wrapper (same idea as overheadTest/gpu_rank.sh): the
# solver always offloads to OpenMP device 0 within what CUDA_VISIBLE_DEVICES
# exposes, so restricting each rank to one physical card is how the mapping is
# chosen from outside the code.
WRAP="$RES/pin_gpu.sh"
cat > "$WRAP" <<'PIN'
#!/bin/bash
if [ -n "${GPUS:-}" ]; then
    IFS=',' read -r -a devs <<< "$GPUS"
    lr="${OMPI_COMM_WORLD_LOCAL_RANK:-${SLURM_LOCALID:-0}}"
    export CUDA_VISIBLE_DEVICES="${devs[$(( lr % ${#devs[@]} ))]}"
fi
exec "$@"
PIN
chmod +x "$WRAP"

for i in "${!VAR_NAME[@]}"; do
    name="${VAR_NAME[$i]}"; envs="${VAR_ENV[$i]}"; gpus="${VAR_GPUS[$i]}"
    run="$RES/${CONFIG}_r${RANKS}_N${NODES}_$name"
    [ -f "$run/run.log" ] && { echo "--- skip $name (done)"; continue; }
    rm -rf "$run"; mkdir -p "$run"
    cp "$CFG/$CONFIG.ini" "$run/config.ini"
    if grep -qE '^[[:space:]]*profile[[:space:]]*=' "$run/config.ini"; then
        sed -i 's/^profile *=.*/profile = true/' "$run/config.ini"
    else
        printf '\n[output]\nprofile = true\n' >> "$run/config.ini"
    fi
    sed -i -e "s/^nsteps *=.*/nsteps = $NSTEPS/" \
           -e "s/^runtime_interval *=.*/runtime_interval = $(( NSTEPS/4 ))/" "$run/config.ini"
    echo "$envs GPUS=${gpus:-<default>}" > "$run/env.txt"

    echo "=== $name : $envs  GPUS=${gpus:-default}  ($(date '+%F %T'))"
    # stdin from /dev/null: mpirun inherits and DRAINS the loop's stdin, which is
    # what silently reduced the first attempt at this sweep to a single variant.
    ( cd "$run" && env -u UCX_MEMTYPE_CACHE -u UCX_IB_GPU_DIRECT_RDMA -u GPUS \
        $envs ${gpus:+GPUS=$gpus} \
        mpirun -n "$RANKS" --map-by "ppr:${per_node}:node" --bind-to core \
        -x GPUS -x UCX_MEMTYPE_CACHE -x UCX_IB_GPU_DIRECT_RDMA \
        --display-map "$WRAP" "$EXE" config.ini > run.log 2>&1 < /dev/null )
    rc=$?
    if [ $rc -ne 0 ]; then
        echo "    FAILED (exit $rc)"; mv "$run/run.log" "$run/run.FAILED.log"; continue
    fi
    awk '/Data for node:/ {print $4}' "$run/run.log" | sort -u > "$run/hosts.txt"
    got=$(wc -l < "$run/hosts.txt")
    [ "$got" -ne "$NODES" ] && echo "PLACEMENT MISMATCH: wanted $NODES got $got" > "$run/VOID"
    printf "    hosts(%s) " "$got"
    grep -m1 "^timing: nsteps" "$run/run.log" | awk '{printf "s/step %s  ", $NF}'
    grep -m1 "^exch_timing: mpi_wait" "$run/run.log" | awk '{printf "wait/round %.1f us\n", $8/$4*1e6}'
    grep -m1 "exchange balance" "$run/run.log" | sed 's/^/      /'
    rm -f "$run"/overhead_*.h5
done
echo "=== ucx sweep done $(date '+%F %T')"
