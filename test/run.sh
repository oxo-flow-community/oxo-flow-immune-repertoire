#!/usr/bin/env bash
# Acceptance test: validate + lint + dry-run for every workflow the repository
# ships. Exits 0 only when the engine accepts all of them.
#
#   ./test/run.sh              # uses $OXO or `oxo-flow` from PATH
#
# qc_metrics.oxoflow is data-gated (it reads results/) and is validated only.
set -euo pipefail
cd "$(dirname "$0")/.."
OXO=${OXO:-oxo-flow}

echo "==> validate main.oxoflow"
"$OXO" validate main.oxoflow

echo "==> lint main.oxoflow (errors fail, warnings do not)"
"$OXO" lint main.oxoflow

echo "==> dry-run main.oxoflow"
"$OXO" dry-run main.oxoflow > /tmp/oxo-dryrun-main-$$.txt 2>&1
grep -qE "Summary:|would (execute|run)" /tmp/oxo-dryrun-main-$$.txt

echo "==> validate + dry-run shard.oxoflow (example 2-sample manifests)"
"$OXO" validate shard.oxoflow
"$OXO" dry-run shard.oxoflow > /tmp/oxo-dryrun-shard-$$.txt 2>&1
grep -qE "Summary:|would (execute|run)" /tmp/oxo-dryrun-shard-$$.txt

echo "==> validate setup.oxoflow and qc_metrics.oxoflow"
"$OXO" validate setup.oxoflow
"$OXO" validate qc_metrics.oxoflow

echo "==> all workflows accepted"
