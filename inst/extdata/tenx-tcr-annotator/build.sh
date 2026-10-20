#!/bin/sh
set -eu
task_root=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
cd "$task_root"
cargo build --release --locked -p tenx-tcr-annotate
cp "${CARGO_TARGET_DIR:-target}/release/tenx-tcr-annotate" "$task_root/tenx-tcr-annotate"
printf '%s\n' "Built $task_root/tenx-tcr-annotate"
