#!/usr/bin/env bash
set -eu
AUDIT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
ROOT_DIR="$(cd "$AUDIT_DIR/.." && pwd)"
SOURCE_DIR="$AUDIT_DIR/nauty2_9_3"
RUNTIME_DIR="$AUDIT_DIR/runtime_293"
# Official release already extracted from archive hashed in upstream_manifest.json.
# Exact prior library build commands, executed in SOURCE_DIR:
# ./configure --enable-generic > ../configure_2_9_3.log 2>&1
# make -j2 nauty.a > ../build_2_9_3.log 2>&1
gcc -O3 -std=c11 -I "$SOURCE_DIR" -c "$ROOT_DIR/computation/nauty_bridge.c" -o "$RUNTIME_DIR/nauty_bridge.o"
g++ -O3 -std=c++17 "$ROOT_DIR/computation/graph_engine.cpp" "$RUNTIME_DIR/nauty_bridge.o" "$SOURCE_DIR/nauty.a" -o "$RUNTIME_DIR/graph_engine"
cp "$SOURCE_DIR/COPYRIGHT" "$RUNTIME_DIR/COPYRIGHT"
cp "$SOURCE_DIR/nauty.a" "$RUNTIME_DIR/nauty.a"
