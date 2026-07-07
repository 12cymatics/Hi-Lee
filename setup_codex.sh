#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$PWD}"
WORK_DIR="${WORK_DIR:-$ROOT_DIR/work_lrx}"
ENGINE_ZIP="${ENGINE_ZIP:-$ROOT_DIR/lrx_vedic_engine_v5.zip}"
ENGINE_DIR="${ENGINE_DIR:-$ROOT_DIR/lrx_vedic_engine_v5/lrx_vedic_engine_v5}"
SIM_SRC="${SIM_SRC:-$ROOT_DIR/simulator/hybrid_simulator.cpp}"
SIM_BIN="${SIM_BIN:-$WORK_DIR/hybrid_simulator}"
REPORT_PATH="${REPORT_PATH:-$WORK_DIR/simulation_report.json}"
BOOST_DIR="${BOOST_DIR:-$ROOT_DIR/third_party/boost_1_84_0}"
BOOST_TARBALL_URL="${BOOST_TARBALL_URL:-https://boostorg.jfrog.io/artifactory/main/release/1.84.0/source/boost_1_84_0.tar.gz}"
BOOST_DRIVE_ID="${BOOST_DRIVE_ID:-}"

MODE="${MODE:-all}"
QUBITS="${QUBITS:-4}"
SHOTS="${SHOTS:-2048}"

CXX="${CXX:-g++}"
CXXFLAGS="${CXXFLAGS:--O3 -std=c++17 -pthread}"

mkdir -p "$WORK_DIR"

need_cmd() { command -v "$1" >/dev/null 2>&1 || { echo "ERROR: Missing required command: $1" >&2; exit 1; }; }

need_cmd "$CXX"
need_cmd "unzip"
need_cmd "curl"
need_cmd "tar"

has_boost_headers() {
  echo '#include <boost/multiprecision/cpp_int.hpp>' | "$CXX" -x c++ - -fsyntax-only >/dev/null 2>&1
}

if ! has_boost_headers; then
  if [[ ! -d "$BOOST_DIR/boost" ]]; then
    echo "[0/3] Boost headers not found; downloading to $BOOST_DIR"
    mkdir -p "$ROOT_DIR/third_party"
    tmp_archive="$ROOT_DIR/third_party/boost_1_84_0.tar.gz"
    if [[ -n "$BOOST_DRIVE_ID" ]]; then
      need_cmd "python3"
      python3 -m pip -q install --user gdown
      python3 -m gdown --id "$BOOST_DRIVE_ID" -O "$tmp_archive" --quiet
    else
      curl -L "$BOOST_TARBALL_URL" -o "$tmp_archive"
    fi
    tar -xzf "$tmp_archive" -C "$ROOT_DIR/third_party"
  fi
fi

BOOST_INCLUDE_FLAGS=()
if [[ -d "$BOOST_DIR/boost" ]]; then
  BOOST_INCLUDE_FLAGS=(-I"$BOOST_DIR")
elif ! has_boost_headers; then
  echo "ERROR: Boost headers are required and could not be found or downloaded." >&2
  exit 1
fi

if [[ ! -f "$ENGINE_ZIP" ]]; then
  echo "ERROR: Engine zip not found at $ENGINE_ZIP" >&2
  exit 1
fi

if [[ ! -d "$ENGINE_DIR" ]]; then
  echo "[1/3] Unpacking engine to $ROOT_DIR/lrx_vedic_engine_v5"
  rm -rf "$ROOT_DIR/lrx_vedic_engine_v5"
  mkdir -p "$ROOT_DIR/lrx_vedic_engine_v5"
  unzip -q "$ENGINE_ZIP" -d "$ROOT_DIR/lrx_vedic_engine_v5"
fi

if [[ ! -f "$SIM_SRC" ]]; then
  echo "ERROR: Simulator source not found at $SIM_SRC" >&2
  exit 1
fi

INCLUDE_DIR="$ENGINE_DIR/references"
if [[ ! -f "$INCLUDE_DIR/vedic_sutras_complete.hpp" ]]; then
  echo "ERROR: Missing Vedic sutra header at $INCLUDE_DIR/vedic_sutras_complete.hpp" >&2
  exit 1
fi

echo "[2/3] Building hybrid simulator"
"$CXX" $CXXFLAGS -I"$INCLUDE_DIR" "${BOOST_INCLUDE_FLAGS[@]}" "$SIM_SRC" -o "$SIM_BIN"

echo "[3/3] Running hybrid simulator"
"$SIM_BIN" --mode "$MODE" --qubits "$QUBITS" --shots "$SHOTS" --report "$REPORT_PATH"

echo "DONE ✅"
echo "Report: $REPORT_PATH"
