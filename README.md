Create two files in the repo root: drive_manifest.json and setup_codex.sh, then run: bash setup_codex.sh

FILE 1: drive_manifest.json  (copy/paste exactly)
{
  "engine_zip_id": "1Qn_kPUHAIYWmZdu9WheuhowEyxwPF1K5",
  "test_csv_id": "PUT_TEST_CSV_FILE_ID_HERE",
  "engine_zip_sha256": "",
  "test_csv_sha256": ""
}

Replace PUT_TEST_CSV_FILE_ID_HERE with the Google Drive file ID of your Kaggle test.csv (from a link like https://drive.google.com/file/d/FILE_ID/view).
Your Drive files must be shareable (“Anyone with link can view”) unless you switch to OAuth/service-account.

FILE 2: setup_codex.sh  (copy/paste exactly; then chmod +x setup_codex.sh)

#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="${ROOT_DIR:-$PWD}"
WORK_DIR="${WORK_DIR:-$ROOT_DIR/work_lrx}"
MANIFEST_PATH="${MANIFEST_PATH:-$ROOT_DIR/drive_manifest.json}"

PYTHON_BIN="${PYTHON_BIN:-python3}"
VENV_DIR="${VENV_DIR:-$WORK_DIR/.venv}"

TIME_PER_PERM="${TIME_PER_PERM:-1.2}"
BEAM="${BEAM:-4096}"
STEPS="${STEPS:-28}"
BACKDEPTH="${BACKDEPTH:-7}"

OUT_SUBMISSION="${OUT_SUBMISSION:-$WORK_DIR/submission.csv}"

die() { echo "ERROR: $*" >&2; exit 1; }
need_cmd() { command -v "$1" >/dev/null 2>&1 || die "Missing required command: $1"; }

read_json_value() { local jf="$1"; local sel="$2"; jq -r "$sel" "$jf"; }
sha256_file() { sha256sum "$1" | awk '{print $1}'; }

gdrive_download_file() {
  local file_id="$1"
  local out_path="$2"
  mkdir -p "$(dirname "$out_path")"
  "$PYTHON_BIN" -m gdown --id "$file_id" -O "$out_path" --quiet
  test -s "$out_path" || die "Downloaded file is empty: $out_path"
}

need_cmd "$PYTHON_BIN"
need_cmd "jq"
need_cmd "unzip"

test -f "$MANIFEST_PATH" || die "Missing manifest: $MANIFEST_PATH"

ENGINE_ZIP_ID="$(read_json_value "$MANIFEST_PATH" '.engine_zip_id')"
TEST_CSV_ID="$(read_json_value "$MANIFEST_PATH" '.test_csv_id')"
ENGINE_ZIP_SHA256_EXPECT="$(read_json_value "$MANIFEST_PATH" '.engine_zip_sha256 // empty')"
TEST_CSV_SHA256_EXPECT="$(read_json_value "$MANIFEST_PATH" '.test_csv_sha256 // empty')"

[[ -n "$ENGINE_ZIP_ID" && "$ENGINE_ZIP_ID" != "null" ]] || die "engine_zip_id missing in manifest"
[[ -n "$TEST_CSV_ID" && "$TEST_CSV_ID" != "null" && "$TEST_CSV_ID" != "PUT_TEST_CSV_FILE_ID_HERE" ]] || die "test_csv_id missing in manifest (paste your Drive file ID for test.csv)"

mkdir -p "$WORK_DIR"

ENGINE_ZIP_PATH="$WORK_DIR/engine.zip"
TEST_CSV_PATH="$WORK_DIR/test.csv"

echo "[1/6] Create venv: $VENV_DIR"
"$PYTHON_BIN" -m venv "$VENV_DIR"
# shellcheck disable=SC1091
source "$VENV_DIR/bin/activate"

echo "[2/6] Install python deps (gdown)"
pip -q install --upgrade pip
pip -q install gdown

echo "[3/6] Download engine zip (Drive id=$ENGINE_ZIP_ID)"
gdrive_download_file "$ENGINE_ZIP_ID" "$ENGINE_ZIP_PATH"
if [[ -n "$ENGINE_ZIP_SHA256_EXPECT" && "$ENGINE_ZIP_SHA256_EXPECT" != "null" ]]; then
  got="$(sha256_file "$ENGINE_ZIP_PATH")"
  [[ "$got" == "$ENGINE_ZIP_SHA256_EXPECT" ]] || die "Engine zip SHA256 mismatch: got=$got expected=$ENGINE_ZIP_SHA256_EXPECT"
fi

echo "[4/6] Download test.csv (Drive id=$TEST_CSV_ID)"
gdrive_download_file "$TEST_CSV_ID" "$TEST_CSV_PATH"
if [[ -n "$TEST_CSV_SHA256_EXPECT" && "$TEST_CSV_SHA256_EXPECT" != "null" ]]; then
  got="$(sha256_file "$TEST_CSV_PATH")"
  [[ "$got" == "$TEST_CSV_SHA256_EXPECT" ]] || die "test.csv SHA256 mismatch: got=$got expected=$TEST_CSV_SHA256_EXPECT"
fi

echo "[5/6] Unzip engine"
rm -rf "$WORK_DIR/engine_unpack"
mkdir -p "$WORK_DIR/engine_unpack"
unzip -q "$ENGINE_ZIP_PATH" -d "$WORK_DIR/engine_unpack"

ENGINE_SUBDIR="$(find "$WORK_DIR/engine_unpack" -type f -name "lrx_solver_v5.py" -print -quit | xargs -r dirname)"
[[ -n "$ENGINE_SUBDIR" ]] || die "Could not find lrx_solver_v5.py after unzip"

SOLVER="$ENGINE_SUBDIR/lrx_solver_v5.py"
VERIFY="$ENGINE_SUBDIR/verify.py"

test -f "$SOLVER" || die "Missing solver: $SOLVER"
test -f "$VERIFY" || die "Missing verify: $VERIFY"

test -f "$ENGINE_SUBDIR/lrx_n8_table.bin" || die "Missing lrx_n8_table.bin"
test -f "$ENGINE_SUBDIR/lrx_n9_table.bin" || die "Missing lrx_n9_table.bin"
test -f "$ENGINE_SUBDIR/lrx_n8_cyclicmindist.u16" || die "Missing lrx_n8_cyclicmindist.u16"
test -f "$ENGINE_SUBDIR/lrx_n9_cyclicmindist.u16" || die "Missing lrx_n9_cyclicmindist.u16"

echo "[6/6] Run solver + verify"
"$PYTHON_BIN" "$SOLVER"   --input "$TEST_CSV_PATH"   --output "$OUT_SUBMISSION"   --time "$TIME_PER_PERM"   --beam "$BEAM"   --steps "$STEPS"   --backdepth "$BACKDEPTH"   --selfcheck

"$PYTHON_BIN" "$VERIFY" --input "$TEST_CSV_PATH" --moves "$OUT_SUBMISSION"

echo "DONE ✅"
echo "Submission: $OUT_SUBMISSION"

RUN:
chmod +x setup_codex.sh
bash setup_codex.sh

TUNING (optional):
TIME_PER_PERM=2.0 BEAM=8192 STEPS=34 BACKDEPTH=9 bash setup_codex.sh
