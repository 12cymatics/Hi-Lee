# Hybrid Quantum-Classical Vedic Simulator

This repo now runs entirely from local assets. The Vedic engine zip is already in the repo, and the simulator compiles against the complete 29-sutra header to drive a hybrid quantum-classical simulation across serial, concurrent, and parallel modes.

## Quick start

```bash
chmod +x setup_codex.sh
bash setup_codex.sh
```

The script will:
1. Unpack the engine zip if needed.
2. Download Boost headers (if missing) into `third_party/boost_1_84_0`.
3. Build the hybrid simulator (C++17 + Boost headers).
4. Run the simulation in all modes and emit a report.

The report is written to:
```
work_lrx/simulation_report.json
```

## Configuration

You can override settings via environment variables:

```bash
MODE=serial QUBITS=5 SHOTS=4096 REPORT_PATH=work_lrx/custom_report.json bash setup_codex.sh
```

Available options:
- `MODE`: `serial`, `concurrent`, `parallel`, or `all` (default: `all`)
- `QUBITS`: number of qubits (1–8, default: `4`)
- `SHOTS`: recorded in the report metadata (default: `2048`)
- `REPORT_PATH`: output path for the JSON report
- `BOOST_DIR`: override the Boost header directory
- `BOOST_TARBALL_URL`: override the Boost source tarball URL used for header download
- `BOOST_DRIVE_ID`: Google Drive file ID for a Boost tarball (used if external downloads are blocked)

## Direct build/run

```bash
g++ -O3 -std=c++17 -pthread \
  -I lrx_vedic_engine_v5/lrx_vedic_engine_v5/references \
  simulator/hybrid_simulator.cpp \
  -o work_lrx/hybrid_simulator

./work_lrx/hybrid_simulator --mode all --qubits 4 --shots 2048 --report work_lrx/simulation_report.json
```

## Notes

- The simulator requires Boost headers (for `boost::multiprecision` and `boost::rational`). `setup_codex.sh` will download a local header-only copy if they are not installed system-wide. If external downloads are blocked, place a Boost tarball on Google Drive and set `BOOST_DRIVE_ID` so the script can fetch it.
- The simulator uses all 29 Vedic sutras (16 main + 13 sub-sutras) to generate classical outputs that parameterize a quantum state-vector simulation.
- The serial, concurrent, and parallel execution modes are all reported to the same JSON output.
