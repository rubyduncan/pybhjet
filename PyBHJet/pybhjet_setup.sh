#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [[ -n "${PYBHJET_PYTHON:-}" ]]; then
    PYTHON_EXECUTABLE="$PYBHJET_PYTHON"
elif [[ -n "${VIRTUAL_ENV:-}" ]]; then
    PYTHON_EXECUTABLE="$VIRTUAL_ENV/bin/python"
else
    PYTHON_EXECUTABLE="$(command -v python3)"
fi

if [[ ! -x "$PYTHON_EXECUTABLE" ]]; then
    echo "Python executable not found: $PYTHON_EXECUTABLE" >&2
    exit 1
fi

if ! "$PYTHON_EXECUTABLE" -c "import pybind11"; then
    echo "pybind11 is missing from $PYTHON_EXECUTABLE" >&2
    echo "Activate the intended environment and run: python -m pip install pybind11" >&2
    exit 1
fi

BUILD_DIR="$SCRIPT_DIR/build"
BACKUP_DIR="$SCRIPT_DIR/build-python39-backup"
if [[ -d "$BUILD_DIR" ]] && find "$BUILD_DIR" -maxdepth 1 -name 'pybhjet.cpython-39-*.so' -print -quit | grep -q .; then
    if [[ -e "$BACKUP_DIR" ]]; then
        echo "Refusing to overwrite existing backup: $BACKUP_DIR" >&2
        exit 1
    fi
    mv "$BUILD_DIR" "$BACKUP_DIR"
    echo "Saved the Python 3.9 build to $BACKUP_DIR"
fi

echo "Building PyBHJet with $PYTHON_EXECUTABLE..."
cmake -S "$SCRIPT_DIR" -B "$BUILD_DIR" -DPython3_EXECUTABLE="$PYTHON_EXECUTABLE"
cmake --build "$BUILD_DIR"
echo "Build complete: $BUILD_DIR"
