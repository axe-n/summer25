#!/usr/bin/env bash

set -euo pipefail

PROJECT_NAME="blackray-deformed-kerr"
PROJECT_DIR="${PWD}/${PROJECT_NAME}"

if [[ -e "${PROJECT_DIR}" ]]; then
    printf 'Error: %s already exists. Remove it or choose another location.\n' "${PROJECT_DIR}" >&2
    exit 1
fi

command -v git >/dev/null 2>&1 || {
    printf 'Error: git is required but was not found in PATH.\n' >&2
    exit 1
}

mkdir -p "${PROJECT_DIR}"
cd "${PROJECT_DIR}"

git init

git clone https://github.com/ABHModels/blackray _legacy_blackray
git clone https://github.com/ABHModels/raytransfer _legacy_raytransfer

mkdir -p \
    src \
    include \
    python \
    tests \
    notebooks \
    data \
    docs

touch CMakeLists.txt Makefile README.md requirements.txt .gitignore

printf '\nCreated project at %s\n' "${PROJECT_DIR}"
printf 'Next steps:\n'
printf '  cd %q\n' "${PROJECT_DIR}"
printf '  git status\n'