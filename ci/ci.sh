#!/bin/bash
set -e
set -o pipefail

OUT_DIR=ci

case "$1" in
  run)
    shift
    echo "Running test pipeline..." >&2
    NXF_VER=22.04.0 nextflow run . -resume -with-docker --dir ${OUT_DIR} $@
    md5sum ${OUT_DIR}/*/*/*
    ;;
  validate)
    echo "Validating test results..." >&2
    md5sum -c ${OUT_DIR}/md5s.txt
    ;;
  cleanup)
    echo "Cleaning up test results..." >&2
    find ${OUT_DIR} -type l -exec rm {} \+
    find ${OUT_DIR} -type d -empty -exec rmdir {} \+
    find ${OUT_DIR} -type d -empty -exec rmdir {} \+
    ;;
  *)
    echo "Usage: ci.sh {run|validate|cleanup}" >&2
    exit 1
esac
