#!/usr/bin/env bash
cd "$( dirname "${BASH_SOURCE[0]}" )"
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:${LD_LIBRARY_PATH:-}"
export SEL_FRAC=1
exec root -l -b -q 'ccpi_selection.C+'
