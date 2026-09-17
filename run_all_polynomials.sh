#!/usr/bin/env bash

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$script_dir"

root -l -b -q \
	-e '.L TLVUtils.cxx' \
	-e '.L AIZ.C' \
	-e 'test = false;' \
	-e 'override = true;' \
	-e 'appendPolynomialOutputs = true;' \
	-e 'FiducialCut = true;' \
	-e 'FiducialCutEtaonly = false;' \
	-e 'FiducialCutCCCF = false;' \
	-e 'FiducialCutCFonly = true;' \
	-e 'AIZ(true, 0);' \
	-e 'AIZ(true, 1);' \
	-e 'AIZ(true, 2);' \
	-e 'AIZ(true, 3);' \
	-e 'AIZ(true, 4);' \
	-e 'AIZ(true, 5);' \
	-e 'AIZ(true, 6);' \
	-e 'AIZ(true, 7);'