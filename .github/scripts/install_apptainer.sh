#!/bin/bash

set -eo pipefail

wget https://raw.githubusercontent.com/apptainer/apptainer/main/tools/install-unprivileged.sh

chmod +x install-unprivileged.sh

mkdir -p /opt/shared/apptainer-1.3.2

./install-unprivileged.sh /opt/shared/apptainer

echo "/opt/apptainer/bin" >> $GITHUB_PATH