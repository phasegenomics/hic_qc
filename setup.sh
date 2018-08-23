#!/usr/bin/sh
COMMIT_VERSION=`git rev-parse HEAD`
sed s/COMMIT_VERSION/$COMMIT_VERSION/g < collateral/HiC_QC_report_template.md > collateral/HiC_QC_report_template_versioned.md
pip install --user -r requirements.txt

