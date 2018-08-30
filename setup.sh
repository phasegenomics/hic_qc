#!/usr/bin/sh
COMMIT_VERSION=`git rev-parse HEAD`
#sed s/COMMIT_VERSION/$COMMIT_VERSION/g < collateral/HiC_QC_report_template.md > collateral/HiC_QC_report_template_versioned.md
echo $COMMIT_VERSION > collateral/commit_id
pip install --user -r requirements.txt

