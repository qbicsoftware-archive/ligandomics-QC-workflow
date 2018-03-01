#!/bin/bash
module load qbic/anaconda2/2.1.0
module load qbic/openms/2.0-44ed56b
module load qbic/netMHCpan/3.0
module load qbic/netMHCIIpan/3.1
module load qbic/comet/2015024
module load qbic/r/qbic-r-3.2.2

workflowDir=$(cat wfdir)
#parse using CTDopts and run workflow
python runWorkflow.py $workflowDir
cp wfdir wfdir2