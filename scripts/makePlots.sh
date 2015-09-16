#!/bin/bash

python plotAcc.py $1
python plotMigration.py $1
python plotNminus1.py $1
python plotCumulativeEff.py $1
python plotIdEff.py $1
python plotIsoEff.py $1
python plotTrigEff.py $1
python plotIdIsoTrigEff.py $1
