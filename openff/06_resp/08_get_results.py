#!/usr/bin/env python

import pathlib
import glob

paths = glob.glob("08_results/*/05*.csv")
molnames = []

for p in paths:
    path = pathlib.Path(p)
    conf10 = path.parent / "01_conf-10.json"
    if conf10.exists():
        molnames.append(f'"hpc3:/pub/lilyw7/pydev/scrapbook/openff/06_resp/{path.parent}"')

with open("copier.sh", "w") as f:

    f.write(f"""
#!/bin/bash

for path in {' '.join(molnames)} ; do scp -r $path . ; done
    """)