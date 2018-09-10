#!/bin/bash
for f in *.txt ; do
name="$(echo $f | sed s'/_SMART_results.txt//')"
cat $f | sed -n '/SANT/{n;p;n;p;}' | cut -d '=' -f 2 | tr -s '\n' '\t' | sed  s"/^/${name}      /" >> coordinates.txt ; done
