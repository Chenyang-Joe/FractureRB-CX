#!/bin/bash

# 运行 FractureRB_vis（注意路径和参数）
../build/FractureRB_vis \
  -i _out/bowl_ \
  -n 100 \
  -o _hi \
  -q 0.001 \
  --vis-obj \
  --vdb-rebuild Arma1 bunny bowl
