#!/usr/bin/env sh

mrview \
  -load T1_in_bzero.nii.gz        \
  -mode 2              \
  -comments 0          \
  -voxelinfo 0         \
  -focus 0 \
  -colourbar 0 \
  -orientationlabel 0 \
  -size 1000,600 \
  -overlay.load x0.nii.gz \
  -overlay.interpolation_off \
  -overlay.colourmap 1 \
  -capture.folder Fig/ \
  -capture.prefix BSM \
  -capture.grab        \
#  -exit 

