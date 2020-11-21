#!/usr/bin/env python3


import nibabel as nib
import numpy as np

import sys


def main(label_file, results_file, x0_outfile, sig_outfile):
    label_nii = nib.load(label_file)
    label_data = label_nii.get_data()

    x0_data = np.zeros(label_nii.shape)
    x0_data[:, :] = np.nan
    sig_data = np.zeros(label_nii.shape)
    sig_data[:, :] = np.nan


    results = np.genfromtxt(results_file, delimiter=",", dtype=None)

    for region, mean, sig in results:
        idxs = label_data[:, :, :] == region + 1
        x0_data[idxs] = mean
        sig_data[idxs] = sig

    x0_nii = nib.Nifti1Image(x0_data, label_nii.affine)
    nib.save(x0_nii, x0_outfile)

    sig_nii = nib.Nifti1Image(sig_data, label_nii.affine)
    nib.save(sig_nii, sig_outfile)


if __name__ == '__main__':
    if len(sys.argv) != 5:
        print("Usage: ")
        print("./gen_x0 LABEL_NIFTI_FILE RESULTS_FILE OUTFILE_X0 OUTFILE_SIGMA")
        sys.exit()

    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
