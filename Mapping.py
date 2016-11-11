# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
from nanoraw import correct_raw as cr
fName = "IMB14_011406_LT_20160928_FNFAB27163_MN17279_mux_scan_GN_003_R9_280916_65776_ch27_read13_strand.fast5"
ref = "pacbio_ref.fa"
gmap_dir = "/home/haotianteng/UQ/BINF7000/Nanopore/graphmap/bin/Linux-x64/graphmap"
b_group = "Basecall_1D_000"
c_group = "RawGenomeCorrected_000"
foo = cr.correct_raw_data(fName, ref, gmap_dir, b_group,
                          c_group, in_place = False)

# names of the columns in the returned numpy array
col_names = np.array(["norm_mean", "norm_stdev", "start",
                        "length", "base"])

print foo