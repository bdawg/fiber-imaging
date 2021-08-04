"""
Script for modifying RSoft .ind files
"""

import shutil
import numpy as np

datadir = '/Users/bnorris/DataAnalysis/fiber-imaging/'
# infilename = '19CorePL_July2021_noOuter_SMtoMM_extraMM_NOLAUNCHTEMPLATE.ind'
# outfilename = '19CorePL_July2021_noOuter_SMtoMM_extraMM_randlaunches04.ind'

infilename = '19CorePL_July2021_noOuter_MMtoSM_extraMM_NOLAUNCHTEMPLATE.ind'
outfilename = '19CorePL_July2021_noOuter_MMtoSM_extraMM_randlaunches04.ind'


# Add more launch modes
nlaunches = 19
power_range = (0.0, 1.0)

LP_modes = np.array([[0,1],
                     [0,2],
                     [0,3],
                     [1,1],
                     [-1,1],
                     [1,2],
                     [-1,2],
                     [2,1],
                     [-2,1],
                     [2,2],
                     [-2,2],
                     [3,1],
                     [-3,1],
                     [3,2],
                     [-3,2],
                     [4,1],
                     [-4,1],
                     [5,1],
                     [-5,1]
                     ])

launch_pathway = 20
for k in range(nlaunches):
    n = k + 1
    power = np.random.uniform(low=power_range[0], high=power_range[1])
    phase = np.random.uniform(low=0, high=360)

    # For the first iteration, read the in file and do find-and-replace on existing
    # launch details, and append field. For subsequent iterations, just append.
    if k==0:
        with open(datadir+infilename, 'r') as infile:
            infiledata = infile.read()
            infiledata = infiledata.replace('launch_power = 9', 'launch_power = %f' % power)
            infiledata = infiledata.replace('launch_phase = 9', 'launch_phase = %f' % phase)
        with open(datadir+outfilename, 'w') as outfile:
            outfile.write(infiledata)

    out_text = []
    out_text.append('launch_field %d' % n)
    out_text.append('   launch_pathway = %d' % launch_pathway)
    out_text.append('   launch_type = LAUNCH_WGMODE')
    out_text.append('   launch_tilt = 1')
    out_text.append('   launch_mode = %d' % LP_modes[k,0])
    out_text.append('   launch_mode_radial = %d' % LP_modes[k,1])
    out_text.append('   launch_normalization = 2')
    out_text.append('   launch_align_file = 1')
    out_text.append('   launch_power = %f' % power)
    out_text.append('   launch_phase = %f' % phase)
    out_text.append('end launch_field')
    out_text.append(' ')

    with open(datadir+outfilename, 'a') as file:
        for out_line in out_text:
            print(out_line)
            file.write(out_line + '\n')



# # Add more launch waveguides
# nlaunches = 19
# power_range = (0.0, 1.0)
#
# for k in range(nlaunches):
#     n = k + 1
#     power = np.random.uniform(low=power_range[0], high=power_range[1])
#     phase = np.random.uniform(low=0, high=360)
#
#     # For the first iteration, read the in file and do find-and-replace on existing
#     # launch details, and append field. For subsequent iterations, just append.
#     if k==0:
#         with open(datadir+infilename, 'r') as infile:
#             infiledata = infile.read()
#             infiledata = infiledata.replace('launch_power = 9', 'launch_power = %f' % power)
#             infiledata = infiledata.replace('launch_phase = 9', 'launch_phase = %f' % phase)
#         with open(datadir+outfilename, 'w') as outfile:
#             outfile.write(infiledata)
#
#     out_text = []
#     out_text.append('launch_field %d' % n)
#     out_text.append('   launch_pathway = %d' % n)
#     out_text.append('   launch_type = LAUNCH_WGMODE')
#     out_text.append('   launch_tilt = 1')
#     out_text.append('   launch_align_file = 1')
#     out_text.append('   launch_power = %f' % power)
#     out_text.append('   launch_phase = %f' % phase)
#     out_text.append('end launch_field')
#     out_text.append(' ')
#
#     with open(datadir+outfilename, 'a') as file:
#         for out_line in out_text:
#             print(out_line)
#             file.write(out_line + '\n')









