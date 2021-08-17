import os
import numpy as np
import logging
import scipy.io as sio
import warnings
import matplotlib.pyplot as plt
import matplotlib.cbook
from abc import ABCMeta, abstractmethod
import errno


class RSOFTPSF(metaclass=ABCMeta):
    psf = None
    wf = None

    osys_obj = None
    pupil_wf = None

    def makePSF(self, makePSFInputDict, makePSFOptions):
        pass

    def saveToRSoft(self, outfile="PSFOut", size_data=400, savemat=False, normtomax=False):
        """
        From Theo's FITS2rsoft script.
        - The parameter 'size_data' is the physical size of the array in um so make
        sure you have the dimensions correct.

        - The format of the .fld file is 2 columns for each single column in the
        fits files, the real then the imaginary. There is a header at the top of the
        file which shouldn't need changing, it is just how the data is interpreted
        by the program
        """

        psf_real = self.complex_psf.real
        psf_imag = self.complex_psf.imag

        len_data = psf_real.shape[0]

        # Empty array to take all data 3d
        whole_psf = np.zeros([len_data, 2 * len_data])

        # RSoft takes data in columns of real and imaginary parts. (make odd ones imag, even real)
        whole_psf[:, ::2] = psf_real
        whole_psf[:, 1::2] = psf_imag
        if normtomax:
            whole_psf = whole_psf / whole_psf.max()

        # -------------------------- Writing FLD file ----------------------------------
        outfile = outfile + "_inputfield"
        print("Writing field to file " + outfile + ".fld")
        header = (
            "/rn,a,b/nx0\n/rn,qa,qb\n{0} -{1} {1} 0 OUTPUT_REAL_IMAG_3D\n{0} -{1} {1}"
        ).format(len_data, size_data)
        np.savetxt(outfile + ".fld", whole_psf, fmt="%.18E", header=header, comments="")

        if self.wf is not None:
            allData = []
            psf_ampl = self.wf.amplitude
            psf_phase = self.wf.phase
            pupil_phase = self.pupil_wf
            cur_wf = [psf_ampl, psf_phase, pupil_phase]
            allData.append(cur_wf)
            if savemat:
                sio.savemat(
                    outfile + ".mat",
                    mdict={
                        "psf_real": psf_real,
                        "psf_imag": psf_imag,
                        "allData": allData,
                        "complex_psf": self.complex_psf,
                    },
                )

    def makeMultiplePSFs(
        self,
        makePSFInputDictList,
        makePSFOptions,
        outpath="./",
        filePrefix="PSFs",
        extraPlots=False,
        size_data=400,
        makeBatFile=False,
        saveAllData=False,
        indFile="bptmp.ind",
        outPrefix="BPScan",
        numBatfiles=1,
        trimLeadingCoeffnames=None,
    ):
        allOutfiles = []
        allData = []
        try:
            os.makedirs(outpath)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise

        for ind, makePSFInputDict in enumerate(makePSFInputDictList):
            coeffs = makePSFInputDict["coeffs"]
            print("Making PSF with coeffs " + str(coeffs))
            self.makePSF(makePSFInputDict, makePSFOptions)
            plt.pause(0.001)
            # filestr = outpath + filePrefix
            filestr = filePrefix
            cocount = 0
            if trimLeadingCoeffnames is "seq":
                filestr = "set_" + str(ind)
            elif trimLeadingCoeffnames is None:
                for c in coeffs:
                    filestr = filestr + "_" + "%.1f" % c
            else:
                if cocount >= trimLeadingCoeffnames:
                    filestr = filestr + "_" + "%.1f" % c
                cocount = cocount + 1
            # print(filestr)
            self.saveToRSoft(outfile=outpath + filestr, size_data=size_data)
            allOutfiles.append(filestr)
            if saveAllData:
                psf_ampl = self.wf.amplitude
                psf_phase = self.wf.phase
                pupil_phase = self.pupil_wf
                cur_wf = [psf_ampl, psf_phase, pupil_phase]
                allData.append(cur_wf)

        if makeBatFile:
            allOutfilenames = []
            progname = "bsimw32"

            nf = len(allOutfiles)
            batfileLength = nf // numBatfiles
            allBatfileNames = []

            for k in range(numBatfiles):
                startInd = k * batfileLength
                endInd = (k + 1) * batfileLength
                if k == (numBatfiles - 1):
                    curOutfiles = allOutfiles[startInd:]
                else:
                    curOutfiles = allOutfiles[startInd:endInd]
                print(
                    "Making .bat file: " + outpath + outPrefix + "_" + str(k) + ".bat"
                )
                batfile = open(outpath + outPrefix + "_" + str(k) + ".bat", "w")
                allBatfileNames.append(outPrefix + "_" + str(k) + ".bat")
                for launch_file in curOutfiles:
                    cmdStr = (
                        progname
                        + " "
                        + indFile
                        + " prefix="
                        + outPrefix
                        + launch_file
                        + " launch_file="
                        + launch_file
                        + "_inputfield"
                        + ".fld wait=0\n"
                    )
                    print(cmdStr)
                    batfile.write(cmdStr)
                    allOutfilenames.append(outPrefix + launch_file)
                batfile.close()
            metadatafile = outpath + outPrefix + "_metadata"

            coeffsList = []
            for coeff in makePSFInputDictList:
                coeffsList.append(coeff["coeffs"])
            np.savez(
                metadatafile + ".npz",
                allOutfilenames=allOutfilenames,
                coeffsList=coeffsList,
                allData=allData,
            )
            sio.savemat(
                metadatafile + ".mat",
                mdict={
                    "allOutfilenames": allOutfilenames,
                    "coeffsList": coeffsList,
                    "allData": allData,
                },
            )

            if numBatfiles > 1:
                superbatfile = open(outpath + "runAllBatfiles.bat", "w")
                for batfilename in allBatfileNames:
                    cmdStr = "start cmd /k call " + batfilename + "\n"
                    superbatfile.write(cmdStr)
                superbatfile.close()
