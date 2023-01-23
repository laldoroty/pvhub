import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
import pandas as pd
import os

c = 299792.458
abspath = os.path.dirname(__file__)
data_dir = os.path.join(abspath, 'data')

class pv_object(object):
    """
    LNA 20230122
    Modified version of Erik Peterson's peculiar velocity
    corrections code (arXiv:2110.03487) to be object-oriented.

    """

    def __init__(self):
        self.model_spec, self.map_inside, self.map_outside = None, None, None

    def choose_model(self,modelflag):
        """Choose and load one of the four models:
        0: 2M++_SDSS (Said et al. 2020, Peterson et al. 2021, Carr et al. 2021)
        1: 2M++_SDSS_6dF (Said et al. 2020)
        2: 2MRS (Lilow & Nusser 2021)
        3: 2M++ (Carrick et al. 2015)"""
        if modelflag == 0:
            modelname = "2M++_SDSS"
            model = "2MPP_SDSS.txt"
            model_ext = "2MPP_SDSS_out.txt"
        elif modelflag == 1:
            modelname = "2M++_SDSS_6dF"
            model = "2MPP_SDSS_6dF.txt"
            model_ext = "2MPP_SDSS_6dF_out.txt"
        elif modelflag == 2:
            modelname = "2MRS"
            model = "2MRS_redshift.txt"
            model_ext = "2MRS_redshift_out.txt"
        elif modelflag == 3:
            modelname = "2M++"
            model = "2MPP_redshift.txt"
            model_ext = "2MPP_redshift_out.txt"
        else:
            raise ValueError("Unknown Model")

        print(f"Loading model {modelflag} ({modelname}).")

        # Model inside reconstruction
        self.map_inside = pd.read_csv(os.path.join(data_dir,model))

        # Model beyond reconstruction
        self.map_outside = pd.read_csv(os.path.join(data_dir,model_ext),delim_whitespace=True)

        self.model_spec = modelflag

        return model, model_ext

    def calculate_pv(self, RA, DEC, z_cmb_in, extrapolation=True):
        """Get peculiar velocities from a peculiar velocity map."""

        if self.model_spec is None:
            print("No model specified with self.choose_model(); will load default model.")
            self.choose_model(0)
        else:
            print(f"Using model {self.model_spec}.")

        vproj = self.map_inside["vproj_2MPP"]

        zcmb_m = self.map_outside["z"]
        vx = self.map_outside["Vsgx"]
        vy = self.map_outside["Vsgy"]
        vz = self.map_outside["Vsgz"]

        dmin = -20000.0
        dmax = 20000.0
        nbins = 129
        bsz = (dmax - dmin) / float(nbins - 1.0)

        cz = c * np.array(z_cmb_in)
        zcmb = z_cmb_in
        # Note that cz is *not* a distance, but treating it as such is self-consistent with 2M++
        ccc = SkyCoord(
            RA * u.degree, DEC * u.degree, distance=cz * u.km / u.s, frame="icrs"
        )
        sgc = ccc.transform_to("supergalactic")
        sgc.representation_type = "cartesian"
        xbin = np.round(((sgc.sgx.value - dmin) / bsz), 0)
        ybin = np.round(((sgc.sgy.value - dmin) / bsz), 0)
        zbin = np.round(((sgc.sgz.value - dmin) / bsz), 0)
        binindex = (
            xbin.astype(int) * nbins * nbins + ybin.astype(int) * nbins + zbin.astype(int)
        )  # calculate bin index even if coords outside 2M++

        try:
            binindex[
                np.where((binindex < 0) | (binindex >= len(vproj)))
            ] = 0  # set indices outside 2M++ to 0
        except TypeError:  # For single input
            pass

        k = np.searchsorted(zcmb_m, zcmb)  # calculate bin index even if coords inside 2M++

        in2MPP = (
            (cz < 19942)  # precise 2M++ boundary
            & ((dmin < sgc.sgx.value) & (sgc.sgx.value < dmax))
            & ((dmin < sgc.sgy.value) & (sgc.sgy.value < dmax))
            & ((dmin < sgc.sgz.value) & (sgc.sgz.value < dmax))
        )
        if extrapolation:
            r = np.sqrt(
                np.sum(np.square([sgc.sgx.value, sgc.sgy.value, sgc.sgz.value]), axis=0)
            )
            vdot = (
                (sgc.sgx.value * vx[k]) + (sgc.sgy.value * vy[k]) + (sgc.sgz.value * vz[k])
            )
            pv = np.where(
                in2MPP,
                vproj.loc[binindex],
                vdot / r,
            )
            # 2M++ is not completely spherical, so extrapolation must be extended inwards in some
            # parts of the sky. The PV of the centre of the reconstruction is undefined, so we only
            # use the extrapolation above a redshift securely inside 2M++ but outside the central cell
            pv = np.where((np.isnan(pv)) & (zcmb > 0.01), vdot / r, pv)
            pv = np.round(pv, 0)
        else:
            pv = np.where(in2MPP, np.round(vproj.loc[binindex], 0), np.nan)
        return pv
