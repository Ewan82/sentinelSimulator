# !/home/db903833//dataLand01/enthoughtDistros/epd-7.2-2-rh5-x86_64/bin/python
# !/usr/bin/env python

# core python modules:
import subprocess
# 3rd party modules:
import numpy as np
# local modules:
from py_julesNML import *


class julesAllNML:
    """This class is populated by the contents
    of a module which contains templates
    of all the required JULES namelist files

    """
    def __init__(self):
        self.triffid_params_nml = triffid_params_nml
        self.triffid_params_txt = triffid_params_txt
        self.jules_surface_nml = jules_surface_nml
        self.jules_surface_txt = jules_surface_txt
        self.pft_params_nml = pft_params_nml
        self.pft_params_txt = pft_params_txt
        self.output_nml = output_nml
        self.output_txt = output_txt
        self.nveg_params_nml = nveg_params_nml
        self.nveg_params_txt = nveg_params_txt
        self.model_grid_nml = model_grid_nml
        self.model_grid_txt = model_grid_txt
        self.crop_params_nml = crop_params_nml
        self.crop_params_txt = crop_params_txt
        self.urban_nml = urban_nml
        self.urban_txt = urban_txt
        self.jules_vegetation_nml = jules_vegetation_nml
        self.jules_vegetation_txt = jules_vegetation_txt
        self.jules_soil_nml = jules_soil_nml
        self.jules_soil_txt = jules_soil_txt
        self.drive_nml = drive_nml
        self.drive_txt = drive_txt
        self.prescribed_data_nml = prescribed_data_nml
        self.prescribed_data_txt = prescribed_data_txt
        self.jules_radiation_nml = jules_radiation_nml
        self.jules_radiation_txt = jules_radiation_txt
        self.jules_rivers_nml = jules_rivers_nml
        self.jules_rivers_txt = jules_rivers_txt
        self.timesteps_nml = timesteps_nml
        self.timesteps_txt = timesteps_txt
        self.jules_surface_types_nml = jules_surface_types_nml
        self.jules_surface_types_txt = jules_surface_types_txt
        self.imogen_nml = imogen_nml
        self.imogen_txt = imogen_txt
        self.jules_hydrology_nml = jules_hydrology_nml
        self.jules_hydrology_txt = jules_hydrology_txt
        self.initial_conditions_nml = initial_conditions_nml
        self.initial_conditions_txt = initial_conditions_txt
        self.fire_nml = fire_nml
        self.fire_txt = fire_txt
        self.ancillaries_nml = ancillaries_nml
        self.ancillaries_txt = ancillaries_txt
        self.jules_snow_nml = jules_snow_nml
        self.jules_snow_txt = jules_snow_txt
        self.jules_soil_biogeochem_nml = jules_soil_biogeochem_nml
        self.jules_soil_biogeochem_txt = jules_soil_biogeochem_txt

    def writeNML(self):
        self.triffid_params_nml.write()
        self.jules_surface_nml.write()
        self.pft_params_nml.write()
        self.output_nml.write()
        self.nveg_params_nml.write()
        self.model_grid_nml.write()
        self.crop_params_nml.write()
        self.urban_nml.write()
        self.jules_vegetation_nml.write()
        self.jules_soil_nml.write()
        self.drive_nml.write()
        self.prescribed_data_nml.write()
        self.jules_radiation_nml.write()
        self.jules_rivers_nml.write()
        self.timesteps_nml.write()
        self.jules_surface_types_nml.write()
        self.imogen_nml.write()
        self.jules_hydrology_nml.write()
        self.initial_conditions_nml.write()
        self.fire_nml.write()
        self.ancillaries_nml.write()
        self.jules_snow_nml.write()
        self.jules_soil_biogeochem_nml.write()


class jules(julesAllNML):
    """Class to run JULES.

    :param jules_exe: location of JULES executable.
    :type jules_exe: str

    .. note:: You must have JULES installed on local system with a version of 4.8 or higher.

    """

    def __init__(self, jules_exe='/home/if910917/jules/models/jules4.8/build/bin/jules.exe'):
        """Class to run JULES.

        :param jules_exe: location of JULES executable.
        :type jules_exe: str

        .. note:: You must have JULES installed on local system with a version of 4.8 or higher.

        """

        julesAllNML.__init__(self)

        self.jules = jules_exe

    def runJules(self):
        """Write all NML files to disk.
        Run JULES in a subprocess.
        Check output for fatal errors.

        :return: stdout and stderr output from JULES model run.
        :rtype: str
        """

        # write all the nml files here so the
        # user doesn't have to remember to...
        self.writeNML()

        # run JULES
        cmd = []
        cmd.append(self.jules)

        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out = p.stdout.readlines()
        err = p.stderr.readlines()
        p.wait()

        # catch "fatal" errors
        for line in out:
            if len(line.split()) == 0: continue
            if line.split()[0] == "[FATAL":
                print >> sys.stderr, "*** runJules: caught fatal error in JULES run:"
                print >> sys.stderr, line,
                sys.exit()

        # catch anything in stderr
        if len(err) > 0:
            for line in err:
                print >> sys.stderr, "*** runJules: caught output on stderr in JULES run:"
                print >> sys.stderr, line,
                sys.exit()

        return out, err


def crop_run(sow_date=110, b=6.631, smwilt=0.1866, neff=5.70e-4):
    """
    Function that runs JULES with crop model turned on and given user defined parameters at Wallerfing site. Output is
    saved in folder and file specified within function.

    :param sow_date: Sow date, between 90 and 150.
    :type sow_date: int.
    :param b: Brooks-Corey exponent factor.
    :type b: float.
    :param smwilt: Soil moisture wilting point.
    :type smwilt: float.
    :param neff: Nitrogen use efficiency of crop (Vcmax).
    :type neff: float.
    :return: 'Done' to notify used JULES run has finished.
    :rtype: str
    """
    j = jules()
    # j.drive_nml.mapping['file']='path/to/your/drivers/metData'+n+'.dat'  # unnecessary here as using WFD for jules
    j.output_nml.mapping["jules_output_1_run_id"] = "'crp_g_" + str(sow_date) + "_" + str(b)[0:6] + "_" + str(smwilt)[0:7] +\
                                                    "_" + str(neff)[0:8] + "',"
    print j.output_nml.mapping["jules_output_1_run_id"]
    output_nml.mapping["jules_output_1_output_dir"] = "'./output/sensitivity_runs',"
    j.timesteps_nml.mapping["jules_time_1_main_run_start"] = " '2012-01-01 00:00:00',"
    j.timesteps_nml.mapping["jules_spinup_1_max_spinup_cycles"] = " 4"
    j.ancillaries_nml.mapping["jules_crop_props_1_const_val"] = " 510.11138916 501.136169434 " + str(sow_date)
    j.ancillaries_nml.mapping["jules_soil_props_1_const_val"] = str(b)+", 0.3967309, 0.0027729999, 0.45809999, " \
                                                              "0.3283205, "+str(smwilt)+", 1185786.0, 0.2269195, 0.17,"
    j.pft_params_nml.mapping["jules_pftparm_1_neff_io"] = "8.00e-4,8.00e-4,8.00e-4,4.00e-4,8.00e-4," + str(neff) + \
                                                          ", 8.00e-4,4.00e-4,8.00e-4,"
    j.runJules()
    return 'Done'


if __name__ == "__main__":

    # j = jules()
    # j.timesteps_nml.mapping["jules_time_1_main_run_start"] = " '2011-01-01 00:00:00',"
    # j.runJules()

    crop_run()
    for x in xrange(100):
        print x
        sow_date = int(np.random.normal(110, 0.1*110))
        b = np.random.normal(6.631, 0.1*6.631)
        smwilt = np.random.normal(0.1866, 0.1*0.1866)
        neff = np.random.normal(5.7e-4, 0.1*5.7e-4)
        crop_run(sow_date, b, smwilt, neff)
