#!/usr/bin/python
"""
Description:
created 2011-01-06 by Dennis Mackin to kick off jobs on the MD Anderson Cluster.
modified 2016-07-11 by Dennis Mackin to run jobs on a single server

"""
__author__ = "Dennis Mackin <dsmackin@mdanderson.org>"
__date__ = "2011-01-06"
__credits__ = """  """
__version__ = "$Revision: 0.1.1$"

# -------------------------------------------------------------------------------
#  IMPORT STATEMENTS
# -------------------------------------------------------------------------------
import sys, os, re, subprocess, multiprocessing, functools
import string
import collections

# -------------------------------------------------------------------------------
#  RUN PARAMETERS
# -------------------------------------------------------------------------------
WORKING_DIR = '.'
OUTPUT_FOLDER = "/scratch1/2025"
EXECUTABLE = '%s/core' % WORKING_DIR

#EVENT_FILE_PATH = '/y_drive/CCData/pj2button/centered_2x.csv'
#EVENT_FILE_PATH = '/y_drive/CCData/uct_2x/all_05_-05_11_2x.csv'
EVENT_FILE_PATH = '%s/run3.csv' % WORKING_DIR


PLOT_RECON_SCRIPT = "%s/scripts/recon_plotter.py" % WORKING_DIR
PLOT_EVENTS_SCRIPT = "%s/scripts/plot_events.py" % WORKING_DIR

WWW_FOLDER = "/home/dsmackin/public_html/2025"
DELETE_IMAGES = True

TEMPLATE_FILE = 'template.cfg'
CONFIG_FILE = 'generated.cfg'


SIMULTANEOUS_JOBS = 3 #Not sure multiple jobs can share the GPU

SCRIPT_PARAMETERS_DICT = {
    'YYY_RANDOM_SEED': 817,
    #'YYY_MONIKER': 'thor',
    'YYY_EVENT_FILE_PATH': EVENT_FILE_PATH,
    'YYY_OUTPUT_FOLDER_PATH': OUTPUT_FOLDER,

    'YYY_IMAGE_ALGORITHM': 'OCTANE',
    'YYY_DCA_CENTER_X': 5.0,
    'YYY_DCA_CENTER_Y': -5.0,
    'YYY_DCA_CENTER_Z': 11.0,

    'YYY_INVERSE_SQUARE_PARAM': 1.3,

    #POINT
    'YYY_TEMPERATURE': 0.9,
    'YYY_NUMBER_TRIES_FOR_RANDOM': '1000',
    'YYY_ITERATIONS': 700,
    'YYY_EVENT_MULTIPLIER': 1,
    'YYY_DENSITY_ESTIMATOR_TYPE': 2,
    'YYY_NUMBER_OF_SHIFTS': 50,
    'YYY_MAXIMUM_NUM_CONES': 2000000,
    'YYY_NUM_CONES_OFFSET': 0,

    'YYY_USE_PARABOLAS': 1,
    'YYY_MIN_GAMMA_ENERGY': 0.642,
    'YYY_MAX_GAMMA_ENERGY': 0.682,
    'YYY_MIN_SCATTERING_ANGLE': 0,
    'YYY_MAX_SCATTERING_ANGLE': 180,
    'YYY_DATA_FILE_FORMAT': 3,
    'YYY_SCATTER_DISTANCE': 5,
    'YYY_DCA_CUT': 50,
    'YYY_X_DCA_CUT': 5.0,
    'YYY_Y_DCA_CUT': -5.0,
    'YYY_Z_DCA_CUT': 11.0,
    'YYY_MAX_ENERGY_LOST': 100.0,
    #DCA LINE CUT
    'YYY_KNOWN_GAMMA_ENERGIES': '1.17, 1.33',

    #DCA for line is segment between BEAM_LINE_POINT1 and BEAM_LINE_POINT2
    'YYY_BEAM_LINE_POINT1': '0, 0, -200',
    'YYY_BEAM_LINE_POINT2': '0, 0, 200',
    'YYY_DCA_LINE_CUT': '15.0',
    'YYY_MIN_ENERGY_SCATTER': '1.0',
    'YYY_MIN_ENERGY_EVENT': '8.0',
    'YYY_X_MIN': -20.0,
    'YYY_X_MAX': 450,

    'YYY_Y_LENGTH': 60,
    'YYY_Y_MIN': -50,
    'YYY_Y_MAX': 50,

    'YYY_Z_MIN': -50,
    'YYY_Z_MAX': 50,
    'YYY_BINS': 100,
    
    'YYY_OUTPUT_BINS_X' : '470',
    'YYY_OUTPUT_BINS_Y' : '100',
    'YYY_OUTPUT_BINS_Z' : '100',

    'YYY_OUTPUT_X_MIN': -20,
    'YYY_OUTPUT_X_MAX': 350,

    'YYY_OUTPUT_Y_MIN': -50,
    'YYY_OUTPUT_Y_MAX': 50,

    'YYY_OUTPUT_Z_MIN': -50,
    'YYY_OUTPUT_Z_MAX': 50,

    #Kernel Backprojection Parameters
    'YYY_KERNEL_BANDWIDTH' : 8.0,
    'YYY_SYSTEM_MATRIX_SCALAR' : 4,
    
    #OCTANE
    'YYY_INTERCEPT_DCA': 2,
    'YYY_PHANTOM_CENTER_X': 0,
    'YYY_PHANTOM_CENTER_Y': 0,
    'YYY_PHANTOM_CENTER_Z': 0,
    'YYY_PHANTOM_BINS': 32,
    'YYY_SOURCE_AXIS_DISTANCE': 400,
    'YYY_CONE_LENGTH_CORRECTION': 1.03,
    'YYY_NUMBER_OF_THREADS': 5,
    'YYY_NUMBER_OF_ITERATIONS' : 2,
    'YYY_PHANTOM_LENGTH' : 100,
    'YYY_BIN_WIDTH' : 6.0,
    'YYY_OCTANE_ITERATIONS' : 1,    
}

CONFIG_FILE_TEMPLATE = "%s/%s" % (WORKING_DIR, TEMPLATE_FILE)

# -------------------------------------------------------------------------------
#    FUNCTIONS
# -------------------------------------------------------------------------------



# def setJobParameters(jobRunnerObj, params_dict):
#     for key in params_dict.keys():
#         jobRunnerObj.set_parameter('YYY_%s' % key, params_dict[key])


def setParametersInString(string_with_parameters, parameters_dict):
    ''' replaces YYY_ parameters in config file with
            values from self._EMMU_CONFIG_PARAMS '''

    for param in parameters_dict.keys():
        # print "%s => %s" % ( param, self._EMMU_CONFIG_PARAMS[param] )
        string_with_parameters = string_with_parameters.replace(param, "%s" % parameters_dict[param])

    xxxSearchString = "YYY_\w+"
    xxxRE = re.compile(xxxSearchString)
    match = xxxRE.search(string_with_parameters)
    if match:
        print "ERROR: failed to set value for %s in \n -------- \n %s\n" % (match.group(), string_with_parameters)
        sys.exit(-100)
    assert (string_with_parameters.find("YYY_") == -1)

    return string_with_parameters


def make_folder(folder_path):
    def get_parent_folder(file_path):
        parts = re.split(r'/|\\', file_path)
        if len(parts) > 1:
            return "/".join(parts[:-1])
        else:
            return file_path

    if os.path.exists(folder_path):
        return folder_path

    parent = get_parent_folder(folder_path)
    if parent != folder_path:
        make_folder(parent)

    os.makedirs(folder_path)

    return folder_path


def getJobDescription(params_dict):
    description = ""
    delimiter = ""
    keys = params_dict.keys()
    for key in keys:
        description = "%s%s%s-%s" % (description, delimiter, key, params_dict[key])
        delimiter = "__"

    return description


def run_exe(exe, config_file_path, output_folder, job_name):
    cmd = "time %s %s" % (exe, config_file_path)
    results = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).stdout.read()
    ofile = open("%s/%s.log" % (output_folder, job_name), "w")
    ofile.write(results)
    ofile.close()

    return results


def run_plot_reconstruction(script_path, output_folder, www_folder, job_name):
    cmd = "python %s %s/%s/output.dat %s/%s" % (script_path, output_folder, job_name, output_folder, job_name)
    results = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).stdout.read()
    ofile = open("%s/%s_recon.log" % (output_folder, job_name), "w")
    ofile.write(results)
    ofile.close()

    cmd = "cp %s/%s/x-21.0.png %s/x-21,0_%s.png" % (output_folder, job_name, www_folder, job_name)
    print subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).stdout.read()
    #cmd = "cp %s/%s/x-16.0.png %s/x-16.0_%s.png" % (output_folder, job_name, www_folder, job_name)
    #print subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).stdout.read()
    cmd = "cp %s/%s/y6.0.png %s/y6.0_%s.png" % (output_folder, job_name, www_folder, job_name)
    print subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).stdout.read()
    cmd = "cp %s/%s/z11.0.png %s/z11.0_%s.png" % (output_folder, job_name, www_folder, job_name)
    print subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).stdout.read()
    #cmd = "cp %s/%s/z16.0.png %s/z16.0_%s.png" % (output_folder, job_name, www_folder, job_name)
    #print subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).stdout.read()

    cmd = "cp %s/%s/profile_z_0.0_0.0.png %s/profile_z_0.0_0.0_%s.png" % (output_folder, job_name, www_folder, job_name)
    print subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).stdout.read()
    cmd = "cp %s/%s/profile_y_0.0_0.0.png %s/profile_y_0.0_0.0_%s.png" % (output_folder, job_name, www_folder, job_name)
    print subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).stdout.read()
    cmd = "cp %s/%s/profile_x_0.0_0.0.png %s/profile_x_0.0_0.0_%s.png" % (output_folder, job_name, www_folder, job_name)
    print subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).stdout.read()
    cmd = "cp %s/%s/x-0.0.png %s/x0.0_%s.png" % (output_folder, job_name, www_folder, job_name)
    print subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).stdout.read()
    cmd = "cp %s/%s/y-0.0.png %s/y0.0_%s.png" % (output_folder, job_name, www_folder, job_name)
    print subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).stdout.read()
    cmd = "cp %s/%s/x0.0.png %s/x0.0_%s.png" % (output_folder, job_name, www_folder, job_name)
    print subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).stdout.read()
    cmd = "cp %s/%s/y0.0.png %s/y0.0_%s.png" % (output_folder, job_name, www_folder, job_name)
    print subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).stdout.read()
    cmd = "cp %s/%s/z0.0.png %s/z0.0_%s.png" % (output_folder, job_name, www_folder, job_name)
    print subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).stdout.read()
    cmd = "cp %s/%s/z-0.0.png %s/z0.0_%s.png" % (output_folder, job_name, www_folder, job_name)
    print subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).stdout.read()

    return results


def run_plot_event_data(script_path, output_folder, www_folder, job_name):
    cmd = "python %s %s/%s/events.dat %s/%s" % (script_path, output_folder, job_name, output_folder, job_name)
    results = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).stdout.read()
    ofile = open("%s/%s_events.log" % (output_folder, job_name), "w")
    ofile.write(results)
    ofile.close()

    return results


def process_combo(moniker, params_dict, varying_params_dict, template_path, output_folder, www_folder, combo):
    description = getJobDescription(combo)
    job_name = "%s-%s" % (moniker, description)
    print "running: %s " % job_name
    # myJobRunner = JobRunner(job_name)

    # myJobRunner.set_parameter('YYY_RUN_NAME', job_type)
    params_dict['YYY_RUN_NAME'] = job_name
    for key in combo.keys():
        params_dict["YYY_%s" % key] = combo[key]

    template = open(template_path, "r").read()
    cfg = setParametersInString(template, params_dict)

    config_file_path = "%s/%s.cfg" % (output_folder, job_name)
    ofile = open(config_file_path, "w")
    ofile.write(cfg)
    ofile.close()

    results = run_exe(EXECUTABLE, config_file_path, output_folder, job_name)
    results += run_plot_reconstruction(PLOT_RECON_SCRIPT, output_folder, www_folder, job_name)
    results += run_plot_event_data(PLOT_EVENTS_SCRIPT, output_folder, www_folder, job_name)


def getParameterCombos(parameters_dict, keys, run_dict, combinations, description):
    param_values = parameters_dict[keys[0]]
    print("KEYS:", keys)

    for param in param_values:
        run_dict[keys[0]] = param
        if len(keys) > 1:
            getParameterCombos(parameters_dict, keys[1:], \
                               run_dict, combinations, \
                               "%s_%s-%s" % (description, keys[0], param))
        else:
            combinations.append(collections.OrderedDict(run_dict))
    return combinations


def run_jobs(moniker, params_dict, varying_params_dict, template_path, output_folder, www_folder):

    make_folder(output_folder)
    combos = getParameterCombos(varying_params_dict, varying_params_dict.keys(), collections.OrderedDict(), [], "")

    partial_process_combo = functools.partial(process_combo, moniker, params_dict,
                                              varying_params_dict, template_path, output_folder, www_folder)
    # for combo in combos:
    #     results = partial_process_combo(combo)
    #     print results

    pool = multiprocessing.Pool(SIMULTANEOUS_JOBS)
    results = pool.map(partial_process_combo, combos)

    print results


def make_www_folder(folder, delete=False):
    make_folder(folder)
    print subprocess.Popen("cp /home/dsmackin/public_html/upenn/*.php %s" % (folder), shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).stdout.read()
    if delete:
        subprocess.Popen("rm %s/*.png" % (folder), shell=True, stdout=subprocess.PIPE,stderr=subprocess.STDOUT).stdout.read()


# -------------------------------------------------------------------------------
#    MAIN PROCESSING
# -------------------------------------------------------------------------------
def process_arguments():
    USAGE = """

    Usage: python run_batch.py


    """
    if len(sys.argv) == 2:
        job_name = sys.argv[1]

    else:
        print USAGE
        sys.exit(-10)
    return job_name


def main():
    # run_params_dict is for the parameters we wish to iterate over

    run_params_dict = collections.OrderedDict([
        ('SYSTEM_MATRIX_SCALAR', [3]),
        ('BINS', [200]),
        #('MAXIMUM_NUM_CONES', [20E7]),
        #('MAX_SCATTERING_ANGLE', [150]),
        #('MIN_SCATTERING_ANGLE', [50]),
        ('KERNEL_BANDWIDTH', [8]),
        ('IMAGE_ALGORITHM', ['KEM']),
        #('EVENT_FILE_PATH', ['20250917_run4.csv', '20250917_run6.csv', '20250917_run1.csv'])
        ('EVENT_FILE_PATH', ['run3.csv'])
    ])
    make_www_folder(WWW_FOLDER, delete=DELETE_IMAGES)
    run_jobs('dose', SCRIPT_PARAMETERS_DICT, run_params_dict, CONFIG_FILE_TEMPLATE, OUTPUT_FOLDER, WWW_FOLDER)

if __name__ == "__main__":
    print "Running as main . . ."
    # import profile
    # profile.run('main()')
    main()
