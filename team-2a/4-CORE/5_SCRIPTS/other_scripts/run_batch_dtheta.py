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

# -------------------------------------------------------------------------------
#  RUN PARAMETERS
# -------------------------------------------------------------------------------
WORKING_DIR = '/home/dsmackin/projects/CORE'
OUTPUT_FOLDER = "%s/dtheta" % WORKING_DIR
EXECUTABLE = '%s/core' % WORKING_DIR
INPUT_FILE_FOLDER = '%s/input' % WORKING_DIR
# EVENT_FILE_PATH = '/y_drive/CCData/UPenn20150725/fwbw.csv'
# EVENT_FILE_PATH = '/y_drive/CCData/UPenn20150725/small.csv'
EVENT_FILE_PATH = '/y_drive/CCData/gold/Co60_cc2csv_2x.csv'
EVENT_FILE_PATH = '/home/dsmackin/public_html/Co60_gold/theta_YYY_TYPE.csv'
EVENT_FILE_PATH = '/home/dsmackin/public_html/dtheta/theta_YYY_TYPE.csv'

PLOT_RECON_SCRIPT = "/home/dsmackin/projects/CORE/scripts/recon_plotter.py"
PLOT_EVENTS_SCRIPT = "/home/dsmackin/projects/CORE/scripts/plot_events.py"

WWW_FOLDER = "/home/dsmackin/public_html/dtheta"
DELETE_IMAGES = True

TEMPLATE_FILE = 'template.cfg'
CONFIG_FILE = 'generated.cfg'

SIMULTANEOUS_JOBS = 4

SCRIPT_PARAMETERS_DICT = {
    'YYY_MONIKER': 'Co60',
    'YYY_EVENT_FILE_PATH': EVENT_FILE_PATH,
    'YYY_OUTPUT_FOLDER_PATH': OUTPUT_FOLDER,
    'YYY_OUTPUT_BINS_X' : '100',
    'YYY_OUTPUT_BINS_Y' : '100',
    'YYY_OUTPUT_BINS_Z' : '100',
    'YYY_IMAGE_ALGORITHM': 'OCTANE',
    'YYY_DCA_CENTER_X': 0,
    'YYY_DCA_CENTER_Y': 0,
    'YYY_DCA_CENTER_Z': 0.0,

    #OCTANE
    'YYY_INTERCEPT_DCA': 2,
    'YYY_PHANTOM_LENGTH': 256,
    'YYY_PHANTOM_CENTER_X': 0,
    'YYY_PHANTOM_CENTER_Y': 0,
    'YYY_PHANTOM_CENTER_Z': 0,
    'YYY_PHANTOM_BINS': 32,
    'YYY_SOURCE_AXIS_DISTANCE': 400,
    'YYY_CONE_LENGTH_CORRECTION': 1.03,
    'YYY_NUMBER_OF_THREADS': 5,

    'YYY_TEMPERATURE': 0.9,
    'YYY_INVERSE_SQUARE_PARAM': 0.9,

    'YYY_NUMBER_TRIES_FOR_RANDOM': '100',
    'YYY_ITERATIONS': 500,
    'YYY_EVENT_MULTIPLIER': 1,
    'YYY_DENSITY_ESTIMATOR_TYPE': 2,
    'YYY_NUMBER_OF_SHIFTS': 200,
    'YYY_MAXIMUM_NUMBER_CONES': 10000,
    'YYY_NUM_CONES_OFFSET': 0,
    'YYY_X_MIN': -50,
    'YYY_X_MAX': 50,
    'YYY_X_BINS': 20,
    'YYY_Y_MIN': -50,
    'YYY_Y_MAX': 50,
    'YYY_Y_BINS': 20,
    'YYY_Z_MIN': -50,
    'YYY_Z_MAX': 50,
    'YYY_Z_BINS': 20,
    'YYY_USE_PARABOLAS': 1,
    'YYY_MIN_GAMMA_ENERGY': 0.6,
    'YYY_MAX_GAMMA_ENERGY': 1.5,
    'YYY_MAX_SCATTERING_ANGLE': 180,
    'YYY_DATA_FILE_FORMAT': 3,
    'YYY_SCATTER_DISTANCE': 10,
    'YYY_DCA_CUT': 200,
    'YYY_X_DCA_CUT': 0,
    'YYY_Y_DCA_CUT': 0,
    'YYY_Z_DCA_CUT': 0,
    'YYY_MAX_ENERGY_LOST': 100.0,
    #DCA LINE CUT
    'YYY_KNOWN_GAMMA_ENERGIES': '1.17, 1.33',

    #DCA for line is segment between BEAM_LINE_POINT1 and BEAM_LINE_POINT2
    'YYY_BEAM_LINE_POINT1': '0, 0, -200',
    'YYY_BEAM_LINE_POINT2': '0, 0, 200',
    'YYY_DCA_LINE_CUT': '15.0',
    'YYY_MIN_ENERGY_SCATTER': '0.30',
    'YYY_MIN_ENERGY_EVENT': '0.6',
    'YYY_X_MIN': -200,
    'YYY_X_MAX': 200,
    'YYY_X_BINS': 10,
    'YYY_Y_LENGTH': 200,
    #'YYY_Y_MAX': 50,
    'YYY_Y_BINS': 10,
    'YYY_Z_MIN': -200,
    'YYY_Z_MAX': 200,
    'YYY_Z_BINS': 10,
}

CONFIG_FILE_TEMPLATE = "%s/%s" % (WORKING_DIR, TEMPLATE_FILE)

# -------------------------------------------------------------------------------
#    FUNCTIONS
# -------------------------------------------------------------------------------
def getParameterCombos(parameters_dict, keys, run_dict, combinations, description):
    param_values = parameters_dict[keys[0]]

    for param in param_values:
        run_dict[keys[0]] = param
        if len(keys) > 1:
            getParameterCombos(parameters_dict, keys[1:], \
                               run_dict, combinations, \
                               "%s_%s-%s" % (description, keys[0], param))
        else:
            combinations.append(dict(run_dict))
    return combinations


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

    cmd = "cp %s/%s/profile_z_0.0_0.0.png %s/profile_z_0.0_0.0_%s.png" % (output_folder, job_name, www_folder, job_name)
    print subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).stdout.read()
    cmd = "cp %s/%s/profile_y_0.0_0.0.png %s/profile_y_0.0_0.0_%s.png" % (output_folder, job_name, www_folder, job_name)
    print subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).stdout.read()
    cmd = "cp %s/%s/profile_x_0.0_0.0.png %s/profile_x_0.0_0.0_%s.png" % (output_folder, job_name, www_folder, job_name)
    print subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).stdout.read()
    cmd = "cp %s/%s/x0.0.png %s/x0.0_%s.png" % (output_folder, job_name, www_folder, job_name)
    print subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).stdout.read()
    cmd = "cp %s/%s/y0.0.png %s/y-0.0_%s.png" % (output_folder, job_name, www_folder, job_name)
    print subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).stdout.read()
    cmd = "cp %s/%s/z0.0.png %s/z-0.0_%s.png" % (output_folder, job_name, www_folder, job_name)
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

    return results


def run_jobs(moniker, params_dict, varying_params_dict, template_path, output_folder, www_folder):

    make_folder(output_folder)
    combos = getParameterCombos(varying_params_dict, varying_params_dict.keys(), {}, [], "")

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
    run_params_dict = {
        #'IMAGE_ALGORITHM': ['OCTANE'],
        # 'TEMPERATURE': ['1.0', '0.75', '0.5'],
        # 'INVERSE_SQUARE_PARAM': ['1.0', '0.8', '0.75', '0.7'],
        # 'r': ['01', '02', '03', '04', '05'],
        #'TEMPERATURE': ['1.0', '0.8', '0.6'],
        #'TEMPERATURE': ['0.8'],
        #'INVERSE_SQUARE_PARAM': ['0.9','0.8','0.7'],
        #'INVERSE_SQUARE_PARAM': ['0.85'],
        'MAXIMUM_NUMBER_CONES': ['5000', '10000'],
        # 'TYPE': [ 'm', 'c_dE', 'c_dTheta', 'c_dca', 'm_dca'],
        # 'TYPE': ['c_dca', 'm', 'm_dE', 'm_dTheta'],
        # 'TYPE': ['m_dE'],
        # 'MIN_GAMMA_ENERGY': [1.0, 0.0],
        # 'MAX_GAMMA_ENERGY': [1.5, 10.0],
        # 'INTERCEPT_DCA': [4, 6, 8, 16],
        # 'PHANTOM_BINS': [32,64,128]
        'TYPE': [ 'm', 'm_dE', 'm_dTheta'],
        'CONE_LENGTH_CORRECTION': ['0.90', '1.0'],
        #'TYPE': [ 'm_cutSinDelta', 'm_cutDE', 'm'],
        # 'SCATTER_DISTANCE': ['10', '05'],
        #'X_BINS':['10','14'],#
        #'Y_BINS':['10','14'],
        #'Z_BINS':['10','14'],
        #'EVENT_MULTIPLIER': ['01', '10', '20']
    }
    make_www_folder(WWW_FOLDER, delete=DELETE_IMAGES)
    run_jobs('dtheta', SCRIPT_PARAMETERS_DICT, run_params_dict, CONFIG_FILE_TEMPLATE, OUTPUT_FOLDER, WWW_FOLDER)


if __name__ == "__main__":
    print "Running as main . . ."
    # import profile
    # profile.run('main()')
    main()
