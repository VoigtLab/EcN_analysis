"""Flow cytometry import helper module

Gathers useful functions for importing fcs files into cytoflow-based analysis
of flow data. Overall functions include: finding the flow cytometry files in
provided directory; associating provided condition types to each file; create
a cytoflow-compatible dataframe (experiment) for downstream analysis.

Conditions shoudl be in 'conditions.csv' file and should be laid out as:
    Row 1: condition name
    Row 2: condition type
    Row 3 - ... : condition for given tube

Requires cytoflow to be installed in the environment.
"""

import glob
import csv
import cytoflow as flow
import pprint

def dir_fix(dir_str):
    """Fixes windows directory strings so they can be used in python

    Parameters
    ----------
    dir_str : str
        RAW directory string (can be copied and pasted from windows explorer).

    Returns
    -------
    fix_str : str
        Fixed directory string.
    """

    # First replace backslashes with forward slashes
    dir_str = dir_str.replace('\\','/')
    # add forware slash at the end
    if dir_str[-1] != '/':
        dir_str = dir_str + '/'

    return dir_str



def parse_fcs(directory):
    """Parses given directory for fcs files and returns as ordered list.

    Parameters
    ----------
    directory : str
        Full path directory to containing folder with .fcs files.
    """
    # Get all the .fcs files
    fcs_files = glob.glob(directory + '*.fcs')

    # Sort files (gotta use the last 3 characters of fortessa file outputs)
    fcs_list = sorted(fcs_files, key = lambda item: item[-7:-4])

    # Print list to make sure values are sorted correctly
    pprint.pprint(fcs_list)
    return fcs_list

def parse_conds(conds_file):
    """Parses conditions .csv file and returns list of dictionaries for each tube
        as well as the condition types as a list.

    Parameters
    ----------
    conds_file : str
        File name of .csv file with full path.

    Returns
    -------
    conds_list: list of dict
        List of dictionaries where each item has conditions corresponding to a tubeself.
    conds_type: list
        List of the condition types that each condition is.
    """

    # Create list of dictionaries for each tube.
    conds_list = []
    reader = csv.DictReader(open(conds_file))
    for line in reader:
        conds_list.append(line)

    # Pop the first item in the list off (this should be the conditions)
    conds_type = conds_list.pop(0)

    return conds_list, conds_type

def create_tube(tube_file, tube_cond):
    """Creates a cytoflow tube object from a .fcs file and the tubes conditions.

    Parameters
    ----------
    tube_file : str
        Full path filename of .fsc file
    tube_cond : dict
        Dictionary with the tube's conditions.

    Returns
    -------
    tube_obj : Cytoflow tube object
    """

    tube_obj = flow.Tube(file = tube_file, conditions = tube_cond)
    return tube_obj

def assoc_conds(fcs_files, conditions):
    """Associate fcs files and conditions file and return a cytoflow experiment.

    Parameters
    ----------
    fcs_files : list of str
        List of fcs files in string format.
    conditions : list of dict
        Object where each "row" corresponsd to a tube and is a dictionary with
        the key being the condition and the value being the condition's value
        for that tube.

    Returns
    -------
    exp: dataframe
        Pandas dataframe compatible with cytoflow experiment method.
    """
    return exp

def exp_from_dirs(fcs_dir,conds_dir=None, event_num = 0):
    """Take directories with .fcs files and experiment conditions and generate
        experiment output for cytoflow. Main function to be called in this module.

    Parameters
    ----------
    fcs_dir : str
        Directory with fcs files in it.
    conds_dir : str, optional
        Directory with conditions .csv file within it. Defaults to same directory
        as .fcs data files if not specified.
    event_num : int, optional
        Number of events to import. Defaults to all events.

    Returns
    -------
    exp: dataframe
        Pandas dataframe compatible with cytoflow experiment method.
    """

    # Grab list of all .fcs files
    fcs_list = parse_fcs(fcs_dir)

    # Get and parse conditions (use same folder as .fcs if unspecified)
    if conds_dir is None:
        conds_dir = fcs_dir
    conds_file = conds_dir + 'conditions.csv'
    conds_list, conds_type = parse_conds(conds_file)

    # Create tubes list to make experiment from
    tube_objs = []
    for i in range(0,len(fcs_list)):
        tube = create_tube(fcs_list[i], conds_list[i])
        tube_objs.append(tube)

    import_op = flow.ImportOp(conditions = conds_type, tubes = tube_objs,
                                events = event_num)
    exp = import_op.apply()
    return exp
