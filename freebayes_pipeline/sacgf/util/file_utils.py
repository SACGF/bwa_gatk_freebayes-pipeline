import glob
import os
import random
import re
import string

from sacgf.util.subprocess_utils import run_and_handle_error


def file_or_file_name(f, mode='r'):
    if isinstance(f, basestring): # Works on unicode
        if 'w' in mode: # Create path if writing
            mk_path_for_file(f)
            
        return open(f, mode)
    elif isinstance(f, file):
        return f # Already a File object
    else:
        raise ValueError("'%s' (%s) not a file or string" % (f, type(f)))

def mk_path(path):
    if path and not os.path.exists(path):
        os.makedirs(path)

def mk_path_for_file(f):
    mk_path(os.path.dirname(f))

def name_from_file_name(file_name):
    '''Gets file name with removing extension and directory'''
    return os.path.splitext(os.path.basename(file_name))[0]


def nice_name_from_file_name(file_name):
    ''' strips sequencing crap from name '''
    
    name = name_from_file_name(file_name)
    name = re.sub("_R[12]", "", name    )
    name = re.sub("_L0\d\d_.*", "", name) # Lane
    name = re.sub("_[GATC]+_.*", "", name) # Barcodes
    name_no_underscores = name.replace("_", "")
    if not name_no_underscores:
        msg = "Stripped too much from '%s' - left with nothing..." % file_name
        raise ValueError(msg)

    return name

def file_to_hash(f, sep=None):
    data = {}
    for line in file_or_file_name(f):
        line = line.rstrip()
        key = line
        value = None
        if sep:
            items = line.split(sep)
            if len(items) >= 2:
                (key, value) = items[0:2]
        data[key] = value
    return data

def file_column_to_hash(f, key_column_id, value_column_id, sep=','):
    data = {}
    for line in file_or_file_name(f):
        line = line.rstrip()
        cols = line.split(sep)
        key = cols[key_column_id]
        value = cols[value_column_id]
        data[key] = value
    return data

def find_files_in_dir(dir_name, start='', end='', full_path=True):
    """
    returns all files in dir_name with starting with 'start' and ending with 'end' (by default,
    returning all files in that dir).
    Set full_path to False to return only file names rather than absolute_file_path
    """
    dir_name = os.path.expanduser(dir_name)
    files = [f for f in os.listdir(dir_name) if f.startswith(start) and f.endswith(end)]
    if not full_path:
        return files
    else:
        return [os.path.join(dir_name, f) for f in files]
    
def single_file_from_wildcard(pathname):
    files = glob.glob(pathname)
    num_files = len(files)
    if num_files == 1:
        return files[0]

    problem = "No matches"
    if num_files > 1:
        problem = "More than 1 (%d) matches" % num_files
    raise ValueError("%s for filename pattern '%s'" % (problem, pathname))

def line_count_wc(filename):
    '''
    Use wc, the unix tool, to count and return the number of lines of a file.
    if the file is .gz, count the line number of the uncompressed file
    '''
    if filename.endswith('.gz'):
        cmd = 'zcat %s | wc -l -' % filename
    else:
        cmd = 'wc -l %s' % filename

    out, _ = run_and_handle_error(cmd, print_cmd=False)
    return int(out.split(' ')[0])

def line_count(file_name):
    lines = 0
    with open(file_name, "r") as f:
        for _ in f:
            lines += 1
    return lines

def random_file_name():
    """
    returns a random string of 8 characters which can be used as a 'unique' file name
    """
    return ''.join(random.choice(string.ascii_uppercase) for _ in range(8))

def random_file_at_same_dir(file_path, prefix='', extension=''):
    """
    returns a file path with random file name and in the same dir as file_path. 
    prefix and extension will be added to the random file name, if available.
    """
    fname = random_file_name()
    if prefix:
        fname = prefix + '.' + fname
    if extension:
        fname = fname + '.' + extension
    return os.path.join(os.path.dirname(file_path), fname)

def remove_a_file(path):
    """remove a file if that file is existing and removable"""
    try:
        os.remove(path)
        return True
    except:
        return False

def remove_gz_if_exists(filename):
    if filename.endswith(".gz"):
        filename = os.path.splitext(filename)[0] # remove .gz
    return filename

def absolute_file_path(file_path):
    """return the absolute file path, with converting '~' to user's home directory.
    if file_path is not a string, just return file_path per se """
    if type(file_path) is not str:
        return file_path
    else:
        return os.path.abspath(os.path.expanduser(file_path))


def array_to_file(file_name, array):
    mk_path_for_file(file_name)
    with open(file_name, "w") as f:
        for line in array:
            f.write(line + "\n")

def file_to_array(file_name, comment=None):
    array = []
    with open(file_name) as f:
        for line in f:
            if comment and line.startswith(comment):
                continue
            array.append(line.rstrip())
    return array

def bioinformatics_dirname():
    ''' Returns the full path name of the bioinformatics project '''
    return up_dir(os.path.dirname(__file__), 2) # Relative to current file (change if you ever move!)

def up_dir(path, n):
    ''' Basically like os.path.join(path, "../" * n) but without the dots '''
    assert n >= 0
    path_components = path.split(os.path.sep)
    return os.path.sep.join(path_components[:-n])

def file_name_insert(insert, file_name):
    ''' Return "insert" just before the file extension '''
    (file_path, extension) = os.path.splitext(file_name)
    out_filename = "%s.%s%s" % (file_path, insert, extension)
    return out_filename
