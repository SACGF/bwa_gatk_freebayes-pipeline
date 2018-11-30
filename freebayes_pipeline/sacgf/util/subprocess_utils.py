#!/usr/bin/env python

"""
For subprocess
"""

import subprocess

def run_subprocess_cmd(command, print_cmd=True, print_stdout_stderr=True, get_returncode=False):
    """
    Run the command (a string), with printing the stand out and stand error.
    Returns standard ouput and error
    returncode will also be returned if that argument is set to True
    """
    if print_cmd:
        print
        print 'Running command:\n%s' % command
        print 

    sp = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, 
            shell=True)
    out, error = sp.communicate() 
    if print_stdout_stderr:
        print
        print out
        print
        print error
        print

    if get_returncode:
        return out, error, sp.returncode
    else:
        return out, error

def run_and_handle_error(cmd, print_cmd=True):
    """
    Runs cmd, and exits with printing the standard errors if there is any.
    Do NOT use this for those applications/softwares which write stdout to stderr!
    """
    stdout, stderr = run_subprocess_cmd(cmd, print_cmd=print_cmd, print_stdout_stderr=False)
    if stderr:
        import sys; sys.exit('Standard Errors:\n%s\n' % stderr)
    return stdout, stderr

def run_and_handle_returncode(cmd, print_cmd=True):
    """
    Runs cmd, and exits if returncode != 0, with printing the standard stdout and stderr.
    Do NOT use this for those applications/softwares which always return 0
    """
    stdout, stderr, returncode = run_subprocess_cmd(cmd, print_cmd=print_cmd, 
            print_stdout_stderr=False, get_returncode=True)
    if returncode != 0:
        import sys; sys.exit('STDOUT:\n%s\nSTDERR:\n%s\n' % (stdout, stderr))
    return stdout, stderr, returncode
