#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This script uses the GSC's internal API to retrieve the file paths for given
library IDs. Currently, this script handles the following library types:
  - genome: Merged genome BAM files
  - mrna: RNA-seq aligned BAM files (JaGUaR pipeline)
  - mirna: miRNA BAM files

Exomes or other unmerged BAM files probably don't work. Once there is a
need, I will look into adding support for these additional library types.

This script also needs your GSC GIN credentials. By default, it looks in
your home directory. The format is shown below:

    → cat ~/.gin_credentials.ini
    [credentials]
    username: bgrande
    password: (Re64(zuqN2.K23++A

N.B. Don't forget to make the file read-only for just the owner (chmod 500).

For more information on the GSC's API, visit the following Wiki page:
https://www.bcgsc.ca/wiki/display/bioapps/Bioapps+API+Documentation
"""

from __future__ import print_function
import sys
import os.path
import datetime
import argparse
from itertools import chain
import ConfigParser
import xmlrpclib
import glob

# Used to track is any library ID returned more than one BAM file path.
# This is used to warn the user at the end about any omitted BAM file paths.
is_multiple = False  # Global variable
IS_MULTIPLE_MSG = ("WARNING: "
    "At least one library ID returned multiple BAM files. "
    "The most recent one was returned by default. "
    "To obtain all of the BAM file paths, enable -a. "
    "To inspect the API results more closely, enable -v or -vv. "
    "To disable this warning, enable -q. "
    "For more information about these options, run with --help.")

# Track verbosity at the global level
verbosity = 0


def main():
    """Main program.

    Returns:
        None
    """
    args = parse_args()
    api = config_api(args.credentials)
    for lib_id in args.library_id:
        bam_paths = get_bam_paths(api, lib_id, args.library_type)
        log_bam_paths(bam_paths)
        if not args.all and len(bam_paths) > 1:
            global is_multiple
            is_multiple = True
            bam_paths = bam_paths[0:1]
        for path in bam_paths:
            oline = "{}\t{}\n".format(lib_id, path)
            args.output.write(oline)
    args.output.close()
    if is_multiple and not args.quiet:
        log(IS_MULTIPLE_MSG)


def parse_args():
    """Parse command-line arguments.

    Returns:
        A namespace containing the arguments
    """
    # Default values
    default_config = os.path.expanduser("~/.gin_credentials.ini")
    # Parse command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("library_id", nargs="+", help="List of library IDs (space separated).")
    parser.add_argument("--library_type", "-t", choices=["genome", "mrna", "mirna"],
                        default="genome", help="Library type for given library IDs. [genome]")
    parser.add_argument("--output", "-o", type=argparse.FileType("w"), default="-",
                        help="Output TSV file with library ID and BAM file path. [stdout]")
    parser.add_argument("--credentials", "-c", default=default_config,
                        help="INI configuration with GIN credentials. See script header for "
                        "example. [{}]".format(default_config))
    parser.add_argument("--all", "-a", action="store_true", help="Output all BAM file paths.")
    parser.add_argument("--verbose", "-v", action="count",
                        help="Increase the script's verbosity. Can use multiple times (-vv)")
    parser.add_argument("--quiet", "-q", action="store_true",
                        help="Mute multiple BAM file warning.")
    args = parser.parse_args()
    # Handle special arguments
    if not os.path.exists(args.credentials):
        raise ValueError("Configuration file for GIN credentials doesn't exist: "
                         "{}".format(args.credentials))
    global verbosity
    verbosity = args.verbose
    # Return arguments
    return args


def config_api(credentials_path):
    """Configure and create the API connection.

    Returns:
        A xmlrpclib.ServerProxy instance
    """
    # Obtain GIN credentials from file
    config = ConfigParser.ConfigParser()
    config.read(credentials_path)
    gin_username = config.get("credentials", "username")
    gin_password = config.get("credentials", "password")
    # Configure and create the API connection
    url_template = 'http://{}:{}@www.bcgsc.ca/data/sbs/viewer/api'
    url = url_template.format(gin_username, gin_password)
    api = xmlrpclib.ServerProxy(url)
    return api


def parse_date(s):
    """Parse timestamps as returned by the GSC API.

    Returns:
        A datetime object
    """
    return datetime.datetime.strptime(s, "%Y-%m-%d %H:%M")


def log(msg):
    """Prints a message to stderr.

    Returns:
        None
    """
    print(msg, file=sys.stderr)


def log_library_id(lib_id):
    """Log library ID.

    Returns:
        None
    """
    global verbosity
    if verbosity >= 1:
        log("Library ID: {}".format(lib_id))


def log_bam_paths(bam_paths):
    """Log BAM file paths.

    Returns:
        None
    """
    global verbosity
    if verbosity >= 1:
        log("Number of BAM paths: {}".format(len(bam_paths)))
    if verbosity >= 2:
        log("All BAM Paths: \n  - {}".format("\n  - ".join(bam_paths)))


def log_api_results(lib_info):
    """Log API results.

    Returns:
        None
    """
    global verbosity
    if verbosity >= 3:
        log("API result: {}".format(lib_info))


def find_bam(dirname, file_glob):
    """Return a BAM file path according to the given directory
    and file name glob.
    """
    if dirname is None:
        return None
    full_glob = os.path.join(dirname, file_glob)
    bam_file = glob.glob(full_glob)
    if len(bam_file) < 1:
        return None
    elif len(bam_file) > 1:
        raise ValueError("More than one file matches the glob {}".format(full_glob))
    return bam_file[0]


def get_bam_paths(api, lib_id, lib_type):
    """Retrieve all BAM file paths using API.

    The appropriate function will be used according to the
    lirbary type.

    Returns:
        List of BAM file paths (strings), ["N/A"] if none are returned
    """
    log_library_id(lib_id)
    if lib_type == "mirna":
        bam_paths = get_mirnaseq_bam(api, lib_id)
    elif lib_type == "mrna":
        bam_paths = get_rnaseq_bam(api, lib_id)
    else:
        bam_paths = get_genome_bam(api, lib_id)
    if bam_paths is None or len(bam_paths) == 0:
        bam_paths = ["N/A"]
    return bam_paths


def get_genome_bam(api, lib_id):
    """Retrieve the path of a merged genome BAM file.

    The list is sorted in order of creation (newest first).

    Returns:
        List of BAM file paths (strings), can be empty
    """
    lib_info = api.getMetaInfo({"library": lib_id})
    if not lib_info:
        return None
    libs = chain.from_iterable([x.values() for x in lib_info.values()])
    libs_success = [lib for lib in libs if lib["success"]]
    libs_success = sorted(libs_success[::-1], key=lambda x: parse_date(x["process_complete"]))
    if len(libs_success) == 0:
        return None
    log_api_results(lib_info)
    bam_paths = []
    for lib in libs_success:
        lib_bam_dir = lib["data_path"]
        lib_bam_file = find_bam(lib_bam_dir, "*.bam")
        if lib_bam_file:
            bam_paths.append(lib_bam_file)
    return bam_paths


def get_rnaseq_bam(api, lib_id):
    """Retrieve the path of a RNA-seq BAM file aligned using
    the JaGUaR pipeline.

    The list is sorted in order of creation (newest first).

    Returns:
        List of BAM file paths (strings), can be empty
    """
    lib_info = api.getAlignedLibcoreInfo({"library": lib_id})
    if not lib_info:
        return None
    libs = chain.from_iterable(lib_info.values())
    libs_success = [lib for lib in libs if lib["successful"]]
    libs_success = sorted(libs_success[::-1], key=lambda x: parse_date(x["created"]))
    if len(libs_success) == 0:
        return None
    log_api_results(lib_info)
    bam_paths = []
    for lib in libs_success:
        lib_bam_dir = lib["bioapps_data_path"]
        lib_bam_file = find_bam(lib_bam_dir, "*withJunctionsOnGenome_dupsFlagged.bam")
        if lib_bam_file:
            bam_paths.append(lib_bam_file)
    return bam_paths


def get_mirnaseq_bam(api, lib_id):
    """Retrieve the path of a miRNA-seq BAM file.

    The list is sorted in order of creation (newest first).

    Returns:
        List of BAM file paths (strings), can be empty
    """
    lib_info = api.getAlignedLibcoreInfo({"library": lib_id})
    if not lib_info:
        return None
    libs = chain.from_iterable(lib_info.values())
    libs_success = [lib for lib in libs if lib["successful"]]
    libs_success = sorted(libs_success[::-1], key=lambda x: parse_date(x["created"]))
    if len(libs_success) == 0:
        return None
    log_api_results(lib_info)
    bam_paths = []
    for lib in libs_success:
        lib_bam_dir = lib["bioapps_data_path"]
        lib_bam_file = find_bam(lib_bam_dir, "*.bam")
        if lib_bam_file:
            bam_paths.append(lib_bam_file)
    return bam_paths


if __name__ == '__main__':
    main()
