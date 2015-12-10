#!/usr/bin/env python2

"""
get_kronos_times.py
===================
This script obtains the runtimes for jobs run using Kronos.

Example Configuration
---------------------
This is an example configuration describing the Strelka
pipeline. Note that parallel tasks can be written on the
same line separated by commas; the max runtime among parallel
tasks will be reported for that line. Make sure that the tasks
you include are run on the cluster (and not locally).

```
TASK_CLIPOVERLAP_TUMOUR,TASK_CLIPOVERLAP_NORMAL
TASK_SAMTOOLS_INDEX_TUMOUR,TASK_SAMTOOLS_INDEX_NORMAL
TASK_STRELKA
TASK_MUTATIONSEQ
TASK_VCF2MAF_SNVS,TASK_VCF2MAF_INDELS
TASK_VCF2MAF_SNVS_MUTATIONSEQ
```
"""

import argparse
import sys
import os
import glob
import re
import subprocess
from datetime import datetime

DEVNULL = open("/dev/null", "w")


def main():
    args = parse_args()
    steps = parse_config(args.config)
    samples = list_samples(args.pipeline_dir)
    with args.output as output:
        output.write("\t".join(["step"] + samples) + "\n")
        for step in steps:
            cols = [",".join(step)]
            for sample in samples:
                job_ids = get_job_ids(step, os.path.join(args.pipeline_dir, sample))
                times = [get_runtime(job_id) for job_id in job_ids] or [None]
                if "IP" in times:
                    maxtime = "IP"
                elif "ERR" in times:
                    maxtime = "ERR"
                else:
                    maxtime = max(times)
                cols.append(str(maxtime))
            output.write("\t".join(cols) + "\n")


def parse_args():
    """Parse command-line arguments, validate them and prepare them."""
    parser = argparse.ArgumentParser()
    parser.add_argument("config", type=argparse.FileType(), help="Configuration (see script header for example)")
    parser.add_argument("pipeline_dir", help="Pipeline run directory")
    parser.add_argument("--output", "-o", type=argparse.FileType("w"), default=sys.stdout, help="Output file")
    args = parser.parse_args()
    # Confirm the existence of the input directory
    assert os.path.exists(args.pipeline_dir), "Pipeline directory doesn't exist"
    return args


def parse_config(config_file):
    """Parse configuration file and return a list of steps."""
    return [l.strip().split(",") for l in config_file]


def list_samples(pipeline_dir):
    """Return a list of samples from a pipeline run"""
    # Get list of samples
    samples = []
    for item in os.listdir(pipeline_dir):
        if os.path.isdir(os.path.join(pipeline_dir, item)):
            samples.append(os.path.basename(item))
    return samples


def get_job_ids(step, sample_dir):
    """Get latest job IDs for given tasks and sample directory"""
    id_pattern = re.compile(r".*\.o(\d+)")
    job_ids = []
    for task in step:
        logs_dir = os.path.join(sample_dir, "logs")
        task_files = glob.glob("{}/{}*.o*".format(logs_dir, task))
        if task_files:
            recent = sorted(task_files, key=lambda x: os.stat(x).st_mtime)[-1]
            job_id = id_pattern.match(recent).group(1)
        else:
            continue
        job_ids.append(job_id)
    return job_ids


def query_qacct(job_id):
    """Query qacct"""
    attrs = {}
    pattern = re.compile(r"(\w+)\s+(.*)")
    cmd = ["qacct", "-j", job_id]
    try:
        cmd_output = subprocess.check_output(cmd, stderr=DEVNULL)
    except subprocess.CalledProcessError:
        return attrs
    for l in cmd_output.split("\n"):
        match = pattern.match(l)
        if match:
            attrs[match.group(1)] = match.group(2).strip()
    return attrs


def get_runtime(job_id):
    """Get runtime for job_id"""
    attrs = query_qacct(job_id)
    if not attrs:
        return "IP"
    elif attrs["exit_status"] != "0":
        return "ERR"
    date_fmt = "%a %b %d %H:%M:%S %Y"
    stime = datetime.strptime(attrs["start_time"], date_fmt)
    etime = datetime.strptime(attrs["end_time"], date_fmt)
    delta = etime - stime
    return int(delta.total_seconds())


if __name__ == '__main__':
    main()
