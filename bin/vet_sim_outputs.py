#! /usr/bin/env python

import sys
import os
import re
import math
import logging
import glob
import argparse

import pycoevolity

import project_util


batch_number_pattern = re.compile(r'batch(?P<batch_number>\d+)')
sim_number_pattern = re.compile(r'-sim-(?P<sim_number>\d+)-')
run_number_pattern = re.compile(r'run-(?P<run_number>\d+)-')
reruns = set()


def is_var_only(path):
    return ("var-only-" in os.path.basename(path))

def log_problem_path(path):
    var_only_str = ""
    if is_var_only(path):
        var_only_str = "var-only-"
    sim_number_matches = sim_number_pattern.findall(path)
    assert(len(sim_number_matches) == 1)
    sim_number = sim_number_matches[0]
    run_number_matches = run_number_pattern.findall(path)
    assert(len(run_number_matches) == 1)
    run_number = run_number_matches[0]
    qsub_path = os.path.join(
            os.path.dirname(path),
            "{0}simcoevolity-sim-{1}-config-run-{2}-qsub.sh".format(
                    var_only_str,
                    sim_number,
                    run_number))
    if not qsub_path in reruns:
        reruns.add(qsub_path)
        sys.stdout.write("{0}\n".format(qsub_path))

def log_missing_runs(paths,
        expected_number_of_runs):
    runs = set()
    sim_numbers = set()
    for p in paths:
        run_number_matches = run_number_pattern.findall(p)
        assert(len(run_number_matches) == 1)
        run_number = int(run_number_matches[0])
        runs.add(run_number)
        sim_number_matches = sim_number_pattern.findall(p)
        assert(len(sim_number_matches) == 1)
        sim_number = sim_number_matches[0]
        sim_numbers.add(sim_number)
    assert len(sim_numbers) == 1
    sim_number = sim_numbers.pop()
    expected_runs = set(range(1, expected_number_of_runs + 1))
    missing_runs = sorted(expected_runs.difference(runs))
    for run_number in missing_runs:
        var_only_str = ""
        if is_var_only(paths[0]):
            var_only_str = "var-only-"
        qsub_path = os.path.join(
                os.path.dirname(paths[0]),
                "{0}simcoevolity-sim-{1}-config-run-{2}-qsub.sh".format(
                        var_only_str,
                        sim_number,
                        run_number))
        if not qsub_path in reruns:
            reruns.add(qsub_path)
            sys.stdout.write("{0}\n".format(qsub_path))
            sys.stderr.write("Missing run {0!r}\n".format(qsub_path))

def line_count(path):
    count = 0
    with open(path) as stream:
        for line in stream:
            count += 1
    return count

def vet_sim_rep(
        posterior_paths,
        stdout_paths,
        expected_number_of_samples = 1501):
    posterior_paths = sorted(posterior_paths)
    stdout_paths = sorted(stdout_paths)
    assert(len(posterior_paths) == len(stdout_paths))

    for pp in posterior_paths:
        if line_count(pp) != (expected_number_of_samples + 1):
            sys.stderr.write("Posterior {0!r} is incomplete\n".format(pp))
            log_problem_path(pp)

    for p in stdout_paths:
        try:
            so = pycoevolity.parsing.EcoevolityStdOut(p)
        except:
            sys.stderr.write("Stdout {0!r} is incomplete\n".format(p))
            log_problem_path(p)
            pass

def vet_simulation_output(
        expected_number_of_runs = 3,
        expected_number_of_samples = 1501):
    val_sim_dirs = glob.glob(os.path.join(project_util.VAL_DIR, '0*'))
    for val_sim_dir in sorted(val_sim_dirs):
        dpp = True
        sim_name = os.path.basename(val_sim_dir)

        batch_dirs = glob.glob(os.path.join(val_sim_dir, "batch*"))
        for batch_dir in sorted(batch_dirs):
            batch_number_matches = batch_number_pattern.findall(batch_dir)
            assert(len(batch_number_matches) == 1)
            batch_number_str = batch_number_matches[0]
            batch_number = int(batch_number_str)

            var_only_present = False
            var_only_path = glob.glob(os.path.join(val_sim_dir, "batch*",
                    "run-*-var-only-simcoevolity-sim-00*-config-state-run-1.log*"))
            if (len(var_only_path) > 0):
                var_only_present = True

            posterior_paths = glob.glob(os.path.join(batch_dir,
                    "run-1-simcoevolity-sim-*-config-state-run-1.log*"))
            for posterior_path in sorted(posterior_paths):
                sim_number_matches = sim_number_pattern.findall(posterior_path)
                assert(len(sim_number_matches) == 1)
                sim_number_str = sim_number_matches[0]
                sim_number = int(sim_number_str)
                post_paths = glob.glob(os.path.join(batch_dir,
                        "run-*-simcoevolity-sim-{0}-config-state-run-1.log*".format(
                                sim_number_str)))
                if (len(post_paths) != expected_number_of_runs):
                    log_missing_runs(post_paths,
                            expected_number_of_runs = expected_number_of_runs)
                stdout_paths = glob.glob(os.path.join(batch_dir,
                        "run-*-simcoevolity-sim-{0}-config.yml.out*".format(
                                sim_number_str)))
                if (len(stdout_paths) != expected_number_of_runs):
                    log_missing_runs(stdout_paths,
                            expected_number_of_runs = expected_number_of_runs)

                vet_sim_rep(
                    posterior_paths = post_paths,
                    stdout_paths = stdout_paths,
                    expected_number_of_samples = expected_number_of_samples)

                if var_only_present:
                    var_only_post_paths = glob.glob(os.path.join(batch_dir,
                            "run-*-var-only-simcoevolity-sim-{0}-config-state-run-1.log*".format(
                                    sim_number_str)))
                    if (len(var_only_post_paths) != expected_number_of_runs):
                        log_missing_runs(var_only_post_paths,
                            expected_number_of_runs = expected_number_of_runs)
                    vo_so_pattern = os.path.join(batch_dir,
                            "run-*-var-only-simcoevolity-sim-{0}-config.yml.out*".format(
                                    sim_number_str))
                    var_only_stdout_paths = glob.glob(vo_so_pattern)
                    if (len(var_only_stdout_paths) != expected_number_of_runs):
                        log_missing_runs(var_only_stdout_paths,
                            expected_number_of_runs = expected_number_of_runs)

                    vet_sim_rep(
                        posterior_paths = var_only_post_paths,
                        stdout_paths = var_only_stdout_paths,
                        expected_number_of_samples = expected_number_of_samples)


def main_cli(argv = sys.argv):
    parser = argparse.ArgumentParser()

    parser.add_argument('-r', '--expected-number-of-runs',
            action = 'store',
            type = int,
            default = 3,
            help = 'Number of MCMC chains that were run for each sim rep.')
    parser.add_argument('-s', '--expected-number-of-samples',
            action = 'store',
            type = int,
            default = 1501,
            help = ('Number of MCMC samples that should be found in the log '
                    'file of each analysis.'))

    if argv == sys.argv:
        args = parser.parse_args()
    else:
        args = parser.parse_args(argv)

    vet_simulation_output(
            expected_number_of_runs = args.expected_number_of_runs,
            expected_number_of_samples = args.expected_number_of_samples)


if __name__ == "__main__":
    main_cli()
