#! /usr/bin/env python

import sys
import os
import tarfile
import re
import glob
import logging

logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.DEBUG)
_LOG = logging.getLogger(os.path.basename(__file__))

import pycoevolity
import project_util

def parse_msbayes_results(sim_posterior_path,
        out_path):
    sim_dir = os.path.dirname(sim_posterior_path)
    sim_file = os.path.basename(sim_posterior_path)
    sim_file_pattern = re.compile(r'.*d1-m1-s(?P<sim_number>\d+)-(?P<nprior>\d+)-.*')
    sim_file_match = sim_file_pattern.match(sim_file) 
    prior_sample_size = int(sim_file_match.group("nprior"))
    data_model_dir = os.path.dirname(sim_dir)
    output_dir = os.path.dirname(data_model_dir)
    results_dir = os.path.dirname(output_dir)
    true_path = os.path.join(results_dir,
            "observed-summary-stats",
            "observed-1.txt")
    post_path_wild_card = os.path.join(
            sim_dir,
            "d1-m1-s*-{prior_sample_size}-posterior-sample.txt.gz".format(
                    prior_sample_size = prior_sample_size))
    posterior_paths = glob.glob(post_path_wild_card)
    number_of_sims = len(posterior_paths)
    number_of_pairs = 0
    for header in pycoevolity.parsing.parse_header_from_path(sim_posterior_path):
        if header.startswith("PRI.t."):
            number_of_pairs += 1

    _LOG.info("Parsing {0} posterior files with {1} pairs in {2!r}...".format(
            number_of_sims,
            number_of_pairs,
            sim_dir))

    column_header_prefixes = {
            'root_height': ("PRI.t.",),
            'pop_size_root': ("PRI.aTheta.",),
            'pop_size': ("PRI.d1Theta.", "PRI.d2Theta."),
            }
    parameter_prefixes = [
            "true",
            "mean",
            "eti_95_lower",
            "eti_95_upper",
            ]
    results = {
            "true_num_events": [],
            "map_num_events": [],
            "true_num_events_cred_level": [],
            }
    for pair_idx in range(number_of_pairs):
        results["num_events_{0}_p".format(pair_idx + 1)] = []
        for parameter_key, header_prefixes in column_header_prefixes.items():
            result_keys = ["{param}_c{comparison}sp1".format(
                            param = parameter_key,
                            comparison = pair_idx + 1)]
            if parameter_key== "pop_size":
                result_keys.append("{param}_c{comparison}sp2".format(
                                param = parameter_key,
                                comparison = pair_idx + 1))
            for suffix in result_keys:
                for prefix in parameter_prefixes:
                    k = "{0}_{1}".format(prefix, suffix)
                    results[k] = []
    true_values = pycoevolity.parsing.get_dict_from_spreadsheets(
            [true_path],
            sep = '\t')
    nsims = len(true_values.values()[0])
    assert (nsims == number_of_sims)
    for sim_idx in range(nsims):
        posterior_path = os.path.join(sim_dir,
                "d1-m1-s{sim_num}-{prior_sample_size}-posterior-sample.txt.gz".format(
                        sim_num = sim_idx + 1,
                        prior_sample_size = prior_sample_size))
        _LOG.info("Parsing {0}".format(posterior_path))
        posterior = pycoevolity.parsing.get_dict_from_spreadsheets(
                [posterior_path],
                sep = '\t')
        for pair_idx in range(number_of_pairs):
            for parameter_key, header_prefixes in column_header_prefixes.items():
                for header_prefix in header_prefixes:
                    header = "{0}{1}".format(header_prefix, pair_idx + 1)
                    pop_number = 1
                    if "d2Theta" in header:
                        pop_number = 2
                    result_key = "{param}_c{comparison}sp{pop}".format(
                            param = parameter_key,
                            comparison = pair_idx + 1,
                            pop = pop_number)
                    if parameter_key.startswith("pop_size"):
                        true_val = float(true_values[header][sim_idx]) / 4.0
                        post_sum = pycoevolity.stats.get_summary((float(x) / 4.0) for x in posterior[header])
                    else:
                        true_val = float(true_values[header][sim_idx])
                        post_sum = pycoevolity.stats.get_summary(float(x) for x in posterior[header])
                    results["true_" + result_key].append(true_val)
                    results["mean_" + result_key].append(post_sum['mean'])
                    results["eti_95_lower_" + result_key].append(post_sum['qi_95'][0])
                    results["eti_95_upper_" + result_key].append(post_sum['qi_95'][1])
        true_nevents = int(true_values["PRI.Psi"][sim_idx])
        nevent_freqs = pycoevolity.stats.get_freqs(int(x) for x in posterior["PRI.Psi"])
        sorted_nevent_freqs = sorted(nevent_freqs.items(), reverse = True,
                key = lambda x: x[1])
        map_nevents = sorted_nevent_freqs[0][0]
        true_nevents_prob = nevent_freqs[true_nevents]
        map_nevents_prob = nevent_freqs[map_nevents]
        results["true_num_events"].append(true_nevents)
        results["map_num_events"].append(map_nevents)
        for pair_idx in range(number_of_pairs):
            results["num_events_{0}_p".format(pair_idx + 1)].append(
                    nevent_freqs.get(pair_idx + 1, 0.0))
        cred_level = 0.0
        for n, p in sorted_nevent_freqs:
            if n == true_nevents:
                break
            cred_level += p
        results["true_num_events_cred_level"].append(cred_level)

    for parameter, values in results.items():
        assert len(values) == number_of_sims

    header = sorted(results.keys())

    with open(out_path, 'w') as stream:
        stream.write("{0}\n".format("\t".join(header)))
        for sim_index in range(number_of_sims):
            stream.write("{0}\n".format(
                    "\t".join(str(results[k][sim_index]) for k in header)))

def posterior_files(members):
    for tar_info in members:
        if (("500000-posterior-sample.txt" in tar_info.name) or
                ("observed-1.txt" in tar_info.name)):
            yield tar_info


def main_cli(argv = sys.argv):
    number_of_sites_pattern = ["00500", "01000", "02000", "10000"]
    for nsites in number_of_sites_pattern:
        tar_dir = os.path.join(project_util.BAKE_OFF_DIR,
                "results-{0}".format(nsites))
        if not os.path.isdir(tar_dir):
            tar_path = tar_dir + ".tar.gz"
            tar = tarfile.open(tar_path, "r:gz")
            tar.extractall(
                    path = os.path.dirname(tar_path),
                    members = posterior_files(tar))
            tar.close()
        sim_posterior_path = os.path.join(project_util.BAKE_OFF_DIR,
                "results-{0}".format(nsites),
                "dpp-msbayes",
                "pymsbayes-results",
                "pymsbayes-output",
                "d1",
                "m1",
                "d1-m1-s1-500000-posterior-sample.txt.gz")
        out_path = os.path.join(project_util.BAKE_OFF_DIR,
                "results-{0}.csv".format(nsites))
        parse_msbayes_results(
                sim_posterior_path = sim_posterior_path,
                out_path = out_path)


if __name__ == "__main__":
    main_cli()
