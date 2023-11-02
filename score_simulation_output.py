#!/usr/bin/env python3
from csv import reader
from os.path import isdir, isfile
from sys import argv
import re

# load calibration data from CSV; calibration_data[description] = (value, weight)
def load_calibration_data(calibration_csv):
    calibration_data = dict()
    for row in reader(open(calibration_csv)):
        d, v, w = [x.strip() for x in row]
        if d.lower() == 'description':
            continue # skip header row
        if d in calibration_data:
            raise ValueError("Duplicate entry in calibration CSV: %s" % d)
        calibration_data[d] = (float(v), float(w))
    return calibration_data

# load calibration data output by the ABM R code; data[metric][year][subgroup] = count
def load_abm_calibration_output(abm_calibration_tsv):
    data = dict()
    for row in reader(open(abm_calibration_tsv), delimiter='\t'):
        metric, month, subgroup, stat = [x.strip() for x in row]
        if metric.lower() == 'metric':
            continue # skip header row
        if metric not in data:
            data[metric] = dict()
        year = int(month) // 12
        if year not in data[metric]:
            data[metric][year] = dict()
        if subgroup not in data[metric][year]:
            data[metric][year][subgroup] = 0
        data[metric][year][subgroup] += int(stat)
    return data

# scoring function we want to optimize
def score(sim_out_folder, calibration_csv, verbose=True):
    calibration_data = load_calibration_data(calibration_csv)
    abm_calibration_data = load_abm_calibration_output('%s/abm_hiv_calibration_data.tsv' % sim_out_folder)
    score = 0
    sim_start_time = int(sorted(k for k in calibration_data.keys() if re.match(r'[0-9]{4}_',k))[0].split('_')[0])
    print("Calibration Key\tReal Value\tSimulation Value")
    for cal_key, cal_tup in calibration_data.items():
        sim_val = None; cal_val, cal_w = cal_tup
        if '_' in cal_key:
            try:
                year, risk = cal_key.split('_')
                year = int(year)
                risk = risk.replace(' & ', 'and')
                if risk.lower() == 'other':
                    risk = 'other' # 'other' is lowercase in the ABM R code output
                sim_val = abm_calibration_data['newinfects_agg'][year-sim_start_time][risk]
            except:
                pass # sim_val will remain None, so ValueError below will be thrown
        if sim_val is None:
            raise ValueError("Unknown calibration key: %s" % cal_key)
        else:
            if verbose:
                print("%s\t%s\t%s" % (cal_key, cal_val, sim_val))
            score += cal_w * ((sim_val-cal_val)**2)
    score = score**0.5
    if verbose:
        print("Overall Calibration Score\tN/A\t%s" % score)
    return score

# main execution
if __name__ == "__main__":
    if len(argv) != 3:
        print("%s <sim_out_folder> <calibration_csv>" % argv[0]); exit(1)
    if not isdir(argv[1]):
        print("Folder not found: %s" % argv[1]); exit(1)
    if not isfile(argv[2]):
        print("File not found: %s" % argv[2]); exit(1)
    s = score(argv[1], argv[2])
