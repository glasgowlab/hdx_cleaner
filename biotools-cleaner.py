import pandas as pd
pd.options.mode.chained_assignment = None
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import scipy.optimize as opt
import argparse
import os

# DEFAULT CONSTS

RTcutoff = 0.5
MZcutoff = 0.1
score_threshold = 150
ppm_threshold = 7
maxfev = 800
plots = 'csvplots'
outpath = 'CLEAN.csv'

# HELPER FNS

def inthelper(n):
    if n == '-':
        return 0
    return int(n.strip())

def invfn(x, a, b, c):
    return c + a / (x - b)

def redflag(topscores, r):
    # return true if they're closer in RT and mass than the cutoffs and if they're nonidentical and nonoverlapping
    redundant = (abs(topscores['Rt(min)'].apply(float) - float(r['Rt(min)'])) < RTcutoff) & (abs(topscores['Meas. M/z'].apply(float) - float(r['Meas. M/z'])) < MZcutoff) & (topscores['Sequence'] != r['Sequence'])
    nonoverlapping = (topscores['Start'].apply(int) > int(r['End'])) | (int(r['Start']) > topscores['End'].apply(int))
    return redundant & nonoverlapping

# CMD LINE ARGS CODE

parser = argparse.ArgumentParser(description='Biotools CSV cleaner: pool csvs and disambiguate duplicates')
parser.add_argument('--t', '--table', dest='table', help="paths to input csvs", nargs='+', required=True)
parser.add_argument('--o', '--out', dest='out', help="path to output csv")
parser.add_argument('--p', '--plots', dest='plots', help="directory for output plots")
parser.add_argument('--RTC', dest='RTC', type=float, help='RT cutoff for duplicates')
parser.add_argument('--MZC', dest='MZC', type=float, help='M/Z cutoff for duplicates')
parser.add_argument('--scoreC', dest='scoreC', type=float, help='Score threshold for initial curve fit')
parser.add_argument('--ppmC', dest='ppmC', type=float, help='ppm cutoff after curve fit')
parser.add_argument('--maxfev', dest='maxfev', type=int, help='maxfev for curve fit')


args = parser.parse_args()

if args.RTC: RTcutoff = args.RTC
if args.MZC: MZcutoff = args.MZC
if args.scoreC: score_threshold = args.scoreC
if args.ppmC: ppm_threshold = args.ppmC
if args.maxfev: maxfev = args.maxfev
if args.out: outpath = args.out
if args.plots: plots = args.plots

print(args)

# MAKE PLOTS DIRECTORY

if not os.path.isdir(plots):
    os.mkdir(plots)

# FILE IN

data = pd.DataFrame()
for table in args.table:
    data = pd.concat([data, pd.read_csv(table)])

data.dropna(inplace=True, subset=['Int.'])
data = data.loc[data['Tree hierarchy'] != 'no peak']
data = data.loc[data['Meas. M/z'] != '-'] # drop rows with no m/z
data['Score'] = data['Score'].apply(inthelper)
data.sort_values('Score', inplace=True, ascending=False)
print('pooled: ' + str(data.shape))

# PLOTTING SECTION

# hist of initial score distribution
figure, (ax1) = plt.subplots(1, 1, figsize=(10,10))
plt.xlabel('score')
plt.ylabel('n')
plt.hist(data['Score'], log=True, bins=20)
plt.savefig(plots + '/hist-scores-all.png')

# scatterplot of full dataset
figure, (ax1) = plt.subplots(1, 1, figsize=(15,15))

plt.scatter(data['Meas. M/z'].apply(float), data['Dev.(ppm)'].apply(float), c=data['Score'], norm=colors.LogNorm())
plt.colorbar(label='Score')
plt.xlabel('M/z')
plt.ylabel('deviation (ppm)')
plt.ylim(-30,30)
plt.savefig(plots + '/scatter-30ppm.png')

# threshold scores
hsdata = data.loc[data['Score'] > score_threshold]

# hist of thresholded score distribution
figure, (ax1) = plt.subplots(1, 1, figsize=(10,10))
plt.xlabel('score')
plt.ylabel('n')
plt.hist(hsdata['Score'], log=True, bins=20)
plt.savefig(plots + '/hist-scores-cut.png')

# scatterplot of score-thresholded data
figure, (ax1) = plt.subplots(1, 1, figsize=(15,15))
popt, pcov = opt.curve_fit(invfn, hsdata['Meas. M/z'].apply(float), hsdata['Dev.(ppm)'].apply(float), maxfev=maxfev,bounds=([-np.inf, -np.inf, -np.inf], [np.inf, min(data['Meas. M/z'].apply(float))-100, np.inf]))
plt.plot(np.sort(data['Meas. M/z'].apply(float)), invfn(np.sort(data['Meas. M/z'].apply(float)), *popt))

plt.scatter(hsdata['Meas. M/z'].apply(float), hsdata['Dev.(ppm)'].apply(float), c=hsdata['Score'], norm=colors.LogNorm())
plt.colorbar(label='Score')
plt.xlabel('M/z')
plt.ylabel('deviation (ppm)')
plt.ylim(-30,30)
plt.savefig(plots + '/scatter-hscoring.png')

# cut everything further than ppm_threshold from the trendline
cleaned = data.loc[abs(data['Dev.(ppm)'].apply(float) - invfn(data['Meas. M/z'].apply(float), *popt)) < ppm_threshold]
print('ppm threshold: ' + str(cleaned.shape))

# hist of ppm distribution, before & after cut
figure, (ax1) = plt.subplots(1, 1, figsize=(10,10))
plt.xlabel('|delta-ppm|')
plt.ylabel('n')
plt.hist(abs(data['Dev.(ppm)'].apply(float) - invfn(data['Meas. M/z'].apply(float), *popt)), log=True, bins=20)
plt.hist(abs(cleaned['Dev.(ppm)'].apply(float) - invfn(cleaned['Meas. M/z'].apply(float), *popt)), log=True, bins=20)
plt.savefig(plots + '/hist-ppms.png')

# scatterplot of ppm-thresholded data
figure, (ax1) = plt.subplots(1, 1, figsize=(15,15))

plt.plot(np.sort(data['Meas. M/z'].apply(float)), invfn(np.sort(data['Meas. M/z'].apply(float)), *popt))

plt.scatter(cleaned['Meas. M/z'].apply(float), cleaned['Dev.(ppm)'].apply(float), c=cleaned['Score'], norm=colors.LogNorm())
plt.colorbar(label='Score')
plt.xlabel('M/z')
plt.ylabel('deviation (ppm)')
plt.ylim(-30,30)
plt.savefig(plots + '/scatter-7ppm.png')

# DISAMBIGUATION SECTION (REMOVE DUPLICATE PEPTIDES BASED ON SCORE)

cleaned['Start'] = cleaned['Range'].str.split(expand=True)[0]
cleaned['End'] = cleaned['Range'].str.split(expand=True)[2]
cleaned['Start'] = cleaned['Start'].apply(int)
cleaned['End'] = cleaned['End'].apply(int)

grouped = cleaned.groupby(['Sequence', 'z'])
topscores = pd.DataFrame(columns=cleaned.columns)
for name, group in grouped:
    topscoring = group.loc[group['Score'] == max(group['Score'])]
    topscores = pd.concat([topscores, topscoring.loc[topscoring['Int.'] == max(topscoring['Int.'])]])
topscores.reset_index(inplace=True, drop=True)
flagged = pd.DataFrame()
for i in topscores.index:
    r = topscores.loc[i]
    grp = topscores.loc[redflag(topscores, r)]
    r = topscores.loc[topscores.index == i]
    if not grp.empty:
        if r['Score'].iloc[0] <= max(grp['Score']):
            flagged = pd.concat([flagged, r])

print('topscores: ' + str(topscores.shape))
print('flagged: ' + str(flagged.shape))

# scatterplot of flagged data
figure, (ax1) = plt.subplots(1, 1, figsize=(15,15))
plt.plot(np.sort(data['Meas. M/z'].apply(float)), invfn(np.sort(data['Meas. M/z'].apply(float)), *popt))

plt.scatter(flagged['Meas. M/z'].apply(float), flagged['Dev.(ppm)'].apply(float), c=flagged['Score'], norm=colors.LogNorm())
plt.colorbar(label='Score')
plt.xlabel('M/z')
plt.ylabel('deviation (ppm)')
plt.ylim(-30,30)
plt.savefig(plots + '/scatter-flagged.png')

# SORT & SAVE TO CSV

clean = topscores.drop(flagged.index)
clean.sort_values(['Start', 'End'], inplace=True)
clean.to_csv(outpath, index=False)
print('clean: ' + str(clean.shape))

# scatterplot of cleaned data
figure, (ax1) = plt.subplots(1, 1, figsize=(15,15))
plt.plot(np.sort(data['Meas. M/z'].apply(float)), invfn(np.sort(data['Meas. M/z'].apply(float)), *popt))

plt.scatter(clean['Meas. M/z'].apply(float), clean['Dev.(ppm)'].apply(float), c=clean['Score'], norm=colors.LogNorm())
plt.colorbar(label='Score')
plt.xlabel('M/z')
plt.ylabel('deviation (ppm)')
plt.ylim(-30,30)
plt.savefig(plots + '/scatter-clean.png')
