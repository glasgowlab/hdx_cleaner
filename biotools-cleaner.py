import matplotlib.colors
import pandas as pd
pd.options.mode.chained_assignment = None
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import scipy.optimize as opt
import argparse
import os

# DEFAULT PARAMS

RTcutoff = 0.5
MZcutoff = 0.1
score_threshold = 150
ppm_threshold = 7
maxfev = 800
plots = 'csvplots'
outpath = 'CLEAN.csv'
overlap_method = 'select'
skip_res = 2
ymin=0.9

# HELPER FNS

def inthelper(n):
    if n == '-':
        return 0
    return int(n.strip())

def invfn(x, a, b, c):
    return c + a / (x - b)

def formatfill(data):
    return (len(data), np.mean(data['Score']), 100 * data['Score'].value_counts()[0] / len(data))

def redflag(topscores, r):
    # return true if they're closer in RT and mass than the cutoffs and if they're nonidentical and nonoverlapping
    RTflag = (abs(topscores['Rt(min)'].apply(float) - float(r['Rt(min)'])) < RTcutoff)
    MZflag = (abs(topscores['Meas. M/z'].apply(float) - float(r['Meas. M/z'])) < MZcutoff)
    notsame = (topscores['Sequence'] != r['Sequence'])
    nonoverlapping = (topscores['Start'].apply(int) + skip_res > int(r['End'])) | (int(r['Start']) + skip_res > topscores['End'].apply(int))
    return RTflag & MZflag & notsame & nonoverlapping

def yellowflag(topscores, r):
    # return true if they're closer in RT and mass than the cutoffs and if they're nonidentical but overlapping
    RTflag = (abs(topscores['Rt(min)'].apply(float) - float(r['Rt(min)'])) < RTcutoff)
    MZflag = (abs(topscores['Meas. M/z'].apply(float) - float(r['Meas. M/z'])) < MZcutoff)
    notsame = (topscores['Sequence'] != r['Sequence'])
    overlapping = (topscores['Start'].apply(int) + skip_res <= int(r['End'])) | (int(r['Start']) + skip_res <= topscores['End'].apply(int))
    return RTflag & MZflag & notsame & overlapping

# CMD LINE ARGS CODE

parser = argparse.ArgumentParser(description='Biotools CSV cleaner: pool csvs and disambiguate duplicates')
parser.add_argument('--t', '--table', dest='table', help="paths to input csvs", nargs='+', required=True)
parser.add_argument('--o', '--out', dest='out', help="path to output csv")
parser.add_argument('--p', '--plots', dest='plots', help="directory for output plots")
parser.add_argument('--ov', '--overlap', dest='overlap', help="how to treat overlapping duplicates: keep/drop/select?")
parser.add_argument('--r', '--rangeslist', dest='rangeslist', help="destination for rangeslist")
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
if args.overlap: overlap_method = args.overlap
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
data['Score'] = data['Score'].apply(inthelper)
data.sort_values('Score', inplace=True, ascending=False)
print('pooled: ' + str(data.shape))
print('mean score: ' + str(np.mean(data['Score'])))
ymax=len(data)
xmax = max(data['Score'])

# TEMP scatterplot of full dataset

def floathelper(n):
    if n == '-':
        return 0
    return float(n.strip())
figure, (ax1) = plt.subplots(1, 1, figsize=(5,5))
data['MSMS Cov.(%)'] = data['MSMS Cov.(%)'].apply(floathelper)
plt.yscale("log")
plt.scatter(data['MSMS Cov.(%)'], data['Score'].apply(float)+1)
plt.savefig(plots + '/scatter-score-coverage.png')




# PLOTTING SECTION

# hist of initial score distribution
figure, (ax1) = plt.subplots(1, 1, figsize=(6,6))
plt.xlabel('score')
plt.ylabel('n')
plt.ylim(ymin,ymax)
plt.title("pooled: %d\n mean score: %1.5f\n %%0: %1.5f" %formatfill(data))
plt.hist(data['Score'], log=True, bins=100, range=[0,xmax])
print('% 0: ' + str(100 * data['Score'].value_counts()[0] / len(data)))
plt.savefig(plots + '/hist-scores-all.png')

# scatterplot of full dataset
figure, (ax1) = plt.subplots(1, 1, figsize=(5,5))

plt.scatter(data['Meas. M/z'].apply(float), data['Dev.(ppm)'].apply(float), c=data['Score'], norm=colors.LogNorm(vmin=1, vmax=xmax))
plt.colorbar(label='Score')
plt.xlabel('M/z')
plt.ylabel('deviation (ppm)')
plt.ylim(-30,30)
plt.savefig(plots + '/scatter-30ppm.png')

print()

# threshold scores
hsdata = data.loc[data['Score'] > score_threshold]

# hist of thresholded score distribution
figure, (ax1) = plt.subplots(1, 1, figsize=(6,6))
plt.xlabel('score')
plt.ylabel('n')
plt.ylim(ymin,ymax)
plt.title("high scores: %d\n mean score: %1.5f" %(len(hsdata), np.mean(hsdata['Score'])))
plt.hist(hsdata['Score'],log=True, bins=100, range=[0,xmax])
print('highscoring: ' + str(hsdata.shape))
print('mean score: ' + str(np.mean(hsdata['Score'])))
plt.savefig(plots + '/hist-scores-cut.png')

# scatterplot of score-thresholded data
figure, (ax1) = plt.subplots(1, 1, figsize=(5,5))
popt, pcov = opt.curve_fit(invfn, hsdata['Meas. M/z'].apply(float), hsdata['Dev.(ppm)'].apply(float), maxfev=maxfev,bounds=([-np.inf, -np.inf, -np.inf], [np.inf, min(data['Meas. M/z'].apply(float))-100, np.inf]))
plt.plot(np.sort(data['Meas. M/z'].apply(float)), invfn(np.sort(data['Meas. M/z'].apply(float)), *popt))

plt.scatter(hsdata['Meas. M/z'].apply(float), hsdata['Dev.(ppm)'].apply(float), c=hsdata['Score'], norm=colors.LogNorm(vmin=1, vmax=xmax))
plt.colorbar(label='Score')
plt.xlabel('M/z')
plt.ylabel('deviation (ppm)')
plt.ylim(-30,30)
plt.savefig(plots + '/scatter-hscoring.png')

print()

# cut everything further than ppm_threshold from the trendline
cleaned = data.loc[abs(data['Dev.(ppm)'].apply(float) - invfn(data['Meas. M/z'].apply(float), *popt)) < ppm_threshold]
print('ppm threshold: ' + str(cleaned.shape))
print('mean score: ' + str(np.mean(cleaned['Score'])))

# hist of ppm distribution, before & after cut
figure, (ax1) = plt.subplots(1, 1, figsize=(6,6))
plt.xlabel('|delta-ppm|')
plt.ylabel('n')
plt.hist(abs(data['Dev.(ppm)'].apply(float) - invfn(data['Meas. M/z'].apply(float), *popt)), log=True, bins=20)
plt.hist(abs(cleaned['Dev.(ppm)'].apply(float) - invfn(cleaned['Meas. M/z'].apply(float), *popt)), log=True, bins=20)
plt.savefig(plots + '/hist-ppms.png')

# scatterplot of ppm-thresholded data
figure, (ax1) = plt.subplots(1, 1, figsize=(5,5))

plt.plot(np.sort(data['Meas. M/z'].apply(float)), invfn(np.sort(data['Meas. M/z'].apply(float)), *popt))

plt.scatter(cleaned['Meas. M/z'].apply(float), cleaned['Dev.(ppm)'].apply(float), c=cleaned['Score'], norm=colors.LogNorm(vmin=1, vmax=xmax))
plt.colorbar(label='Score')
plt.xlabel('M/z')
plt.ylabel('deviation (ppm)')
plt.ylim(-30,30)
plt.savefig(plots + '/scatter-7ppm.png')

# hist of ppm-thresholded score distribution
figure, (ax1) = plt.subplots(1, 1, figsize=(6,6))
plt.xlabel('score')
plt.ylabel('n')
plt.ylim(ymin,ymax)
plt.title("ppm: %d\n mean score: %1.5f\n %%0: %1.5f" %formatfill(cleaned))
plt.hist(cleaned['Score'],log=True, bins=100, range=[0,xmax])
print('% 0: ' + str(100 * cleaned['Score'].value_counts()[0] / len(cleaned)))
plt.savefig(plots + '/hist-scores-ppm.png')

print()

# DISAMBIGUATION SECTION (REMOVE DUPLICATE PEPTIDES BASED ON SCORE)

cleaned['Start'] = cleaned['Range'].str.split(expand=True)[0]
cleaned['End'] = cleaned['Range'].str.split(expand=True)[2]
cleaned['Start'] = cleaned['Start'].apply(int)
cleaned['End'] = cleaned['End'].apply(int)

grouped = cleaned.groupby(['Sequence', 'z'])
topscores = pd.DataFrame(columns=cleaned.columns)
ties = 0
for name, group in grouped:
    topscoring = group.loc[group['Score'] == max(group['Score'])]
    if(len(topscoring['Score']) > 1):
        ties = ties + 1
    topscores = pd.concat([topscores, topscoring.loc[topscoring['Int.'] == max(topscoring['Int.'])]])
topscores.reset_index(inplace=True, drop=True)

print('ties: ' + str(ties))

print('topscores: ' + str(topscores.shape))
print('mean score: ' + str(np.mean(topscores['Score'])))

# hist of topscore score distribution
figure, (ax1) = plt.subplots(1, 1, figsize=(6,6))
plt.xlabel('score')
plt.ylabel('n')
plt.ylim(ymin,ymax)
plt.title("topscores: %d\n mean score: %1.5f\n %%0: %1.5f" %formatfill(topscores))
plt.hist(topscores['Score'],log=True, bins=100, range=[0,xmax])
print('% 0: ' + str(100 * topscores['Score'].value_counts()[0] / len(topscores)))
plt.savefig(plots + '/hist-scores-topscores.png')

print()

flagged = pd.DataFrame()
yellow = pd.DataFrame()
ties = 0
overlaps = 0
for i in topscores.index:
    r = topscores.loc[i]
    grp = topscores.loc[redflag(topscores, r)]
    overlap = topscores.loc[yellowflag(topscores, r)]
    r = topscores.loc[topscores.index == i]
    if not grp.empty:
        if r['Score'].iloc[0] == max(grp['Score']):
            ties = ties + 1
        if r['Score'].iloc[0] <= max(grp['Score']):
            flagged = pd.concat([flagged, r])
            continue
    if not overlap.empty:
        if overlap_method == 'drop' or overlap_method != 'keep' and r['Score'].iloc[0] <= max(overlap['Score']):
            yellow = pd.concat([yellow, r])
        overlaps = overlaps + 1



print('ties: ' + str(ties))
print('overlaps: ' + str(overlaps))

# HANDLE REDFLAG PEPTIDES

if not flagged.empty:

    print('flagged: ' + str(flagged.shape))
    print('mean score: ' + str(np.mean(flagged['Score'])))


    # scatterplot of flagged data
    figure, (ax1) = plt.subplots(1, 1, figsize=(5,5))
    plt.plot(np.sort(data['Meas. M/z'].apply(float)), invfn(np.sort(data['Meas. M/z'].apply(float)), *popt))
    plt.scatter(flagged['Meas. M/z'].apply(float), flagged['Dev.(ppm)'].apply(float), c=flagged['Score'], norm=colors.LogNorm(vmin=1, vmax=xmax))
    plt.colorbar(label='Score')
    plt.xlabel('M/z')
    plt.ylabel('deviation (ppm)')
    plt.ylim(-30,30)
    plt.savefig(plots + '/scatter-flagged.png')

    # hist of flagged score distribution
    figure, (ax1) = plt.subplots(1, 1, figsize=(6,6))
    plt.xlabel('score')
    plt.ylabel('n')
    plt.ylim(ymin, ymax)
    plt.title("flagged: %d\n mean score: %1.5f\n %%0: %1.5f" % formatfill(flagged))
    plt.hist(flagged['Score'],log=True, bins=100, range=[0,xmax])
    print('% 0: ' + str(100 * flagged['Score'].value_counts()[0] / len(flagged)))
    print()
    plt.savefig(plots + '/hist-scores-flagged.png')

# HANDLE YELLOW FLAG PEPTIDES

if not yellow.empty:

    print('yellow: ' + str(yellow.shape))
    print('mean score: ' + str(np.mean(yellow['Score'])))


    # scatterplot of yellow flag data
    figure, (ax1) = plt.subplots(1, 1, figsize=(5,5))
    plt.plot(np.sort(data['Meas. M/z'].apply(float)), invfn(np.sort(data['Meas. M/z'].apply(float)), *popt))
    plt.scatter(yellow['Meas. M/z'].apply(float), yellow['Dev.(ppm)'].apply(float), c=yellow['Score'], norm=colors.LogNorm(vmin=1, vmax=xmax))
    plt.colorbar(label='Score')
    plt.xlabel('M/z')
    plt.ylabel('deviation (ppm)')
    plt.ylim(-30,30)
    plt.savefig(plots + '/scatter-yellow.png')

    # hist of yellow flag score distribution
    figure, (ax1) = plt.subplots(1, 1, figsize=(6,6))
    plt.xlabel('score')
    plt.ylabel('n')
    plt.ylim(ymin, ymax)
    plt.title("yellow: %d\n mean score: %1.5f\n %%0: %1.5f" % formatfill(yellow))
    plt.hist(flagged['Score'],log=True, bins=100, range=[0,xmax])
    print('% 0: ' + str(100 * yellow['Score'].value_counts()[0] / len(yellow)))
    print()
    plt.savefig(plots + '/hist-scores-yellow.png')

clean = topscores.drop(flagged.index.union(yellow.index))

# SORT & SAVE TO CSV

clean.sort_values(['Start', 'End'], inplace=True)
clean.to_csv(outpath, index=False)
if args.rangeslist:
    clean[['Start', 'End']].to_csv(args.rangeslist, index=False)
print('clean: ' + str(clean.shape))
print('mean score: ' + str(np.mean(clean['Score'])))

# scatterplot of cleaned data
figure, (ax1) = plt.subplots(1, 1, figsize=(5,5))
plt.plot(np.sort(data['Meas. M/z'].apply(float)), invfn(np.sort(data['Meas. M/z'].apply(float)), *popt))

plt.scatter(clean['Meas. M/z'].apply(float), clean['Dev.(ppm)'].apply(float), c=clean['Score'], norm=colors.LogNorm(vmin=1, vmax=xmax))
plt.colorbar(label='Score')
plt.xlabel('M/z')
plt.ylabel('deviation (ppm)')
plt.ylim(-30,30)
plt.savefig(plots + '/scatter-clean.png')

# hist of cleaned score distribution
figure, (ax1) = plt.subplots(1, 1, figsize=(6,6))
plt.xlabel('score')
plt.ylabel('n')
plt.ylim(ymin,ymax)
plt.title("clean: %d\n mean score: %1.5f\n %%0: %1.5f" %formatfill(clean))
plt.hist(clean['Score'],log=True, bins=100, range=[0,xmax])
print('% 0: ' + str(100 * clean['Score'].value_counts()[0] / len(clean)))
plt.savefig(plots + '/hist-scores-clean.png')