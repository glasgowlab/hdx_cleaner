import warnings
warnings.filterwarnings("ignore")
import pandas as pd
pd.options.mode.chained_assignment = None
import matplotlib
from matplotlib import rcParams
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import scipy.optimize as opt
import argparse
import os

# DEFAULT PARAMS

RTcutoff = 0.5
MZcutoff = 0.1
score_floor = 0
score_threshold = 150
ppm_threshold = 7
fit = 'inv'
maxfev = 800
plots = 'csvplots'
overlap_method = 'select'
skip_res = 2
select_score = 0
ov_threshold = 0.8
ymin=0.9


# COLUMN LABELS

rtlabel = 'Rt(min)'
scorelabel = 'Score'
zlabel = 'z'
mzlabel = 'Meas. M/z'
mrlabel = 'Meas. Mr'
seqlabel = 'Sequence'
rangelabel = 'Range'
startlabel = 'Start'
endlabel = 'End'
ppmlabel = 'Dev.(ppm)'

# HELPER FNS

def inthelper(n):
    if n == '-':
        return 0
    if type(n) is str:
        return int(n.strip())
    return int(n)

def fitfn(x, a, b, c):
    if fit == 'inv':
        return c + a / (x - b)
    else:
        return c + a * x + b * (x ** 2)

def formatfill(data):
    if 0 in data[scorelabel].values:
        return (len(data), np.mean(data[scorelabel]), 100 * data[scorelabel].value_counts()[0] / len(data))
    return (len(data), np.mean(data[scorelabel]), 0.)

def redflag(topscores, r):
    # return true if they're closer in RT and mass than the cutoffs and if they're nonidentical and nonoverlapping
    RTflag = (abs(topscores[rtlabel].apply(float) - float(r[rtlabel])) < RTcutoff)
    MZflag = (abs(topscores[mzlabel].apply(float) - float(r[mzlabel])) < MZcutoff) | (abs(topscores[mrlabel].apply(float) - float(r[mrlabel])) < MZcutoff)
    notsame = (topscores[seqlabel] != r[seqlabel])
    nonoverlapping = (topscores[startlabel].apply(int) + skip_res > int(r[endlabel])) | (int(r[startlabel]) + skip_res > topscores[endlabel].apply(int))
    return RTflag & MZflag & notsame & nonoverlapping

def yellowflag(topscores, r):
    # return true if they're closer in RT and mass than the cutoffs and if they're nonidentical but overlapping
    RTflag = (abs(topscores[rtlabel].apply(float) - float(r[rtlabel])) < RTcutoff)
    MZflag = (abs(topscores[mzlabel].apply(float) - float(r[mzlabel])) < MZcutoff) | (abs(topscores[mrlabel].apply(float) - float(r[mrlabel])) < MZcutoff)
    notsame = (topscores[seqlabel] != r[seqlabel])
    overlapping = (topscores[startlabel].apply(int) + skip_res <= int(r[endlabel])) & (int(r[startlabel]) + skip_res <= topscores[endlabel].apply(int))
    return RTflag & MZflag & notsame & overlapping

def scorehist(data, title, path):
    figure, (ax1) = plt.subplots(1, 1, figsize=(8, 8))
    plt.xlabel(scorelabel)
    plt.ylabel('n')
    plt.ylim(ymin, ymax)
    plt.title(title)
    plt.hist(data[scorelabel], log=True, bins=100, range=[0, xmax])
    #print('% 0: ' + str(100 * data[scorelabel].value_counts()[0] / len(data)))
    plt.tight_layout()
    plt.savefig(path)

def mzscatter(data, path, popt=False):
    rcParams['font.sans-serif'][7] = 'DejaVu Sans'
    rcParams['font.sans-serif'][0] = 'Arial'
    # rcParams['mathtext.fontset'] = 'Arial'
    # plt.rcParams['text.usetex'] = True

    afont = {'size': 6}
    matplotlib.rc('font', **afont)

    pt = 72
    figure, (ax1) = plt.subplots(1, 1, figsize=(109.7 / pt, 109.7 / pt))
    if type(popt) != bool:
        plt.plot(np.sort(data[mzlabel].apply(float)), fitfn(np.sort(data[mzlabel].apply(float)), *popt), linewidth=0.5)
    plt.scatter(data[mzlabel].apply(float), data[ppmlabel].apply(float), c=data[scorelabel]+1,
                norm=colors.LogNorm(vmin=1, vmax=xmax), s=50/pt, linewidth=0.05)
    plt.colorbar(label=scorelabel)
    plt.xlabel('$\it{{m}}$/$\it{{z}}$ (D/e)')
    plt.ylabel('mass error (ppm)')
    plt.ylim(-30, 30)
    plt.tight_layout()
    plt.savefig(path)

# CMD LINE ARGS CODE

parser = argparse.ArgumentParser(description='Biotools CSV cleaner: pool csvs and disambiguate duplicates')
parser.add_argument('--t', '--table', dest='table', help="paths to input csvs", nargs='+', required=True)
parser.add_argument('--o', '--out', dest='out', help="path to output csv")
parser.add_argument('--p', '--plots', dest='plots', help="directory for output plots")
parser.add_argument('--ov', '--overlap', dest='overlap', help="how to treat overlapping duplicates: keep/drop/select?")
parser.add_argument('--sf', '--scorefloor', dest='scorefloor', type=float, help='Score floor for assignments')
parser.add_argument('--r', '--rangeslist', dest='rangeslist', help="destination for rangeslist")
parser.add_argument('--RTC', dest='RTC', type=float, help='RT cutoff for duplicates')
parser.add_argument('--MZC', dest='MZC', type=float, help='M/Z cutoff for duplicates')
parser.add_argument('--fit', dest='fit', help='type of curve fit: inv, or quad for quadratic')
parser.add_argument('--scoreC', dest='scoreC', type=float, help='Score threshold for initial curve fit')
parser.add_argument('--ppmC', dest='ppmC', type=float, help='ppm cutoff after curve fit')
parser.add_argument('--maxfev', dest='maxfev', type=int, help='maxfev for curve fit')

args = parser.parse_args()

if args.RTC: RTcutoff = args.RTC
if args.MZC: MZcutoff = args.MZC
if args.scorefloor: score_floor = args.scorefloor
if args.scoreC: score_threshold = args.scoreC
if args.ppmC: ppm_threshold = args.ppmC
if args.fit: fit = args.fit
if args.maxfev: maxfev = args.maxfev
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

data = data[data[mzlabel] != '-']
data = data[data[ppmlabel] != '-']
data.dropna(subset=[seqlabel], inplace=True)
data[scorelabel].fillna(0, inplace=True)
data[scorelabel] = data[scorelabel].apply(inthelper)

#data = data.loc[data['Tree hierarchy'] != 'no peak']
data.sort_values(scorelabel, inplace=True, ascending=False)

print('pooled: ' + str(data.shape))
print('mean score: ' + str(np.mean(data[scorelabel])))
print('% 0: ' + str(100 * data[scorelabel].value_counts()[0] / len(data)))
print()

ymax=len(data)
xmax = max(data[scorelabel])

plt.rc('font', size=25)

# hist of initial score distribution

scorehist(data, "pooled: %d\n mean score: %1.5f\n %%0: %1.5f" %formatfill(data), plots + '/hist-scores-all.pdf')

# scatterplot of full dataset

mzscatter(data, plots + '/scatter-30ppm.pdf')

# threshold scores
hsdata = data.loc[data[scorelabel] > score_threshold]

# hist of thresholded score distribution

scorehist(hsdata, "high scores: %d\n mean score: %1.5f" %(len(hsdata), np.mean(hsdata[scorelabel])), plots + '/hist-scores-cut.pdf')

print('highscoring: ' + str(hsdata.shape))
print('mean score: ' + str(np.mean(hsdata[scorelabel])))
print()

popt, pcov = opt.curve_fit(fitfn, hsdata[mzlabel].apply(float), hsdata[ppmlabel].apply(float), maxfev=maxfev,
                           bounds=([-np.inf, -np.inf, -np.inf], [np.inf, min(data[mzlabel].apply(float))-100 if fit == 'inv' else np.inf, np.inf]))

# scatterplot of score-thresholded data

mzscatter(hsdata, plots + '/scatter-hscoring.pdf', popt)

# hist of ppm distribution, before & after cut

figure, (ax1) = plt.subplots(1, 1, figsize=(6,6))
plt.xlabel('|delta-ppm|')
plt.ylabel('n')
plt.hist(abs(data[ppmlabel].apply(float) - fitfn(data[mzlabel].apply(float), *popt)), log=True, bins=100)
plt.vlines(ppm_threshold, 0, len(data), linestyles='dashed', colors='k')
plt.savefig(plots + '/hist-ppms.pdf')

print('m/z deviation (30ppm): ' + str(np.sqrt(np.mean((data[ppmlabel].apply(float) - fitfn(data[mzlabel].apply(float), *popt)) ** 2))))
print('m/z deviation (hsdata): ' + str(np.sqrt(np.mean((hsdata[ppmlabel].apply(float) - fitfn(hsdata[mzlabel].apply(float), *popt)) ** 2))))


# cut everything further than ppm_threshold from the trendline
data = data.loc[abs(data[ppmlabel].apply(float) - fitfn(data[mzlabel].apply(float), *popt)) < ppm_threshold]

print('ppm threshold: ' + str(data.shape))
print('mean score: ' + str(np.mean(data[scorelabel])))
print('% 0: ' + str(100 * data[scorelabel].value_counts()[0] / len(data)))
print()

print('m/z deviation (7ppm): ' + str(np.sqrt(np.mean((data[ppmlabel].apply(float) - fitfn(data[mzlabel].apply(float), *popt)) ** 2))))

# PLOTS
mzscatter(data, plots + '/scatter-7ppm.pdf', popt)
scorehist(data, "ppm: %d\n mean score: %1.5f\n %%0: %1.5f" %formatfill(data), plots + '/hist-scores-ppm.pdf')


# DISAMBIGUATION SECTION (REMOVE DUPLICATE PEPTIDES BASED ON SCORE)

if rangelabel in data.columns:
    data[startlabel] = data[rangelabel].str.split(expand=True)[0]
    data[endlabel] = data[rangelabel].str.split(expand=True)[2]
data[startlabel] = data[startlabel].apply(int)
data[endlabel] = data[endlabel].apply(int)

grouped = data.groupby([seqlabel, zlabel])
topscores = pd.DataFrame(columns=data.columns)
ties = 0
for name, group in grouped:
    topscoring = group.loc[group[scorelabel] == max(group[scorelabel])]
    if(len(topscoring[scorelabel]) > 1):
        ties = ties + 1
    topscores = pd.concat([topscores, topscoring.loc[topscoring['Int.'] == max(topscoring['Int.'])]])
topscores.reset_index(inplace=True, drop=True)

print('ties: ' + str(ties))
print('topscores: ' + str(topscores.shape))
print('mean score: ' + str(np.mean(topscores[scorelabel])))
print('% 0: ' + str(100 * topscores[scorelabel].value_counts()[0] / len(topscores)))
print('m/z deviation (topscores): ' + str(np.sqrt(np.mean((topscores[ppmlabel].apply(float) - fitfn(topscores[mzlabel].apply(float), *popt)) ** 2))))
print()

# PLOTS
scorehist(topscores, "topscores: %d\n mean score: %1.5f\n %%0: %1.5f" %formatfill(topscores), plots + '/hist-scores-topscores.pdf')
mzscatter(topscores, plots + '/scatter-topscore.pdf', popt)

if score_floor:
    topscores = topscores.loc[topscores[scorelabel] >= score_floor]

    # PLOTS
    scorehist(topscores, "score-floor: %d\n mean score: %1.5f\n %%0: %1.5f" %formatfill(topscores), plots + '/hist-scores-scorefloor.pdf')
    mzscatter(topscores, plots + '/scatter-scorefloor.pdf', popt)

    print('scorefloor: ' + str(topscores.shape))
    print('mean score: ' + str(np.mean(topscores[scorelabel])))
    print()

topscores.sort_values(['Start', 'End'], inplace=True)
topscores.reset_index(inplace=True, drop=True)

flagged = pd.DataFrame()
yellow = pd.DataFrame()
ties = 0
overlaps = 0
for i in topscores.index:
    r = topscores.loc[i]
    grp = topscores.loc[redflag(topscores, r)]
    overlap = topscores.loc[yellowflag(topscores, r)]
    r = topscores.loc[topscores.index == i]
    # if not grp.empty:
    #     if r[scorelabel].iloc[0] == max(grp[scorelabel]):
    #         ties = ties + 1
    #     if r[scorelabel].iloc[0] <= max(grp[scorelabel]):
    #         flagged = pd.concat([flagged, r])
    #         continue
    # if not overlap.empty:
    #     if overlap_method == 'drop' or overlap_method != 'keep' and r[scorelabel].iloc[0] <= max(overlap[scorelabel]):
    #         yellow = pd.concat([yellow, r])
    #     overlaps = overlaps + 1

    if not grp.empty:
        if r[scorelabel].iloc[0] == max(grp[scorelabel]):
            ties = ties + 1
        if r[scorelabel].iloc[0] <= select_score or (overlap_method != 'keep'): # max(grp[scorelabel]) > select_score and
            print('redflag')
            print(r[[mzlabel, mrlabel, scorelabel, startlabel, endlabel, seqlabel]])
            print(grp[[mzlabel, mrlabel, scorelabel, startlabel, endlabel, seqlabel]])
            print('\n\n\n')
            flagged = pd.concat([flagged, r])
            continue
    if not overlap.empty:
        lens_overlap = []
        for j in overlap.index:
            ovrow = overlap.loc[j]
            lens_overlap.append(min(np.abs(int(ovrow[startlabel]) + skip_res - int(r[endlabel])), np.abs(
                    int(r[startlabel]) + skip_res - int(ovrow[endlabel]))))
        ovfrac = min(lens_overlap) / (int(ovrow[endlabel]) - (int(ovrow[startlabel]) + skip_res))
        if r[scorelabel].iloc[0] <= select_score or (max(overlap[scorelabel]) > select_score and (overlap_method == 'drop' or overlap_method != 'keep' and ovfrac < ov_threshold)):
            print('yellowflag')
            print(r[[mzlabel, mrlabel, scorelabel, startlabel, endlabel, seqlabel]])
            print(overlap[[mzlabel, mrlabel, scorelabel, startlabel, endlabel, seqlabel]])
            print('\n\n\n')
            yellow = pd.concat([yellow, r])
        overlaps = overlaps + 1

print('ties: ' + str(ties))
print('overlaps: ' + str(overlaps))

# HANDLE REDFLAG PEPTIDES

if not flagged.empty:

    print('flagged: ' + str(flagged.shape))
    print('mean score: ' + str(np.mean(flagged[scorelabel])))
    if 0 in flagged[scorelabel].values:
        print('% 0: ' + str(100 * flagged[scorelabel].value_counts()[0] / len(flagged)))
    else:
        print('% 0: ' + str(0))
    print('m/z deviation (flagged): ' + str(np.sqrt(np.mean((flagged[ppmlabel].apply(float) - fitfn(flagged[mzlabel].apply(float), *popt)) ** 2))))
    print()

    # PLOTS
    mzscatter(flagged, plots + '/scatter-flagged.pdf', popt)
    scorehist(flagged, "flagged: %d\n mean score: %1.5f\n %%0: %1.5f" % formatfill(flagged),
              plots + '/hist-scores-flagged.pdf')

# HANDLE YELLOW FLAG PEPTIDES

if not yellow.empty:

    print('yellow: ' + str(yellow.shape))
    print('mean score: ' + str(np.mean(yellow[scorelabel])))
    if 0 in yellow[scorelabel].values:
        print('% 0: ' + str(100 * yellow[scorelabel].value_counts()[0] / len(yellow)))
    else:
        print('% 0: ' + str(0))
    print('m/z deviation (yellow): ' + str(np.sqrt(np.mean((yellow[ppmlabel].apply(float) - fitfn(yellow[mzlabel].apply(float), *popt)) ** 2))))

    print()

    # PLOTS
    mzscatter(yellow, plots + '/scatter-yellow.pdf', popt)
    scorehist(yellow, "yellow: %d\n mean score: %1.5f\n %%0: %1.5f" % formatfill(yellow),
              plots + '/hist-scores-yellow.pdf')

clean = topscores.drop(flagged.index.union(yellow.index))

# SORT & SAVE TO CSV

clean.sort_values([startlabel, endlabel], inplace=True)
if args.out:
    clean.to_csv(args.out, index=False)
if args.rangeslist:
    clean[[startlabel, endlabel]].to_csv(args.rangeslist, index=False)

print('clean: ' + str(clean.shape))
print('mean score: ' + str(np.mean(clean[scorelabel])))
if not yellow.empty:
    print('ambiguous: ' + str(overlaps - yellow.shape[0]))
if 0 in clean[scorelabel].values:
    print('% 0: ' + str(100 * clean[scorelabel].value_counts()[0] / len(clean)))
else:
    print('% 0: ' + str(0))
print('m/z deviation (clean): ' + str(np.sqrt(np.mean((clean[ppmlabel].apply(float) - fitfn(clean[mzlabel].apply(float), *popt)) ** 2))))

# PLOTS
mzscatter(clean, plots + '/scatter-clean.pdf', popt)
scorehist(clean, "clean: %d\n mean score: %1.5f\n %%0: %1.5f" %formatfill(clean), plots + '/hist-scores-clean.pdf')