import pandas as pd 
import matplotlib.pyplot as plt
import os 

def ms_data_plot_file(filename: str):
    num_plots: int = 0
    files = []

    for (dirpath, dirnames, filenames) in os.walk(filename):
        num_plots = len(filenames)
        print(filenames)
        files.extend(filenames)
        break

    data = []
    for file in files:
        print(file)
        frame = pd.read_csv(filename+"/" + file)
        print(frame)
        data.append(pd.read_csv(filename+"/" + file))

    fig, axes = plt.subplots(num_plots)
    i: int = 0
    for timepoint in data:
        axes[i].plot(timepoint.iloc[:, 0], timepoint.iloc[:, 1])
        i = i + 1
    fig.tight_layout(pad=5.0)
    plt.show()

ms_data_plot_file('/Users/pranav/programming/PIGEON/example/0012-0017-YAGVSY')

