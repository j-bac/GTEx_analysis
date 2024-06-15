import matplotlib.pyplot as plt

def barplot(series,save_path="../figures/number_of_samples_barplot.png",figsize=(12,4)):
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111)
    x = range(len(series))
    ax.bar(x, series.values, width=0.8)
    plt.xticks(x, series.index, rotation=-45, ha='left', fontsize=10)
    ax.set_xlim([x[0]-1,x[-1]+1])
    ax.set_ylabel('Samples', fontsize=12)
    if save_path is not None:
        plt.savefig(save_path,dpi=300,bbox_inches='tight')
    plt.show()