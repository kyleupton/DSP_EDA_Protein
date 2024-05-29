import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


titleFont = {'fontsize': 32}
# ,
#  'fontweight': rcParams['axes.titleweight'],
#  'color': rcParams['axes.titlecolor'],
#  'verticalalignment': 'baseline',
#  'horizontalalignment': loc}

labelFont = {'fontsize': 20}
# ,
#  'fontweight': rcParams['axes.titleweight'],
#  'color': rcParams['axes.titlecolor'],
#  'verticalalignment': 'baseline',
#  'horizontalalignment': loc}

namedColourList = ['#377eb8', '#ff7f00', '#4daf4a',
                  '#f781bf', '#a65628', '#984ea3',
                  '#999999', '#e41a1c', '#dede00']



def plot_SA_Hist(surfArea):
    
        
    fig, ax = plt.subplots(figsize=(10,5))
    ax.hist(surfArea, bins=50)


    ax.set_xlabel('AOI Surface Area µm', fontdict=labelFont)
    ax.set_ylabel('Count', fontdict=labelFont)
    ax.set_title('AOI Surface Area distribution', fontdict=titleFont)

    print('Min SA')
    print(min(surfArea))

    print('Max SA')
    print(max(surfArea))
    
    return(fig)


# Plot log2 transformed raw data before any normalisation
# def draw_probe_plot_2(dataRaw, dataSortedRaw, namedColourList, sampleInfo, selectedInfo, subSelection=None, title='Title', exp=False, violin=False):
def draw_probe_plot(dataRaw, namedColourList, sampleInfo, selectedInfo, subSelection=None, title='Title', exp=False, violin=False):

    if not (type(subSelection) == 'NoneType'):
        selectedInfo = selectedInfo.loc[subSelection]
    if (type(selectedInfo) == pd.core.series.Series):
        selectedInfo = pd.DataFrame(selectedInfo).T
    
    # Sort data according to dataSortedRaw
    dataRaw = dataRaw.drop(labels=['mean','probeClass'], axis=1).reindex(labels=dataRaw.index)
    # dataRaw = dataRaw[dataSortedRaw.drop(labels=['mean','probeClass'], axis=1).columns]
    # sampleInfo = sampleInfo[dataRaw.drop(labels=['mean','probeClass'], axis=1).columns]
    selectedInfo = selectedInfo[dataRaw.columns]
    # print('sampleInfo.columns')
    # print(sampleInfo.columns)
    # print('selectedInfo')
    # print(selectedInfo)

    fig, ax = plt.subplots(figsize=(15,8))
    
    if exp:
        # ax.boxplot(np.exp2(dataRaw.drop(labels=['mean','probeClass'], axis=1).reindex(labels=dataSortedRaw.index).T) -1, sym='-', labels=dataSortedRaw.index)
        ax.boxplot(np.exp2(dataRaw.T) -1, sym='-', labels=dataRaw.index)
    else:
        # ax.boxplot(dataRaw.drop(labels=['mean','probeClass'], axis=1).reindex(labels=dataSortedRaw.index).T, sym='-', labels=dataSortedRaw.index)
        ax.boxplot(dataRaw.T, sym='-', labels=dataRaw.index)

    if violin:
        if exp:
            # ax.violinplot(np.exp2(dataRaw.drop(labels=['mean','probeClass'], axis=1).reindex(labels=dataSortedRaw.index).T) -1)
            ax.violinplot(np.exp2(dataRaw.T) -1)
        else:
            # ax.violinplot(dataRaw.drop(labels=['mean','probeClass'], axis=1).reindex(labels=dataSortedRaw.index).T)
            ax.violinplot(dataRaw.T)
    else:
        sampleInfo, my_cmap, colours = get_colour_mapping(sampleInfo, selectedInfo)
        my_cmap = plt.get_cmap("nipy_spectral")(colours)
        # print('colours')
        # print(colours)
        # print('my_cmap')
        # print(my_cmap)
        for i,j in enumerate(dataRaw.index):
            # y = dataRaw.drop(labels=['mean','probeClass'], axis=1).loc[j]
            y = dataRaw.loc[j]
            y = y
            if exp:
                y = np.exp2(y.values)-1
            else:
                y = y.values
            x = np.random.normal(i+1, 0.1, len(y))
            for i in range(len(x)): 
                # ax.plot(x[i], y[i], color=my_cmap[i], marker='.', alpha=0.25)
                ax.plot(x[i], y[i], color=my_cmap[i], marker='.')#

    ax.set_xticks(np.arange(1,len(dataRaw.index)+1,1))
    ax.set_xlabel=list(dataRaw.index)
    # print(len(np.arange(0,len(dataSortedRaw.index),1)))
    # print(len(list(dataSortedRaw.index)))
    ax.tick_params(axis='x', labelrotation = 90)

    if exp:
        ax.semilogy()
        ax.set_title(title + ' (untransformed)', size=36)
        ax.set_ylabel('Probe value', size=24)
    else:
        ax.set_title(title + ' (Log2 transformed)', size=36)
        ax.set_ylabel('Log2 probe value', size=24)
#     plt.show()
    
    return(fig)




def probe_GeoMean_Plots(plotData, title=''):
    rows=1
    cols=2
    colours = [namedColourList[2] if x.split('_')[-1] == 'Tumour' else namedColourList[5] if x.split('_')[-1] == 'TME' else namedColourList[1] for x in plotData.index]

    fig,ax = plt.subplots(rows,cols, sharey=True, gridspec_kw={'width_ratios': [4,1]}, figsize=(15,5))
    ax[0].bar(np.linspace(1,len(plotData),len(plotData)), plotData, color=colours)
    ax[1].hist(plotData, bins=int(len(plotData)/10),orientation='horizontal', color='k')
    ax[0].set_xlim(0,len(plotData))
    
    ax[0].text(2,max(plotData)*.95,'Tumour', size=20, c=namedColourList[2])
    ax[0].text(2,max(plotData)*.825,'TME', size=20, c=namedColourList[5])
    ax[0].text(2,max(plotData)*.7,'Other', size=20, c=namedColourList[1])

    fig.suptitle(title, size=36)
    ax[0].set_ylabel('Probe Value', size=18)
    ax[0].set_xlabel('Probes', size=18)
    ax[1].set_xlabel('Count', size=18)

    fig.tight_layout()



class threshold_probes:
    def __init__(self, data, bins):
        
        self.data = data.drop(labels=['mean','probeClass'], axis=1)
        self.bins = bins
        self.thisHist = plt.hist(self.data.values.flatten(), bins = self.bins)

    def zoom_plot(self, start, end):

        maxY = max(self.thisHist[0][2:])*1.2
        # print(maxY)
        plt.hist(self.data.values.flatten(), bins = self.bins)
        plt.xlim(start,end)        
        plt.ylim(0,maxY)        
        
    def check_threshold(self, start, end):
        print(self.thisHist[0][start:end])
        print(self.thisHist[1][start:end])

    def set_threshold_idx(self, idx):
        print(self.thisHist[0][idx])
        print(self.thisHist[1][idx])
        
        self.threshold_idx = idx
        self.threshold = self.thisHist[1][idx]

    def get_filter(self):
        self.ETfilter = self.data >= self.threshold
        return(self.ETfilter)



def get_colour_mapping(sampleInfoExternal, selectedInfo):

    # Get parameters of unique combinations for colour mapping space
    comboUniques = []
    comboColourDictRev = {}
    for c in selectedInfo.columns:
        thisCol = selectedInfo[c]
        combined = '_'.join(thisCol.values)
        comboUniques.append(combined)
        comboColourDictRev[c] = combined
    comboUniques = sorted(list(set(comboUniques)))
    print('\nNumber of unique combinations: {}'.format(len(comboUniques)))
    # print(comboColourDictRev)
    gradient = np.linspace(0, 1, len(comboUniques))
    gradDict = dict(zip(comboUniques,gradient))
    
    # sampleInfoExternal.sort_values(by=['Plate', 'Col', 'Row'], axis=1, inplace=True)
    # sampleInfoExternal.sort_values(by=list(selectedInfo.index), axis=1, inplace=True)
    # Binding Density plot:
    plt.figure(figsize=(40,10))
    my_cmap = plt.get_cmap("nipy_spectral")
    
    colours = []
    for c in sampleInfoExternal.columns:
        # When there are only 2 different values the nipy_spectral cmap outputs black and grey. A cmap like rainbow may be better, or here we change the values away from the extremes of the cmap.
        if len(comboUniques) == 2:
            colours.append(abs(gradDict[comboColourDictRev[c]] - 0.15))
        else:
            colours.append(gradDict[comboColourDictRev[c]])
    # print('selectedInfo.index')
    # print(list(selectedInfo.index))
    # print('selectedInfo.columns')
    # print(list(selectedInfo.columns))
    # print()
    return sampleInfoExternal, my_cmap, colours



def binding_density_plot(sampleInfoExternal, selectedInfo, subSelection):
    # print('selectedInfo')
    # print(selectedInfo)
    if not (subSelection == None):
        selectedInfo = selectedInfo.loc[subSelection]
    if (type(selectedInfo) == pd.core.series.Series):
        selectedInfo = pd.DataFrame(selectedInfo).T
        
    sampleInfoExternal, my_cmap, colours = get_colour_mapping(sampleInfoExternal, selectedInfo)
    sampleInfoExternal.sort_values(by=['Plate', 'Col', 'Row'], axis=1, inplace=True)

    fig, ax = plt.subplots(figsize=(20,5))
    
    bar = ax.bar(sampleInfoExternal.columns,
            sampleInfoExternal.loc['BindingDensity'].values.astype(np.float32), 
            color=my_cmap(colours)
           )#, bottom=0)
    ax.set_title('_'.join(selectedInfo.index))
    ax.set_xticklabels(sampleInfoExternal.columns, rotation='vertical')
    plt.show()

# ToDo: Add legend