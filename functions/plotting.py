import os
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
    print(f'Min SA : {min(surfArea)}')
    print(f'Max SA : {max(surfArea)}')
    return(fig)


# Plot log2 transformed raw data before any normalisation
def draw_probe_plot(probeData, sampleInfo, selectedInfo, subSelection=None, title='Title', exp=False, violin=False):
    if not (type(subSelection) == 'NoneType'):
        selectedInfo = selectedInfo.loc[subSelection]
    if (type(selectedInfo) == pd.core.series.Series):
        selectedInfo = pd.DataFrame(selectedInfo).T
    
    # Sort data according to dataSortedRaw
    probeData = probeData.drop(labels=['mean','probeClass'], axis=1).reindex(labels=probeData.index)
    sampleInfo = sampleInfo.loc[:,probeData.columns]
    selectedInfo = selectedInfo[probeData.columns]

    fig, ax = plt.subplots(figsize=(15,8))
    if exp:
        ax.boxplot(np.exp2(probeData.T) -1, sym='-', labels=probeData.index)
    else:
        ax.boxplot(probeData.T, sym='-', labels=probeData.index)

    if violin:
        if exp:
            ax.violinplot(np.exp2(probeData.T) -1)
        else:
            ax.violinplot(probeData.T)
    else:
        sampleInfo, my_cmap, colours = get_colour_mapping(sampleInfo, selectedInfo)
        my_cmap = plt.get_cmap("nipy_spectral")(colours)
        for i,j in enumerate(probeData.index):
            y = probeData.loc[j]
            y = y
            if exp:
                y = np.exp2(y.values)-1
            else:
                y = y.values
            x = np.random.normal(i+1, 0.1, len(y))
            for i in range(len(x)): 
                ax.plot(x[i], y[i], color=my_cmap[i], marker='.')#

    ax.set_xticks(np.arange(1,len(probeData.index)+1,1))
    ax.set_xlabel=list(probeData.index)
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
        plt.title('Thresholding plot')
        plt.xlabel('Probe value (log2 transformed)')
        plt.ylabel('Count')

    def zoom_plot(self, start, end):
        maxY = max(self.thisHist[0][2:])*1.2
        plt.hist(self.data.values.flatten(), bins = self.bins)
        plt.xlim(start,end)        
        plt.ylim(0,maxY)        
        plt.title('Thresholding Zoom plot')
        plt.xlabel('Probe value (log2 transformed)')
        plt.ylabel('Count')
        
    def check_threshold(self, start, end):
        print(f'Histogram[0] range : \t{self.thisHist[0][start:end]}')
        print(f'Histogram[1] range : \tself.thisHist[1][start:end]}')

    def set_threshold_idx(self, idx):
        print(f'Histogram[0][idx] value : \t{self.thisHist[0][idx]}')
        print(f'Histogram[1][idx] value : \t{self.thisHist[1][idx]}')
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
        
    sampleInfoExternal.sort_values(by=['Plate', 'Col', 'Row'], axis=1, inplace=True)
    sampleInfoExternal, my_cmap, colours = get_colour_mapping(sampleInfoExternal, selectedInfo)

    fig, ax = plt.subplots(figsize=(20,5))
    
    bar = ax.bar(sampleInfoExternal.columns,
            sampleInfoExternal.loc['BindingDensity'].values.astype(np.float32), 
            color=my_cmap(colours)
           )#, bottom=0)
    ax.set_title('_'.join(selectedInfo.index))
    ax.set_xticklabels(sampleInfoExternal.columns, rotation='vertical')
    plt.show()

# ToDo: Add legend

def volcanoPlot(dataPath, file, pVal=True, plot=True):
    if pVal:
        sigType = 'PValue'
    else:
        sigType = 'FDR'

    sigGenes = []
    data = pd.read_csv(os.path.join(dataPath,file), index_col = 0)

    colours = ['r' if (abs(data.loc[x, 'logFC'])>1 and data.loc[x, sigType]<0.05)  else 'c' if data.loc[x, sigType]<0.05 else 'k' for x in data.index ]
    # plt.scatter(data['logFC'],np.log(data[sigType])*-1, c=colours)
    
    if plot:
        plt.scatter(data['logFC'],np.log(data[sigType])*-1, c=data['logFC'], cmap='coolwarm', vmin = -max(abs(data['logFC'])), vmax = max(abs(data['logFC'])))
        plt.ylabel('-log10 ' + sigType)

        plt.xlabel('log2 Fold Change')
        plt.axhline(-np.log(0.05))
        plt.axvline(np.log2(2), c='r', dashes=[5,3])
        plt.axvline(-np.log2(2), c='r', dashes=[5,3])
        plt.title(file.split('.')[0][file.rfind('/')+1:]) # Use path and extension trimmed file name as figure title

    for gene in data.index:
        if (data.loc[gene,sigType] < 0.05):
            sigGenes.append(gene)
            label = gene
            if plot:
                plt.annotate(label, # this is the text
                             (data.loc[gene,'logFC'],-np.log(data.loc[gene,sigType])), # these are the coordinates to position the label
                             textcoords="offset points", # how to position the text
                             xytext=(0,10), # distance from text to points (x,y)
                             ha='left') # horizontal alignment can be left, right or center
                
    if plot:
        outfile = file.split('.')[0] + '.png'   
        print(outfile)
        plt.savefig(os.path.join(dataPath,outfile), dpi='figure', format='png')
        outfile = file.split('.')[0] + '.svg'   
        plt.savefig(os.path.join(dataPath,outfile), dpi='figure', format='svg')
        plt.show()
    return(sigGenes)
    