from openpyxl import load_workbook
import numpy as np
import pandas as pd
from copy import copy# as copy


class master_data:
    def __init__(self, dataPath):
        ### import data from excel workbook
        
        self.wb = load_workbook(dataPath)
        self.ws = self.wb['Exported dataset']
        self.values = [[y.value for y in x] for x in self.ws[self.ws.calculate_dimension()]]

        self.dropData = False
        self.threshold = False
        
        ### Convert nested list to a pandas dataFrame and extract expression data with labels
    def get_data(self, fix_zeros=True):
        df = pd.DataFrame(self.values)
        col3 = df.iloc[:,3].tolist()
        self.targIdx = col3.index('Target name (display name)') + 1

        
        rowLabels = df.iloc[self.targIdx:,3]
        # rowLabels = [x.split(' (')[0] for x in rowLabels.values]
        rowLabels = dict(zip([x for x in range(self.targIdx,self.targIdx+len(rowLabels))],rowLabels))
        rowLabels

        colLabels = df.iloc[0,4:]
        colLabels = [x.replace(' | ','_') for x in colLabels.values]
        colLabels = dict(zip([x for x in range(4,4+len(colLabels))],colLabels))
        colLabels

        self.data = df.iloc[self.targIdx:,4:].astype(np.float32)
        self.data.rename(index=rowLabels,columns=colLabels, inplace=True)

        if fix_zeros:
            for x in self.data.columns:
                # print(x)
                for y in self.data.index:
                    # print(y)
                    # Subtract 1 from hyb pos and hyb neg values as we currently don't use RCC files to determine which values have been changed from 0 to 1 for Neg control.
                    # Canging hyb-pos values may help to ameliorate effects of changing true 1 values to zero, average decrease in hyb neg values.
                    if (y == 'HYB-NEG'):
                        self.data.loc[y,x] = self.data.loc[y,x] -1
                    # elif (y == 'HYB-POS'):
                    #     self.data.loc[y,x] = self.data.loc[y,x] -1
                    else:
                        if (self.data.loc[y,x] == 1):
                            print(x,y)
                            self.data.loc[y,x] = 0.0

        self.sampleInfo = pd.DataFrame(df.iloc[0:self.targIdx-1,4:])
        self.sampleInfo.rename(index=df.iloc[0:self.targIdx-1,0], columns=colLabels, inplace=True)
        print('sampleInfo.shape')
        print(self.sampleInfo.shape)
        
        print('data.shape')
        print(self.data.shape)

        self.dataOrig = self.data.copy()
        # Log transform data for QC and analysis steps
        self.dataLog1 = np.log2(self.data+1)
        
        self.probeClass = df.iloc[self.targIdx:,2]      ### Index needs updating here also
        self.probeClass.rename(index=rowLabels, inplace=True)
        self.probeClass.rename(index='ProbeClass', inplace=True)
        self.probeClassDict = {
            'Positive': 'A',
            'Negative': 'B',
            'Control': 'C',
            'Endogenous': 'E'
        }
        return self.data.copy(), self.dataLog1.copy(), self.sampleInfo.copy()
    
    def get_descriptors(self):
        ### Extract descriptions for each sample
        nuclei = sampleInfo.loc['AOI nuclei count']
        surfArea = sampleInfo.loc['AOI surface area']

        print(sampleInfo.shape)
        
    def add_class_mean(self, df):
        ## Add column to data with mean values for each probe (row)
        mean = df.mean(axis = 1)
        df = df.assign(mean=mean.values)

        ## Add column to data with probe class for each probe
        df = df.assign(probeClass=[self.probeClassDict[v] for v in self.probeClass.values])

        ### Extract lists of controls and their values
        self.posCTLs = self.probeClass.index[self.probeClass== 'Positive'].tolist()
        self.negCTLs = self.probeClass.index[self.probeClass== 'Negative'].tolist()
        self.IgCTLs = copy(self.negCTLs)
        try:
            self.IgCTLs.remove('HYB-NEG')
        except ValueError:
            pass
        self.HK = self.probeClass.index[self.probeClass== 'Control'].tolist()
        self.endog = self.probeClass.index[self.probeClass== 'Endogenous'].tolist()

        print('Positive Control count:\t{:d}, {}'.format(len(self.posCTLs), self.posCTLs))
        print('Nagative Control count:\t{:d}, {}'.format(len(self.negCTLs), self.negCTLs))
        print('Ig Control count:\t{:d}, {}'.format(len(self.IgCTLs), self.IgCTLs))
        print('HK Control count:\t{:d}, {}'.format(len(self.HK), self.HK))
        print('Endogenous probe count:\t{:d}, {}'.format(len(self.endog), self.endog))

        return df.copy(), self.sampleInfo.copy()

    def drop_AOIs(self, includes, writeOrig=False):

        dropAOIs = [x for x in list(self.data.columns) if (x in includes)]
        print(dropAOIs)
        if writeOrig:
            self.dataOrig = self.dataOrig.drop(labels=dropAOIs, axis=1, inplace=True)

        # print(self.dataLog1.shape)
        self.dataLog1.drop(labels=dropAOIs, axis=1, inplace=True)
        # print(self.dataLog1.shape)

        try:
            self.ERCCData.drop(labels=dropAOIs, axis=1, inplace=True)
        except:
            pass

        self.sampleInfo.drop(labels=dropAOIs, axis=1, inplace=True)
        return self.dataLog1.copy(), self.sampleInfo.copy()
        
    def set_threshold(self, threshold):
        self.threshold = threshold
        
        
        # ToDo: Check that all values in master data are also included in threshold dataFrame
        # ToDo: Convert threshold data to 0/1 data if needed
        
    def drop_probes(self, labels):
        try:
            assert type(labels) == list
        except:
            print('labels need to be a list')
            return False
            
        if not self.dropData:
            self.dropData = self.dataLog1

        self.dropData.drop(labels=labels)
        
        return self.dropData.copy()
        
    def ERCC_norm(self):
        self.ERCCData = self.dataLog1
            
        try:
            self.ERCCData = self.ERCCData.drop(labels=['mean'], axis=1)
        except:
            pass
        try:
            self.ERCCData = self.ERCCData.drop(labels=['probeClass'], axis=1)
        except:
            pass
            
        # ERCC normalisation. 
        # Divide by individual HYB-POS values, then scale data using the geometric mean of all HYB-POS values. subtract in log space is same as divide in normal space. Add in log space is same as multiply in normal space
        self.ERCCData = self.ERCCData - self.ERCCData.loc['HYB-POS'] + np.mean(self.ERCCData.loc['HYB-POS'])

        # #ToDo: set below threshold values to 0
        # if (self.threshold):
        #     self.ERCCData = self.ERCCData * self.threshold

        #set any negative values following ERCC noirmalisation to 0. These are at or below the limit of detection
        for x in self.ERCCData.columns:
            # print(x)
            for y in self.ERCCData.index:
                if (self.ERCCData.loc[y,x] < 0.0):
                    # print(x,y)
                    self.ERCCData.loc[y,x] = 0.0

# #         zeroFilter = self.ERCCData > 0
#         self.ERCCData = self.ERCCData * (self.ERCCData>0)
        
        return self.ERCCData.copy()





def read_Surf_Areas(wsPath, indexList, columnList):
    with open(wsPath, 'r')as f:
        results = []
        record = False
        lines = f.readlines()
        for line in lines:
            line=line.strip()
            if line.startswith('Area '):
                record = True
            elif line.startswith('Totals '):
                record = False
            elif record:
                fields = line.split()
                results.append([int(x) for x in fields[1:]])
    return(pd.DataFrame(results, index=indexList, columns = columnList))        
