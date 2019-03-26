import loompy
import glob
import velocyto as vcy
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys

def get_loom(x):
    if x == '180504':
        #files = glob.glob('/data/sc-10x/data-runs/171120-scheele-adipose/*/velocyto/*.loom')
        #files = list(filter(lambda x: '180831' not in x, files))
        #files = sorted(files, key=lambda x: x.split('/')[-1])
        #loompy.combine(files, '../output/velocyto/10x-180504.loom', key="Accession")
        vlm = vcy.VelocytoLoom('../output/velocyto/10x-180504.loom')
    else:
        #files = sorted(glob.glob('/data/sc-10x/data-runs/171120-scheele-adipose/180831_10x_s*/velocyto/*.loom'))
        #loompy.combine(files, '../output/velocyto/10x-180831.loom', key="Accession")
        vlm = vcy.VelocytoLoom('../output/velocyto/10x-180831.loom')
    return vlm

def rename_cellid(x):
    if ('BAT8-3000_cells' in x):
        return x.replace('BAT8-3000_cells', 'Supra_4')
    elif ('44BSAT-3000_cells' in x):
        return x.replace('44BSAT-3000_cells', 'Subq_4')
    else:
        return x

def rename_colnames(x, vlm):
    if x == '180504':
        #rename BAT8 and 44BSAT samples
        vlm.ca['CellID'] = np.asarray(list(map(rename_cellid, vlm.ca['CellID'])))
        vlm.ca['sample_name'] = np.asarray(list(map(lambda x: x.split(':')[0], vlm.ca['CellID'])))
        print(set(vlm.ca['sample_name']))
    else:
        vlm.ca['sample_name'] = np.asarray(list(map(lambda x: x.split(':')[0], vlm.ca['CellID'])))
        vlm.ca['CellID'] = np.asarray(list(map(lambda x: x.split(':')[1][:-1] + '-' + x[12], vlm.ca['CellID'])))
    return vlm

def merge_vlm_with_metadata(vlm, metadata, x):
    
    if x == '180504':
        metadata['CellID-old'] = metadata.index
        metadata['CellID'] = metadata[['CellID-old', 'sample_name']].apply(lambda x: x[1] + ':' + x[0].split('-')[0] + 'x', axis=1)
    else:
        metadata['CellID'] = metadata.index

    cells_to_keep = np.array(metadata['CellID'].tolist())
    cellids_vlm = np.array(vlm.ca['CellID'])

    vlm.filter_cells(bool_array=np.isin(cellids_vlm, cells_to_keep))

    ca_df = pd.DataFrame(vlm.ca)
    if x == '180504':
        merged = ca_df.merge(metadata, on='CellID', how='right')
    else:
        merged = ca_df.merge(metadata, on='CellID', how='left')
    merged = merged.to_dict('list')


    vlm.ca['CellID'] = np.asarray([x.encode('utf8') for x in vlm.ca['CellID']])
    vlm.ca['sample_name'] = np.asarray([x.encode('utf8') for x in vlm.ca['sample_name']])
    
    for key in merged:
        try:
            vlm.ca[key] = np.asarray([x.encode('utf8') for x in merged[key]])
        except:
            vlm.ca[key] = np.asarray(merged[key])
            
    return vlm


def main():
    if len(sys.argv) != 2:
        print('usage: python velocyto_preprocessing.py [180504 OR 180831]')
        sys.exit()
    elif sys.argv[1] != '180504' and sys.argv[1] != '180831':
        print('usage: python velocyto_preprocessing.py [180504 OR 180831]')
        sys.exit()
    
    dataset = sys.argv[1]
    print('aggregating loom files')
    vlm = get_loom(dataset)
    print('rename colnames')
    vlm = rename_colnames(dataset, vlm)
    
    print('read metadata')
    if dataset == '180504':
        metadata = pd.read_table('../tables/10x-180504-downsampled36-metadata.txt', sep='\t')
    else:
        metadata = pd.read_table('../tables/10x-180831-metadata.txt', sep='\t')
    
    print('merge vlm with metadata')
    vlm = merge_vlm_with_metadata(vlm, metadata, dataset)
    
    print('save loom as hdf5')        
    if dataset == '180504':
        vlm.to_hdf5('../output/velocyto/10x-180504-downsampled36.hdf5')
    else:
        vlm.to_hdf5('../output/velocyto/10x-180831.hdf5')

if __name__ == '__main__':
    main()
