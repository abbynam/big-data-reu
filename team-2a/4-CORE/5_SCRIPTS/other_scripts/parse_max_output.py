import sys, os
import collections
import pandas
import math

MAX_VALS_FILE = '/scratch1/max_values.txt'
OUTPUT_FILE = '/y_drive/temp/kde_results_20210221.xls'

def get_params_dict(rec):
    
    d = collections.OrderedDict()
   
    pairs = rec.replace('.0__', "__").split('dose-')[-1].split(".")[0].split("__")
    for p in pairs:
        parts = p.split('-')
        print(parts)
        if len(parts) == 3:
            parts[1] = f'{parts[1]}-{parts[2]}'
        d[parts[0]] = parts[1]
    
        
    return d

def get_results_df(txt):
    dictlist = []
    for rec in txt:
        x, y, z = [float(x) for x in rec.split(":")[-1].strip().split(", ")]
        d = get_params_dict(rec)
        d['x'] = x
        d['y'] = y
        d['z'] = z
        d['x0'] = float(d['EVENT_FILE_PATH'].replace('x',""))
        d['y0'] = -6.0
        d['z0'] = 13.0
        d['dx'] = d['x'] - d['x0']
        d['dy'] = d['y'] - d['y0']
        d['dz'] = d['z'] - d['z0']
        d['r'] = math.sqrt(d['dx']*d['dx'] + d['dy']*d['dy'] + d['dz']*d['dz'])
        
        dictlist.append(d)
        
    print(d.keys())
    
    df = pandas.DataFrame(dictlist, columns=d.keys())
    
        
    return df
    
def main():
    txt = open(MAX_VALS_FILE).readlines()
    df = get_results_df(txt)
    df.sort_values(by=['x0', 'r'], ascending=[True, True], inplace=True)
    print(df.head(100))
    df.to_excel(OUTPUT_FILE, index=False, header=True)


    
        

main()