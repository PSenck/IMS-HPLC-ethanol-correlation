import numpy as np
import pandas as pd
from read_mea import read_mea
import os.path
import glob
from scipy.interpolate import interp1d
from BaselineRemoval import BaselineRemoval
from copy import deepcopy

def load_raw_data(path_list, S_ds = 1,  N_ds = 1, start_dict = None , end_dict = None):

    """ load IMS data, it is possible to filter according to Timestamps in filename

    INPUT:
    path_list ... list with folder paths for respective measurement days
    S_ds ... int, optional, factor to downsample spectra according to measurements - every S_ds measurement, default = 1  
    N_ds ... int, optional, factor to downsample spectra according to drift and retention_times - every N_ds drift/ret_time, default = 1
    start_dict ... dictionary with start points for the respective experiment, default: None
    end_dict ... dictionary with end points for the respective experiment, default: None
    Note: elemnts in path_list and in start/end_dict have to match if not None
    If no start or end_dict are passed, whole files will be loaded

    Output:
    dat_wide_as_dict ... dictionary with key : timestamp, value : pandas DataFrame 
    """

    files_list = []
    for path in path_list:
        files = glob.glob(os.path.join(path, '*.mea'))  
        files.sort() 
        files_list.append(files)

    files_df = []
    for files in files_list:
        g = pd.DataFrame(files)             
        t = pd.to_datetime(files, format='%y%m%d_%H%M%S', exact=False )
        g.set_index(t, inplace= True)
        files_df.append(g)     

    file_filter_list = []

    if start_dict == None and end_dict != None:

        for df, [ferm, end] in zip(files_df, end_dict.items()):
            df = df[(df.index <= end)]
            file_filter_list.append(df)

    if end_dict == None and start_dict != None:

        for df, [ferm, start] in zip(files_df, start_dict.items()):
            df = df[(df.index >= start)]
            file_filter_list.append(df)

    if start_dict != None and end_dict != None:
        
        for df, start, [ferm, end] in zip(files_df, start_dict.values(), end_dict.items()):
            df = df[(df.index >= start ) &  (df.index <= end)]
            file_filter_list.append(df)
    
    if start_dict == None and end_dict == None:
        file_filter_list = files_df


    merged_files = pd.concat(file_filter_list)

    dat_wide_as_dict = {}
    for f in merged_files.iloc[::S_ds,0]:
        print('Reading file ' + f)

        # read binary files
        values, meta_attr, ret_time, drift_time = read_mea(f)

        # extract timestamp from metadata
        t = pd.to_datetime(meta_attr['Timestamp'])

        # turn data into pandas dataframe and store in dict
        tmp = pd.DataFrame(values[::N_ds, ::N_ds], columns=drift_time[::N_ds])    
        tmp['ret_time'] = ret_time[::N_ds]   

        dat_wide_as_dict[t] = tmp.set_index('ret_time')
    
    return dat_wide_as_dict



def ret_time_filter(data_dict, ret_times):
    """filtering 2D spectra for desired retention times

    INPUT:
    data_dict ... raw or preprocessed data_dict with key: timestamp, value: 2D-spectrum
    ret_times ... list with desired retention times

    OUTPUT:
    ret_time_filtered_dict ... filtered data_dict with key: timestamp, value: 2D-spectrum
    """


    ret_time_filtered_dict = {}
    for t, df in data_dict.items():
        ret_time_index = [data_dict[t].index.get_loc(r, 'nearest') for r in ret_times]
        ret_time_filtered_dict[t] = df.iloc[ret_time_index]

    return ret_time_filtered_dict



def integrate_spectrum(data_dict):
    """integrate each spectrum in data_dict at each drift_time over all retention times

    INPUT:
    data_dict ... dictionary with key: timestamp, value: pd.DataFrame
    
    OUTPUT:
    col_sum_df : pd.DataFrame with integrated values, each spectrum corresponds to one row. With index: timestamps 
    """

    col_sum_list = []
    for df in data_dict.values():
        col_sum = np.trapz(df, x = df.index, axis = 0)
        col_sum_list.append(col_sum)

    col_sum_df = pd.DataFrame(col_sum_list, columns = df.columns, index = data_dict.keys())

    return col_sum_df


def baseline_correction(dataframe, algorythm = "Zhang", polynomial_degree=2):
    """Perform baseline with three possible algorithms 

    INPUT:
    dataframe ...  pd.DataFrame with integrated values, each spectrum corresponds to one row
    type ... Zhang, Modpoly, Imodpoly choose one of the different algorythms
    polynomial_degree ... the degree the a polynomial, only necessary for Modpoly and Imodpoly

    OUTPUT: baseline_corrected_df
    """

    tmp = []
    for t,row in dataframe.iterrows():
        base_obj=BaselineRemoval(row)
        if algorythm == "Zhang":
            output = base_obj.ZhangFit()

        elif algorythm == "Modpoly":
            output = base_obj.ModPoly(polynomial_degree)

        elif algorythm == "Imodpoly":
            output = base_obj.IModPoly(polynomial_degree)
        else:
            raise ValueError("Choose Zhang, Modpoly or Imodpoly")

        tmp.append(output)

    baseline_corrected_df = pd.DataFrame(tmp, columns = dataframe.columns, index = dataframe.index)

    return baseline_corrected_df



def get_value(array, start, end, desired_val):
    """searching between given borders for the nearest value to desired_val

    INPUT:
    array ... array of values which could contain the nearest value
    start ... start_border
    end ... end_border
        The value of interrest can only be between this borders, the two ethanol peak should lie within the borders 
    desired_val ... value for which the next value in the array between the borders is searched

    OUTPUT:
    nearest_val: if present, nearest value to desired_val, if not return NaN

    Warning: This function is probably prone to undesired outcomes with complex chemical mixtures.
    """
    filtered_array = [i for i in array if start <= i <= end]
    
    if not filtered_array:
        nearest_val = np.nan
    else:
        
        abs_diff_func = lambda val : abs(val - desired_val)
        nearest_val = min(filtered_array, key=abs_diff_func)
    
    return nearest_val



def integrate_peaks(pd_series, left, right):
    """calculate for each spectrum and each peak the peak_area
    INPUT:
    pd_series ... pd.series which contain the spectral values
    left ... left_ips value of peak
    right ... right_ips value of peak

    OUTPUT:
    area ... peak area, if no left or right values are given, return 0.0
    """  

    if np.isnan(left) or np.isnan(right):
        area = 0.0
        

    else:
        filter = (pd_series.index >= left) & (pd_series.index <= right)
        filtered = pd_series[filter]
        area = np.trapz(filtered.values, x = filtered.index.values)
        
    return area
