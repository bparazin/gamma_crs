#Import statements
import pandas as pd
import numpy as np

#Define pb210 decay constant
LAMBDA = 0.03114

#define all of the functions used here
def import_data(path):
    raw_data = pd.read_csv(path)
    data_struct = pd.DataFrame()
    data_struct['depth_low'] = raw_data['depth_low (cm)']
    data_struct['volume'] = raw_data['volume (cm3)']
    data_struct['weight'] = raw_data['weight (g)']
    data_struct['210pb'] = raw_data['210Pb']
    data_struct['210pb_err'] = raw_data['210Pb_err']
    return data_struct

def get_crs_age(depth_low, volume, weight, pb, equilbrium_row, pb_err, pb_bkg = None, p = 0.07, A_0_method = 'default', ref_line = None, ref_year = None):
    
    if A_0_method == 'refdate':
        assert ref_line is not None
        assert ref_year is not None
    
    if pb_bkg == None:
        pb_bkg = np.mean(pb[equilbrium_row:])
    density = (weight / volume)
    c = (pb - pb_bkg)
    m = []
    for i, density_i in enumerate(density):
        if i == 0:
            m_i = density_i * depth_low[i]
        else:
            m_i = m[i-1] + density_i * (depth_low[i] - depth_low[i-1])
        m.append(m_i)
    m = np.asarray(m)
    
    a_hat = []
    a_hat_err = []
    for i, c_i in enumerate(c):
        if i == 0:
            a_hat_i = c_i * m[i] / 2
            a_hat_i_err = np.sqrt((pb_err[i]**2 + c_i**2 * p**2) * m[i]**2)
        else:
            a_hat_i = a_hat[i-1] + (c_i + c[i-1]) * (m[i] - m[i-1]) / 2
            a_hat_i_err = np.sqrt((pb_err[i]**2 + c_i**2 * p**2) * (m[i] - m[i-1])**2)
        a_hat.append(a_hat_i)
        a_hat_err.append(a_hat_i_err)
    a_hat = np.asarray(a_hat)
    a_hat_err = np.asarray(a_hat_err)

    a_0, a_0_err = get_A_0(a_hat, a_hat_err, c, m, pb, pb_err, equilbrium_row, A_0_method, ref_line, ref_year)
    
    t = 1 / LAMBDA * np.log(a_0 / (a_0 - a_hat))
    t_err = 1 / LAMBDA * np.sqrt(np.abs((a_0_err/a_0)**2 + (1 - 2*(a_0 - a_hat)/a_0) * (a_0_err**2 + a_hat_err**2)/(a_0-a_hat)**2))
    
    return t[:equilbrium_row], t_err[:equilbrium_row]

def get_A_0(a_hat, a_hat_err, c, m, pb, pb_err, equilbrium_row, A_0_method, ref_line, ref_year):
    if A_0_method == 'default':
        A_0 = (a_hat[equilbrium_row-1] + c[equilbrium_row-1] * (m[equilbrium_row] - m[equilbrium_row-1]) / 2) * np.ones(np.shape(a_hat))
        A_0_err = np.sqrt(pb_err[equilbrium_row-1]**2 + a_hat_err[equilbrium_row-1]**2) * np.ones(np.shape(a_hat))
        
    if A_0_method == 'refdate':
        A_0 = ((a_hat[ref_line] + a_hat[ref_line]/(np.exp(LAMBDA * ref_year) - 1)) * np.heaviside(ref_line - np.arange(0, len(a_hat), 1), 1)
              + (a_hat[equilbrium_row-1] + c[equilbrium_row-1] * (m[equilbrium_row] - m[equilbrium_row-1]) / 2) * np.heaviside(np.arange(0, len(a_hat), 1) - ref_line, 0))
        A_0_err = (np.sqrt(pb_err[ref_line]**2 + a_hat_err[ref_line]**2) * np.heaviside(ref_line - np.arange(0, len(a_hat), 1), 1)
                   + (np.sqrt(pb_err[equilbrium_row-1]**2 + a_hat_err[equilbrium_row-1]**2) * np.ones(np.shape(a_hat)) * np.heaviside(np.arange(0, len(a_hat), 1) - ref_line, 0)))
    
    return A_0, A_0_err

def plot_years(depths, years, yr_err):
    fig, ax = plt.subplots(figsize = (10,10))
    ax.errorbar(years, depths[:len(years)], xerr = yr_err, fmt = 'o')
    ax.invert_yaxis()
    ax.set_ylabel('Depth (cm)')
    ax.set_xlabel('Year')
    ax.grid(alpha = 0.5)
    plt.show
    
def write_to_csv(depth, age, error, name):
    assert len(age)==len(error)
    if len(age) != len(depth):
        nan_arr = np.empty(len(depth)-len(age))
        nan_arr[:] = np.nan
        age = np.concatenate((age, nan_arr))
        error = np.concatenate((error, nan_arr))
    df = pd.DataFrame()
    df['lower depth (cm)'] = depth
    df['age (years)'] = age
    df['age error (years)'] = error
    df.to_csv(name)

    
#interface with the user to get the needed info
path = input('Hello, please input the name of the csv file set up as described in the readme, including the .csv extension: ')
data_struct = import_data(path)
refdateinput = input('y/n do you want to use a reference date?: ')
while refdateinput not in ['y', 'n']:
    refdateinput = input('please enter y or n, do you want to use a reference date?: ')
refdate = 'refdate' if  refdateinput == 'y' else 'default'
ref_line = None
ref_year = None
if refdate == 'refdate':
    ref_line = int(input('please enter the row number of the reference date: ' )) - 2
    ref_year = float(input('please enter (in years) what the age of this row was determined to be: ' ))
    
eq_row = int(input('please enter the row number of the first row where the pb210 values have reached equilibrium: ')) - 2
#both of these are minus 2 since excel is 2 indexed, since the first row of each column is the header and it starts counting at one
p_in = input('please enter the constant percentage error in dry mass increments, or hit enter to use the default (7%): ')
if p_in == '':
    p = 0.07
else:
    p = float(p_in)
    
pb_bkg_in = input('please enter the pb210 background value or hit enter to use the default calculation method: ')
if pb_bkg_in == '':
    pb_bkg = None
else:
    pb_bkg = float(pb_bkg)

time, time_error = get_crs_age(data_struct['depth_low'], data_struct['volume'], data_struct['weight'], data_struct['210pb'], eq_row, data_struct['210pb_err'],
                               pb_bkg, p, refdate, ref_line, ref_year)

name = input('Please enter the name of the file you want to output the results to, including the .csv extension: ')
write_to_csv(data_struct['depth_low'], time, time_error, name)
print('Successfully completed')