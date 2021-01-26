"""
PSU-WOPWOP Multi Observer Post Processing Script
1/21/21
"""

import wopwop

case_directory = '/Users/danielweitsman/Desktop/Masters_Research/UAM_param_study/Runs/linear_twist/'
f1 = lambda a: wopwop.extract_wopwop_quant(case_directory=a,prefix = 'pressure')
wopwop.apply_to_namelist([f1], cases_directory=case_directory, cases='cases.nam')
