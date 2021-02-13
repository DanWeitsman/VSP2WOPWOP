"""
PSU-WOPWOP Multi Observer Post Processing Script
1/21/21
"""
import os
import wopwop

case_directory = os.getcwd()
f1 = lambda a: wopwop.extract_wopwop_quant(case_directory=a,prefix = 'pressure')
wopwop.apply_to_namelist([f1], cases_directory=case_directory, cases='cases.nam')
