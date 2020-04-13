import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import coordinate_conversions as cc
import plotting as plot

data_files_path = '/Users/admin/Desktop/astro_research/orion/snn/'
filenames = [data_files_path + 'final_result_20191105.csv',data_files_path + 'apogee_query.csv']
df_combined = cc.combine_chen_cam(filenames[0],filenames[1],False)

lambda_ori_labels = [4,26,27]

# Reveals that average velocity of lambda orion is (dx,dy,dz) = (-13.79,0.78,1.48) km/s. Use this for b30 and b35
lori_icrs = cc.grouping_ICRS(df_combined,lambda_ori_labels)
lori = cc.rf_cartesian(lori_icrs) 

average_velocities = (-13.78727137270535, 0.7818539838416585, 1.3838735000868423) 

# Stars in region of b30 gas cloud
ra_region = [82,83.5]
dec_region = [12,14]
region_select = (ra_region,dec_region)
lori_b30_icrs = cc.grouping_ICRS(df_combined,lambda_ori_labels,region_select)
lori_b30 = cc.rf_cartesian(lori_b30_icrs, average_velocities)

# Stars in region of b35 gas cloud
ra_region = [85,87]
dec_region = [8,10]
region_select = (ra_region,dec_region)
lori_b35_icrs = cc.grouping_ICRS(df_combined,lambda_ori_labels,region_select)
lori_b35 = cc.rf_cartesian(lori_b35_icrs, average_velocities)

plot.cartesian_3d([lori,lori_b30,lori_b35], colors=['grey','red','blue'],labels=['Lambda Orion', 'B30 Region', 'B35 Region'], linewidths=[.4,.8,.8])

