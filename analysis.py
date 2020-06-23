import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import coordinate_conversions as cc
import plotting as plot

data_files_path = '/Users/cam/Desktop/astro_research/orion/snn/'
filenames = [data_files_path + 'final_result_20191105.csv',data_files_path + 'apogee_query.csv']
df_combined = cc.combine_chen_cam(filenames[0],filenames[1],False)

lambda_ori_labels = [4,26,27]

# Reveals that average velocity of lambda orion is (dx,dy,dz) = (-13.79,0.78,1.48) km/s. Use this for b30 and b35

for i in range(1,27): 
    lori = cc.grouping_ICRS(df_combined,[i])
    lori_icrs = cc.rf_equitorial(lori) # in the kinematic reference frame
    lori, lori_center = cc.rf_cartesian(lori_icrs) # in the kinematic reference fram
    #plot.cartesian_3d([lori, lori_center], colors=['purple','black'], labels = ['Lambda Orion stars', 'Lambda Orion center' ], linewidths = [.8,1.5])
    cc.implement_ward_2018(lori)


#lori26 = cc.grouping_ICRS(df_combined, [26])
#lori26_icrs = cc.rf_equitorial(lori26) # in the kinematic reference frame
#lori26 = cc.rf_cartesian(lori26_icrs) # in the kinematic reference fram
#
#lori27 = cc.grouping_ICRS(df_combined, [27])
#lori27_icrs = cc.rf_equitorial(lori27) # in the kinematic reference frame
#lori27 = cc.rf_cartesian(lori27_icrs) # in the kinematic reference fram

#df_combined = df_combined.loc[df_combined.probability > 0.2]
#ngc1980 = cc.grouping_ICRS(df_combined.loc[df_combined.probability > 0.2],[5])
#ngc1980_icrs = cc.rf_equitorial(ngc1980) # in the kinematic reference frame
#ngc1980, ngc1980_center = cc.rf_cartesian(ngc1980_icrs) # in the kinematic reference fram
#cc.implement_ward_2018(ngc1980_icrs)

#orionY = cc.grouping_ICRS(df_combined,[2])
#orionY_icrs = cc.rf_equitorial(orionY) # in the kinematic reference frame
#orionY, orionY_center = cc.rf_cartesian(orionY_icrs) # in the kinematic reference fram
#
#plot.cartesian_3d([ngc1980, ngc1980_center, orionY, orionY_center], colors=['blue','black','green','black'], labels=['NGC1980','center', 'Orion Y', 'center'], linewidths = [.5,.8,.5,.8])



#plot.cartesian_3d([lori], colors=['blue'], labels=['L-ori'], linewidths = [.8])


#plot.cartesian_3d([lori,lori26,lori27], colors=['blue', 'green', 'black'], linewidths = [.6,.6,.6], labels = ['Lambda Orion Main', 'LOri26', 'LOri27'])


average_velocities = (-13.78727137270535, 0.7818539838416585, 1.3838735000868423) 


#plot.equitorial_2d(lori_icrs)

#plot.cartesian_3d([lori,lori_b30,lori_b35], colors=['grey','red','blue'],labels=['Lambda Orion', 'B30 Region', 'B35 Region'], linewidths=[.4,.8,.8])



#------ Plotting specific regions of lori ------#

# Stars in region of b30 gas cloud
#ra_region = [82,83.5]
#dec_region = [12,14]
#region_select = (ra_region,dec_region)
#lori_b30_icrs = cc.grouping_ICRS(df_combined,lambda_ori_labels,region_select)
#lori_b30 = cc.rf_cartesian(lori_b30_icrs, average_velocities)
#
## Stars in region of b35 gas cloud
#ra_region = [85,87]
#dec_region = [8,10]
#region_select = (ra_region,dec_region)
#lori_b35_icrs = cc.grouping_ICRS(df_combined,lambda_ori_labels,region_select)
#lori_b35 = cc.rf_cartesian(lori_b35_icrs, average_velocities)
