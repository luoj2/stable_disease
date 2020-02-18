
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import publication_settings
import pandas as pd

from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection

matplotlib.rcParams.update(publication_settings.params)

labels = ['PD','SD responder','minor PR']

def array_fraction(array,condition):
    return np.sum(array==condition)/len(array)

# extract SD data
data = pd.read_csv('SD2_SDfeatures.csv')
#SD_feat = np.array(data['Smoking'])
SD_feat = np.array(data['STK11'])

BOR_SD = np.array(data['BOR'])
PFS_SD = np.array(data['PFS'])

# define SD -0 to -30 percent AND PFS>6mo
mask1 = (BOR_SD <=0) & (PFS_SD>6)
SD_feat = SD_feat[mask1]

known = (SD_feat != 'Unknown')
SD_feat = SD_feat[known]#.astype(np.float64).flatten()
#SD_featfloat = SD_feat[known].astype(np.float64)

never = (SD_feat == 'no_alt')
#never = (SD_feat == 0)
never = SD_feat[never]
SD_ever = (len(SD_feat)-len(never))/len(SD_feat)
SD_never = len(never)/len(SD_feat)

# extract PR data
data = pd.read_csv('SD2_PRfeatures.csv')
#PR_feat = np.array(data['Smoking'])
PR_feat = np.array(data['STK11'])

BOR_PR = np.array(data['BOR'])
# define PR minor
mask = BOR_PR >-57
PR_feat = PR_feat[mask]
known = (PR_feat != 'Unknown')
PR_feat = PR_feat[known]#.astype(np.float64).flatten()

never2 = (PR_feat == 'no_alt')
#never2 = (PR_feat == 0)
never2 = PR_feat[never2]
PR_ever = (len(PR_feat)-len(never2))/len(PR_feat)
PR_never = len(never2)/len(PR_feat)

# extract PD data
data = pd.read_csv('PDfeatures.csv')
#PD_feat = np.array(data['Smoking'])
PD_feat = np.array(data['STK11'])

known = (PD_feat != 'Unknown')
PD_feat = PD_feat[known]#.astype(np.float64).flatten()

never3 = (PD_feat == 'no_alt')
#never3 = (PD_feat == 0)
never3 = PD_feat[never3]
PD_ever = (len(PD_feat)-len(never3))/len(PD_feat)
PD_never = len(never3)/len(PD_feat)


# 5 x 5 inches
total_width = 5
total_height = 2.5
fig = plt.figure(1, figsize=(total_width,total_height))

# give all in inches
top_border = 0.2
bottom_border = 0.5
left_border = 1.4
right_border = 0.4

left = left_border/total_width
bottom = bottom_border/total_height
width = (total_width - left_border - right_border)/total_width
height = (total_height - top_border - bottom_border)/total_height

ax = fig.add_axes([left,bottom,width,height])

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()

ax.set_xlim([0,1.375])
ax.set_ylim([0,len(labels)])

ax.tick_params('y',direction='out',which='both')
#ax.set_yticks(minor_ticks,minor=True)

font = 'Arial'

ax.set_xticks([0,0.25,.5,.75,1.])
ax.set_xticklabels((r"0%",r"25%",r"50%",r"75%",r"100%"),fontsize=12)
y=np.array([0.1,0.5,1.1,1.5,2.1,2.5])
ax.set_yticks(y)
ax.vlines(0,0.5,2.5)
#ax.xaxis.set_tick_params(length=0)
tick_labels = ('n=%i' %len(PD_feat), labels[0], 'n=%i' %len(SD_feat), labels[1], 'n=%i' %len(PR_feat), labels[2])
ax.set_yticklabels(tick_labels, fontsize=12)
plt.setp(ax.get_yticklines()[0],visible=False)
plt.setp(ax.get_yticklines()[4],visible=False)
plt.setp(ax.get_yticklines()[8],visible=False)


#ax.set_xlabel('mut/mb',fontname=font,fontsize=12)

#SD
#blue3 = '#eff3ff'
#blue3 = '#c6dbef' #2
#blue3 = '#9ecae1' #real
#blue3 = '#6baed6'
#blue3 = '#3182bd' #5
blue3 = '#08519c'

#PR
#orange3 = '#fef0d9'
#orange3 = '#fdd49e' #2
#orange3 = '#fdbb84' #real
#orange3 = '#fc8d59'
#orange3 = '#e34a33' #5
orange3 = '#b30000'

#PD
#purp3 = '#f2f0f7'
#purp3 = '#dadaeb' #2
#purp3 = '#bcbddc' #real
#purp3 = '#9e9ac8'
#purp3 = '#756bb1' #5
purp3 = '#54278f'

grey = '#919191' #med
#grey3 = '#f7f7f7' #light
grey3 = '#252525' #dark

print(PD_ever,PD_never,PR_ever,PR_never,SD_ever,SD_never)

ax.add_patch( Rectangle([0,0.3],PD_ever,0.4,edgecolor='k',facecolor=orange3) )
ax.add_patch( Rectangle([PD_ever,0.3],PD_never,0.4,edgecolor='k',facecolor=grey3) )
ax.vlines(PD_ever,0.3,0.7,linewidth=3,edgecolor='k')

ax.add_patch( Rectangle([0,1.3],SD_ever,0.4,edgecolor='k',facecolor=blue3) )
ax.add_patch( Rectangle([SD_ever,1.3],SD_never,0.4,edgecolor='k',facecolor=grey3) )
ax.vlines(SD_ever,1.3,1.7,linewidth=3,edgecolor='k')

ax.add_patch( Rectangle([0,2.3],PR_ever,0.4,edgecolor='k',facecolor=purp3) )
ax.add_patch( Rectangle([PR_ever,2.3],PR_never,0.4,edgecolor='k',facecolor=grey3) )
ax.vlines(PR_ever,2.3,2.7,linewidth=3,edgecolor='k')

ax.text (1.06,2,'p=0.6',fontsize=12,fontname=font)
ax.text (1.11,1,'p=0.002',fontsize=12,fontname=font)

# Legend
ax.add_patch( Rectangle([1.07,2.78],0.06,0.2,facecolor=grey3) )
ax.text (1.15,2.82,'no alteration',fontsize=12,fontname=font)
#ax.text (1.15,2.82,'never',fontsize=12,fontname=font)

ax.add_patch( Rectangle([1.07,2.56],0.02,0.2,facecolor=purp3) )
ax.add_patch( Rectangle([1.09,2.56],0.02,0.2,facecolor=blue3) )
ax.add_patch( Rectangle([1.11,2.56],0.02,0.2,facecolor=orange3) )
ax.text (1.15,2.61,'mutation',fontsize=12,fontname=font)
#ax.text (1.15,2.61,'ever',fontsize=12,fontname=font)

ax.vlines(1.03,1.5,2.5)
ax.hlines(2.5,1.01,1.03)
ax.hlines(1.5,1.01,1.03)

ax.vlines(1.08,0.5,1.5)
ax.hlines(1.5,1.06,1.08)
ax.hlines(0.5,1.06,1.08)

#plt.savefig('SD1_F_bar_smokingnever_SDPR.png',dpi=300)
#plt.savefig('SD1_F_bar_smokingnever_SDPR.eps')
#plt.savefig('SD1_F_bar_stk11_SDPR.png',dpi=300)
plt.savefig('SD1_F_bar_stk11_SDPR.eps')
