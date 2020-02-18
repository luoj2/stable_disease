
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import publication_settings
import pandas as pd

from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection

matplotlib.rcParams.update(publication_settings.params)

labels = ['PD','SD responder','PR minor']

def array_fraction(array,condition):
    return np.sum(array==condition)/len(array)

data = pd.read_csv('SD2_SDfeatures.csv')
SD_feat = np.array(data['dNLR'])
#SD_feat = np.array(data['Line'])
#SD_feat = np.array(data['TMB'])
#SD_feat = np.array(data['PDL1'])
#SD_feat = np.array(data['Smoking'])

BOR_SD = np.array(data['BOR'])
PFS_SD = np.array(data['PFS'])
# define SD -0 to -30 percent AND PFS>6mo
mask1 = (BOR_SD <=0) & (PFS_SD>6)
SD_feat = SD_feat[mask1]
known = (SD_feat != 'Unknown')
SD_featfloat = SD_feat[known].astype(np.float64)

print(SD_featfloat)

data = pd.read_csv('SD2_PRfeatures.csv')
PR_feat = np.array(data['dNLR'])
#PR_feat = np.array(data['Line'])
#PR_feat = np.array(data['TMB'])
#PR_feat = np.array(data['PDL1'])
#PR_feat = np.array(data['Smoking'])

BOR_PR = np.array(data['BOR'])
# define PR minor
mask = BOR_PR >-57
PR_feat = PR_feat[mask]
#print(len(PR_feat))
known = (PR_feat != 'Unknown')
PR_featfloat = PR_feat[known].astype(np.float64)
#PR_featfloat = PR_feat.astype(np.float64)
#print(PR_feat)

data = pd.read_csv('PDfeatures.csv')
PD_feat = np.array(data['dNLR'])
#PD_feat = np.array(data['Line'])
#PD_feat = np.array(data['TMB'])
#PD_feat = np.array(data['PDL1'])
#PD_feat = np.array(data['Smoking'])

known = (PD_feat != 'Unknown')
PD_featfloat = PD_feat[known].astype(np.float64)

#define IQR, median
SD_q75, SD_median, SD_q25 = np.percentile(SD_featfloat, [75 ,50, 25])
print(SD_median, SD_q75, SD_q25)
PR_q75, PR_median, PR_q25 = np.percentile(PR_featfloat, [75 ,50, 25])
print(PR_q75)
PD_q75, PD_median, PD_q25 = np.percentile(PD_featfloat, [75 ,50, 25])

#define max and min, outliners
SD_max = 2.5*(SD_q75-SD_median)+SD_median
if SD_max > np.max(SD_featfloat):
	SD_max = np.max(SD_featfloat)
SD_min = SD_median-2.5*(SD_median-SD_q25)
if SD_min < np.min(SD_featfloat):
	SD_min = np.min(SD_featfloat)

#print(SD_max,SD_min)
SD_outlier_mask = (SD_featfloat > SD_max)
SD_outlier_max = SD_featfloat[SD_outlier_mask].astype(np.float64)
#print(SD_outlier_max)
SD_outlier_mask = (SD_featfloat < SD_min)
SD_outlier_min = SD_featfloat[SD_outlier_mask].astype(np.float64)
#print(SD_outlier_min)

PR_max = 2.5*(PR_q75-PR_median)+PR_median
if PR_max > np.max(PR_featfloat):
	PR_max = np.max(PR_featfloat)
print(PR_median, PR_max)
PR_min = PR_median-2.5*(PR_median-PR_q25)
if PR_min < np.min(PR_featfloat):
	PR_min = np.min(PR_featfloat)
PR_outlier_mask = (PR_featfloat > PR_max)
PR_outlier_max = PR_featfloat[PR_outlier_mask].astype(np.float64)
PR_outlier_mask = (PR_featfloat < PR_min)
PR_outlier_min = PR_featfloat[PR_outlier_mask].astype(np.float64)

PD_max = 2.5*(PD_q75-PD_median)+PD_median
if PD_max > np.max(PD_featfloat):
	PD_max = np.max(PD_featfloat)
PD_min = PD_median-2.5*(PD_median-PD_q25)
if PD_min < np.min(PD_featfloat):
	PD_min = np.min(PD_featfloat)
PD_outlier_mask = (PD_featfloat > PD_max)
PD_outlier_max = PD_featfloat[PD_outlier_mask].astype(np.float64)
PD_outlier_mask = (PD_featfloat < PD_min)
PD_outlier_min = PD_featfloat[PD_outlier_mask].astype(np.float64)

# 5 x 5 inches
total_width = 5
total_height = 2.5
fig = plt.figure(1, figsize=(total_width,total_height))

# give all in inches
top_border = 0.2
bottom_border = 0.5
left_border = 1.5
right_border = 0.3

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

#ax.set_xlim([0,10])
ax.set_xlim([0,7])
ax.set_ylim([0,len(labels)])

ax.tick_params('y',direction='out',which='both')
#ax.set_yticks(minor_ticks,minor=True)

font = 'Arial'

#ax.set_xticks([0,25,50,75,100])
#ax.set_xticks([0,2,4,6,8,10])
ax.set_xticks([0,2,4,6])
#ax.set_xticklabels((r"0%",r"25%",r"50%",r"75%",r"100%"))
y=np.array([0.1,0.5,1.1,1.5,2.1,2.5])
ax.set_yticks(y)
ax.vlines(0,0.5,2.5)
#ax.xaxis.set_tick_params(length=0)
#tick_labels = ('n=%i' %len(PD_featfloat), labels[0], 'n=%i' %len(PR_featfloat), labels[1], 'n=%i' %len(SD_featfloat), labels[2])
#tick_labels = ('n=%i' %len(PD_featfloat), labels[0], 'n=%i' %len(SD_featfloat), labels[1], 'n=%i' %len(PR_featfloat), labels[2])
tick_labels = ('n=%i' %len(PD_featfloat), labels[0], 'n=%i' %len(SD_featfloat), labels[1], 'n=118')
ax.set_yticklabels(tick_labels, fontsize=12)
plt.setp(ax.get_yticklines()[0],visible=False)
plt.setp(ax.get_yticklines()[4],visible=False)
plt.setp(ax.get_yticklines()[8],visible=False)

#ax.set_xlabel('pack-year',fontname=font,fontsize=12)

#SD
blue1 = '#eff3ff' #smoking ever never
#blue3 = '#c6dbef' #2 smoking pack year
blue3 = '#9ecae1' #dNLR
#blue3 = '#6baed6' #PDL1
#blue3 = '#3182bd' #5
blue6 = '#08519c'

#PR
orange1 = '#fef0d9'
#orange3 = '#fdd49e' #2
orange3 = '#fdbb84' #real
#orange3 = '#fc8d59'
#orange3 = '#e34a33' #5
orange6 = '#b30000'

#PD
purp1 = '#f2f0f7'
#purp3 = '#dadaeb' #2
purp3 = '#bcbddc' #real
#purp3 = '#9e9ac8'
#purp3 = '#756bb1' #5
purp6 = '#54278f'

np.random.seed(seed=42)
y_amp = 0.03

ax.add_patch( Rectangle([PD_q25,0.3],(PD_q75-PD_q25),0.4,edgecolor='k',linewidth='2',facecolor=orange3) )
ax.vlines(PD_median,0.3,0.7,lw=2)
ax.scatter(PD_min,0.5,color='k',marker='|',s=50,lw=2)
ax.scatter(PD_max,0.5,color='k',marker='|',s=50,lw=2)
plt.plot([PD_min,PD_max],[0.5,0.5],color='k',ls='-',lw=2,zorder=0)

for outlier in PD_outlier_min:
	ax.scatter(outlier,0.5 + y_amp*np.random.normal(),color=orange3,marker='o',s=10,lw=1,zorder=0,alpha=0.5)
for outlier in PD_outlier_max:
	ax.scatter(outlier,0.5 + y_amp*np.random.normal(),color=orange3,marker='o',s=10,lw=1,zorder=0,alpha=0.5)

ax.add_patch( Rectangle([SD_q25,1.3],(SD_q75-SD_q25),0.4,edgecolor='k',linewidth='2',facecolor=blue3) )
ax.vlines(SD_median,1.3,1.7,lw=2)
ax.scatter(SD_min,1.5,color='k',marker='|',s=50,lw=2)
ax.scatter(SD_max,1.5,color='k',marker='|',s=50,lw=2)
plt.plot([SD_min,SD_max],[1.5,1.5],color='k',ls='-',lw=2,zorder=0)

for outlier in SD_outlier_min:
	ax.scatter(outlier,1.5 + y_amp*np.random.normal(),color=blue3,marker='o',s=10,lw=1,zorder=0,alpha=0.5)
for outlier in SD_outlier_max:
	ax.scatter(outlier,1.5 + y_amp*np.random.normal(),color=blue3,marker='o',s=10,lw=1,zorder=0,alpha=0.5)

ax.add_patch( Rectangle([PR_q25,2.3],(PR_q75-PR_q25),0.4,edgecolor='k',linewidth='2',facecolor=purp3) )
ax.vlines(PR_median,2.3,2.7,lw=2)
ax.scatter(PR_min,2.5,color='k',marker='|',s=50,lw=2)
ax.scatter(PR_max,2.5,color='k',marker='|',s=50,lw=2)
plt.plot([PR_min,PR_max],[2.5,2.5],color='k',ls='-',lw=2,zorder=0)

for outlier in PR_outlier_min:
	ax.scatter(outlier,2.5 + y_amp*np.random.normal(),color=purp3,marker='o',s=10,lw=1,zorder=0,alpha=0.5)
for outlier in PR_outlier_max:
	ax.scatter(outlier,2.5 + y_amp*np.random.normal(),color=purp3,marker='o',s=10,lw=1,zorder=0,alpha=0.5)


ax.text (5.2,2,'p=0.5',fontsize=12,fontname=font)
ax.text (4.9,1,'p<0.0001',fontsize=12,fontname=font)

ax.vlines(5.,1.5,2.5)
ax.hlines(2.5,4.85,5)
ax.hlines(1.5,4.85,5)

ax.vlines(6.55,0.5,1.5)
ax.hlines(1.5,6.45,6.55)
ax.hlines(0.5,6.45,6.55)

plt.savefig('SD1_F_bar_dNLR_SDPR.png',dpi=300)
#plt.savefig('SD2_F_bar_line_SDPR.png',dpi=300)
#plt.savefig('SD2_F_bar_TMB_SDPR.png',dpi=300)
#plt.savefig('SD2_F_bar_PDL1_SDPR.png',dpi=300)
#plt.savefig('SD2_F_bar_PDL1_SDPR.eps')
#plt.savefig('SD2_F_bar_smoking_SDPR.png',dpi=300)


