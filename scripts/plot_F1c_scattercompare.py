
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import publication_settings
import pandas as pd

from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection

matplotlib.rcParams.update(publication_settings.params)

data = pd.read_csv('SD2_F2_BOR_PRCR_data.csv')

BOR_listPR = np.array(data['BOR'])
PFS_listPR = np.array(data['PFS'])

# define PRCR <= -57 median
# Less is -57 to -100
mask = (BOR_listPR <=-57)
Less_BOR = BOR_listPR[mask]
Less_PFS = PFS_listPR[mask]

# define PRCR > -57 median
# Great is -57 to -30
mask2 = (BOR_listPR >-57)
Great_BOR = BOR_listPR[mask2]
Great_PFS = PFS_listPR[mask2]

data = pd.read_csv('SD2_F2_BOR_data.csv')

BOR_list = np.array(data['BOR'])
PFS_list = np.array(data['PFS'])

# 8 x 8 inches
total_width = 8
total_height = 6
fig = plt.figure(1, figsize=(total_width,total_height))

# give all in inches
top_border = 0.5
bottom_border = 0.5
left_border = 1.2
right_border = 0.5

left = left_border/total_width
bottom = bottom_border/total_height
width = (total_width - left_border - right_border)/total_width
height = (total_height - top_border - bottom_border)/total_height

ax = fig.add_axes([left,bottom,width,height])

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()


y_max = 22
y_min = -101
x_min = 0
x_max = 48

y_unit = total_height*height/y_max
x_unit = total_width*width/x_max


ax.set_xlim([x_min,x_max])
ax.set_ylim([y_min,y_max])

ax.set_yticks([20,0,-20,-40,-60,-80,-100])
ax.set_yticklabels((r"20%",r"0%",r"$-20$%",r"$-40$%",
	r"$-60$%",r"$-80$%",r"$-100$%"),fontsize=14)

ax.set_xticks([0,12,24,36])
#ax.set_xticks([0,12,24,36,48,60,72,84])

font = 'Arial'

ax.set_ylabel(r"best overall response (%)",fontname=font,fontsize=14)
ax.set_xlabel('PFS (months)',fontname=font,fontsize=14)

blue = '#3182bd' # blue
purple = '#756bb1' #purple
purple2 = '#54278f' #purple
grey = '#919191' # grey

ax.scatter(PFS_list,BOR_list,color=blue,marker='o',s=30,edgecolors=grey,alpha=0.5)
ax.scatter(Great_PFS,Great_BOR,color=purple2,marker='o',s=30,edgecolors=grey,alpha=0.7)
ax.scatter(Less_PFS,Less_BOR,color=purple,marker='o',s=30,alpha=0.6)

#### calculate fit for SD
fit_params = np.polyfit(BOR_list,PFS_list,1)
fit = np.poly1d(fit_params)
PFS_fit = fit(BOR_list)
mean = np.mean(PFS_list)
error_mean = np.sum((PFS_list-mean)**2)
error_fit  = np.sum((PFS_list-PFS_fit)**2)
R2 = 1 - error_fit/error_mean
plt.plot(PFS_fit,BOR_list, color='k', linewidth=2)

R2_text = r"$R^2 = %.3f$" % R2
fit_text = r"$\mathrm{PFS} = %.2f \mathrm{BOR} + %.2f$" %(fit[1],fit[0])
print(R2_text)
print(fit_text)

ax.text(20,-18, R2_text, horizontalalignment="left",verticalalignment="center",
	fontsize=14,fontname=font)
ax.text(20,-23, fit_text, horizontalalignment="left",verticalalignment="center",
	fontsize=14,fontname=font)
#### calculate fit for SD

#### calculate fit for PR/CR
fit_params = np.polyfit(BOR_listPR,PFS_listPR,1)
fit = np.poly1d(fit_params)
PFS_fit = fit(BOR_listPR)
mean = np.mean(PFS_listPR)
error_mean = np.sum((PFS_listPR-mean)**2)
error_fit  = np.sum((PFS_listPR-PFS_fit)**2)
R2 = 1 - error_fit/error_mean
plt.plot(PFS_fit,BOR_listPR, color='k', linewidth=2)

R2_text = r"$R^2 = %.3f$" % R2
fit_text = r"$\mathrm{PFS} = %.2f \mathrm{BOR} + %.2f$" %(fit[1],fit[0])
print(R2_text)
print(fit_text)

ax.text(32,-88, R2_text, horizontalalignment="left",verticalalignment="center",
	fontsize=14,fontname=font)
ax.text(32,-93, fit_text, horizontalalignment="left",verticalalignment="center",
	fontsize=14,fontname=font)
#### calculate fit for PR/CR

### LEGEND
legend_x = 40
legend_y = 15
legend_dy = 1
legend_dx = 0.8
#
delta_y = 4
colors = [blue,purple]
labels = ["SD","PR/CR"]
for i,(color,label) in enumerate(zip(colors,labels)):
  #ax.add_patch( Rectangle([legend_x,legend_y-i*delta_y],legend_dx,legend_dy,facecolor=color,
  #	edgecolor='w') )
  ax.scatter(legend_x,legend_y-i*delta_y+1,color=color,marker='o',s=80,edgecolors=grey,alpha=0.5)
  ax.text( legend_x + legend_dx,legend_y - i*delta_y + legend_dy, label, 
  	horizontalalignment="left",verticalalignment="center",fontsize=14,fontname=font)


plt.savefig('SD2_F1c_compare.png',dpi=300)
#plt.savefig('SD1_F1c_compare.eps')
