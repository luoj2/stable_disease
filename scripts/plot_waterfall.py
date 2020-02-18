
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from readdata_F1b import load
import publication_settings
import pandas as pd

from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection

matplotlib.rcParams.update(publication_settings.params)

(number_list,change_list,censor_list) = load()

data = pd.read_csv('SD2_F1b_rawsort.csv')

number_list = np.array(data['Index'])
change_list = np.array(data['BOR'])
censor_list = np.array(data['Censor'])

# 11 x 5 inches
total_width = 12
total_height = 3.5
fig = plt.figure(1, figsize=(total_width,total_height))

# give all in inches
top_border = 0.5
bottom_border = 0.5
left_border = 1
right_border = 0.25

left = left_border/total_width
bottom = bottom_border/total_height
width = (total_width - left_border - right_border)/total_width
height = (total_height - top_border - bottom_border)/total_height

ax = fig.add_axes([left,bottom,width,height])

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()
ax.axhline(y=0,xmin=0,xmax=1,color='k',zorder=4)
ax.axhline(y=20,xmin=0,xmax=1,linestyle='--',color='k',zorder=0)
ax.axhline(y=-30,xmin=0,xmax=1,linestyle='--',color='k',zorder=0)

num_entries = len(number_list)

y_max = 32
y_min = -32
x_max = 300
ax.set_xlim([0,num_entries])
ax.set_ylim([y_min,y_max])

y_unit = total_height*height/y_max
x_unit = total_width*width/num_entries

ax.set_yticks([30,20,10,0,-10,-20,-30])
ax.set_yticklabels((r"30%",r"20%",r"10%",r"0%",r"$-10$%",r"$-20$%",r"$-30$%"),fontsize=14)

font = 'Arial'

ax.set_ylabel(r"BOR (% change from baseline)",fontname=font,fontsize=14)
ax.set_xlabel('patient number',fontname=font,fontsize=14)

grey = '#919191' # grey
green = '#CDDCAB' # green
coral = '#FCCBA1' # coral
blue = '#6baed6' # blue
ltgold = '#f7dc6f'

colors_bar = [grey,blue]

for i,(number,change,censor) in enumerate(zip(number_list,change_list,censor_list)):

  # make box
  ax.add_patch( Rectangle([i+0.1,0],0.8,change,facecolor=colors_bar[censor], zorder=3) )

# make legend
legend_x = 240
legend_y = 28
legend_dy = 3
legend_dx = legend_dy/num_entries*(height*total_height)*x_max/(width*total_width)
#
delta_y = 6
colors = [blue,grey]
labels = ["not censored (n=281)","censored (n=32)"]
for i,(color,label) in enumerate(zip(colors,labels)):
  ax.add_patch( Rectangle([legend_x,legend_y-i*delta_y],legend_dx+1,legend_dy,facecolor=color,edgecolor='w') )
  ax.text( legend_x + legend_dx + 2,legend_y - i*delta_y + legend_dy/2.-0.1, label,
   horizontalalignment="left",verticalalignment="center",fontsize=14,fontname=font)

#plt.savefig('SD2_F1b_waterfall.png',dpi=300)
plt.savefig('SD2_F1b_waterfall.eps')

