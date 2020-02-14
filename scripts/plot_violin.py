import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from readdata_F1a_swim import load
import publication_settings

from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection

#seaborn
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats

matplotlib.rcParams.update(publication_settings.params)

data = pd.read_csv('SD2_F2_BOR_data.csv')
PFS_list = np.array(data['PFS'])
censor_list = np.array(data['Censor'])

# 8.5 x 5 inches
total_width = 8.5
total_height = 5
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

num_entries = len(censor_list)

y_max = 70

y_unit = total_height*height/y_max
x_unit = total_width*width/num_entries

font = 'Arial'

grey = '#919191' # grey
green = '#CDDCAB' # green
coral = '#FCCBA1' # coral
blue = '#A3BFD9' # blue

x=PFS_list
sns.violinplot(x=PFS_list, ax=ax, inner="quartile", color=blue)

ax.set_xlabel('months',fontname=font,fontsize=14)
ax.set_xlim([0,54])

ax.set_xticks([0,6,12,18,24,30,36,42,48,54])
ax.set_xticklabels((r"0",r"",r"12",r"",r"24",r"",r"36",r"",r"48",r"",r"60"),fontsize=14)

ax.text(3,-0.425,"25%",fontsize=14,fontname=font)
ax.text(5.9,-0.38,"50%",fontsize=14,fontname=font)
ax.text(9.6,-0.2,"75%",fontsize=14,fontname=font)

#plt.savefig('SD2_F1a_violin.png',dpi=300)
plt.savefig('SD2_F1a_violin.eps')
