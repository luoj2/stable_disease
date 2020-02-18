
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import publication_settings
import pandas as pd

from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection

matplotlib.rcParams.update(publication_settings.params)

# ** define CSV file **
# Stable disease data load
data = pd.read_csv('SD2_F2_BOR_data.csv')

BOR = np.array(data['BOR'])
PFS = np.array(data['PFS'])
censor = np.array(data['Censor'])
OS = np.array(data['OS'])
OScensor = np.array(data['OSCensor'])

# Stable disease data load
data = pd.read_csv('SD2_F2_BOR_alldata.csv')

#BORall = np.array(data['BOR'])
PFSall = np.array(data['PFS'])
censorall = np.array(data['Censor'])
OSall = np.array(data['OS'])
OScensorall = np.array(data['OSCensor'])

# Partial response/ complete response data load
data = pd.read_csv('SD2_F2_BOR_PRCR3040_data.csv')

BOR_PR = np.array(data['BOR'])
PFS_PR = np.array(data['PFS'])
censor_PR = np.array(data['Censor'])
OS_PR = np.array(data['OS'])
OScensor_PR = np.array(data['OSCensor'])

data = pd.read_csv('SD2_F2_BOR_PRCR_data.csv')

BOR_aPR = np.array(data['BOR'])
PFS_aPR = np.array(data['PFS'])
censor_aPR = np.array(data['Censor'])
OS_aPR = np.array(data['OS'])
OScensor_aPR = np.array(data['OSCensor'])

# Progressive disease data load
data = pd.read_csv('SD2_F2_BOR_PD_data.csv')

#BOR_PD = np.array(data['BOR'])
PFS_PD = np.array(data['PFS'])
censor_PD = np.array(data['Censor'])
OS_PD = np.array(data['OS'])
OScensor_PD = np.array(data['OSCensor'])

# define SD-DCB as >6 months
mask = (PFS > 6)
DCB_BOR = BOR[mask]
DCB_PFS = PFS[mask]
DCB_censor = censor[mask]
DCB_OS = OS[mask]
DCB_OScensor = OScensor[mask]

# define SD-NDB as <=6 months
mask8 = (PFS <= 6)
NDB_BOR = BOR[mask8]
NDB_PFS = PFS[mask8]
NDB_censor = censor[mask8]
NDB_OS = OS[mask8]
NDB_OScensor = OScensor[mask8]

# define NDB as <=6 months
mask2 = (PFSall <= 6)
#NDBall_BOR = BOR[mask2]
NDBall_PFS = PFSall[mask2]
NDBall_censor = censorall[mask2]
NDBall_OS = OSall[mask2]
NDBall_OScensor = OScensorall[mask2]

# define DCB as >6 months
mask4 = (PFSall > 6)
#DCBall_BOR = BOR[mask4]
DCBall_PFS = PFSall[mask4]
DCBall_censor = censorall[mask4]
DCBall_OS = OSall[mask4]
DCBall_OScensor = OScensorall[mask4]

# define SD -0 to -30 percent AND PFS>6mo
mask3 = (BOR <=0) & (PFS>6)
SD20P3_BOR = BOR[mask3]
SD20P3_PFS = PFS[mask3]
SD20P3_censor = censor[mask3]
SD20P3_OS = OS[mask3]
SD20P3_OScensor = OScensor[mask3]

# define PRCR <= -57 median
mask8 = (BOR_aPR <=-57)
PRmedLess_BOR = BOR_aPR[mask8]
PRmedLess_PFS = PFS_aPR[mask8]
PRmedLess_censor = censor_aPR[mask8]
PRmedLess_OS = OS_aPR[mask8]
PRmedLess_OScensor = OScensor_aPR[mask8]

# define PRCR > -57 median
mask9 = (BOR_aPR >-57)
PRmedGreat_BOR = BOR_aPR[mask9]
PRmedGreat_PFS = PFS_aPR[mask9]
PRmedGreat_censor = censor_aPR[mask9]
PRmedGreat_OS = OS_aPR[mask9]
PRmedGreat_OScensor = OScensor_aPR[mask9]

# ** define months 0, 0.1, 0.2... max(PFS) **
#months = np.linspace(0,np.max(PFS),int(10*np.max(PFS))+1)
months = np.linspace(0,np.max(OS),int(10*np.max(OS))+1)
prob_of_survival = np.zeros(len(months))
num_at_risk = np.zeros(len(months))

def KaplanMeier(aDuration,aCensor):
	
	mask_censor = (aCensor==0)
	censored_duration = aDuration[mask_censor]

	num_alive = len(aCensor)
	prob_of_survival[0] = 1.

	for j,month in enumerate(months):
		num_died = np.sum( (np.abs(aDuration-month)<0.01) & (aCensor==1) )
		num_censored = np.sum( (np.abs(aDuration-month)<0.01) & (aCensor==0) )

		if num_died == 0 and j>0:
			prob_of_survival[j] = prob_of_survival[j-1]
			num_alive -= num_censored

		if num_died > 0 and j>0:
		
			prob_of_survival[j] = prob_of_survival[j-1]*((num_alive-num_died)
									/(num_alive))
			num_alive = num_alive - num_died - num_censored
			
		num_at_risk[j] = num_alive

	return months,prob_of_survival,censored_duration,len(aCensor),num_at_risk,aDuration,aCensor

def hazard_ratio_censor(group1,group2,group1censor,group2censor):

	#change PFS_censor != 0 if calc PFS hazard ratios
	deaths1 = np.sum(group1censor)
	deaths2 = np.sum(group2censor)

	mo = np.arange(0,np.ceil(max(np.max(group1),np.ceil(np.max(group2)))))

	expected_deaths1 = 0
	expected_deaths2 = 0

	for i,month in enumerate(mo):
		deaths1_mo_x = np.sum( (group1 < (month+1)) & (group1 >= month) & (group1censor == 1) )
		alive1_mo_x = np.sum( group1 >= month )
		deaths2_mo_x = np.sum( (group2 < (month+1)) & (group2 >= month) & (group2censor == 1) )
		alive2_mo_x = np.sum( group2 >= month )

		prob_death_mo_x = (deaths1_mo_x+deaths2_mo_x)/(alive1_mo_x+alive2_mo_x)

		expected_deaths1 += prob_death_mo_x*alive1_mo_x
		expected_deaths2 += prob_death_mo_x*alive2_mo_x

	hazard_ratio=(deaths1/expected_deaths1)/(deaths2/expected_deaths2)

	se = np.sqrt(1/expected_deaths1 + 1/expected_deaths2)
	hazard_ratio_m = hazard_ratio*np.exp(-1.96*se)
	hazard_ratio_p = hazard_ratio*np.exp(+1.96*se)

	return hazard_ratio_m, hazard_ratio, hazard_ratio_p

# ******* DESIGN *******
# ******* DESIGN *******
# 8.5 x 5 inches
total_width = 5
total_height = 4.5
fig = plt.figure(1, figsize=(total_width,total_height))

# give all in inches
top_border = 0.2
bottom_border = 1.4
left_border = 1.4
right_border = 0.2

left = left_border/total_width
bottom = bottom_border/total_height
width = (total_width - left_border - right_border)/total_width
height = (total_height - top_border - bottom_border)/total_height

ax = fig.add_axes([left,bottom,width,height])

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()

y_max = 1.01
x_max = 72
ax.set_xlim([0,48])
ax.set_ylim([0,y_max])

ax.set_yticks([0,.25,.50,.75,1])
ax.set_yticklabels((r"0%",r"25%",r"50%",r"75%",r"100%"))

ax.set_xticks([0,12,24,36,48])

font = 'Arial'

ax.set_ylabel(r"overall survival (%)",fontname=font)
ax.set_xlabel('time since treatment (months)',fontname=font)

grey = '#919191' # grey

green = '#CDDCAB' # green
coral = '#FCCBA1' # coral
blue = '#9ecae1' # blue
grblue = '#3182bd' # bio responder dk blue
ltgold = '#f7dc6f' # light gold

turq = '#06735F' # turquoise
tang = '#fdbb84' # tangerine
redor = '#D96248' # red orange
red = '#A64444' # red

dkpur = '#52428C' # dark purple
purp = '#54278f'
ltpurp = '#bcbddc' # light purple

# ******* DESIGN *******
# ******* DESIGN *******

ax.text(-18,-0.14, 'No. at risk', fontweight='bold', horizontalalignment="left"
	,verticalalignment="center",fontsize=10,fontname=font,clip_on=False)

nor_times = np.array([0,12,24,36,48])
nor_times10 = nor_times*10

#return months,prob_of_survival,censored_duration,len(aCensor),num_at_risk,aDuration,aCensor
# graph KM curve SD
m_SD,pos_SD,c_PFS_SD,t,n_SD,dur_SD,censor1 = KaplanMeier(OS,OScensor)
ax.plot(m_SD, pos_SD,color=blue,linewidth=1.5)
		
for value in c_PFS_SD:
	i = np.argmin(np.abs(months-value))
	ax.scatter(value,pos_SD[i],marker='|',color=blue,zorder=4,s=20)

n_SD_times = n_SD[nor_times10]

print(nor_times, nor_times10, n_SD_times)

ax.text(-18,-0.3, "SD", color=blue,fontweight='bold', horizontalalignment="left"
	,verticalalignment="center",fontsize=10,fontname=font,clip_on=False)

for i,(nor_times_i,n_SDt_i) in enumerate(zip(nor_times,n_SD_times)):
  ax.text(nor_times_i,-0.3, int(n_SDt_i), color=blue, horizontalalignment="center"
  	,verticalalignment="center",fontsize=10,fontname=font,clip_on=False)

# graph KM curve PD
m_PD,pos_PD,c_PFS_PD,t,n_PD,dur_PD,censor1 = KaplanMeier(OS_PD,OScensor_PD)
ax.plot(m_PD, pos_PD,color=tang,linewidth=1.5)
		
for value in c_PFS_PD:
	i = np.argmin(np.abs(months-value))
	ax.scatter(value,pos_PD[i],marker='|',color=tang,zorder=4,s=20)

n_PD_times = n_PD[nor_times10]

ax.text(-18,-0.34, "PD", color=tang,fontweight='bold', horizontalalignment="left"
	,verticalalignment="center",fontsize=10,fontname=font,clip_on=False)

for i,(nor_times_i,n_PDt_i) in enumerate(zip(nor_times,n_PD_times)):
  ax.text(nor_times_i,-0.34, int(n_PDt_i), color=tang, horizontalalignment="center"
  	,verticalalignment="center",fontsize=10,fontname=font,clip_on=False)

# graph KM curve SD 0 to -30 AND PFS >6months
m_SD20P3,pos_SD20P3,c_PFS_SD20P3,t,n_SD20P3,dur_SD20P3,censor1 = KaplanMeier(SD20P3_OS,SD20P3_OScensor)
ax.plot(m_SD20P3, pos_SD20P3,color=grblue,linewidth=1.5,zorder=4)
		
for value in c_PFS_SD20P3:
	i = np.argmin(np.abs(months-value))
	ax.scatter(value,pos_SD20P3[i],marker='|',color=grblue,zorder=4,s=20)

n_SD20P3_times = n_SD20P3[nor_times10]

ax.text(-18,-0.26, "SD responder", color=grblue,fontweight='bold', horizontalalignment="left"
	,verticalalignment="center",fontsize=10,fontname=font,clip_on=False)

for i,(nor_times_i,n_SD20P3t_i) in enumerate(zip(nor_times,n_SD20P3_times)):
  ax.text(nor_times_i,-0.26, int(n_SD20P3t_i), color=grblue, horizontalalignment="center"
  	,verticalalignment="center",fontsize=10,fontname=font,clip_on=False)

# graph KM curve PRCR <= -57 median
m_PRmedLess,pos_PRmedLess,c_PFS_PRmedLess,t,n_PRmedLess,dur_PRmedLess,censor1 = KaplanMeier(PRmedLess_OS,PRmedLess_OScensor)
ax.plot(m_PRmedLess, pos_PRmedLess,color=ltpurp,linewidth=1.5)
		
for value in c_PFS_PRmedLess:
	i = np.argmin(np.abs(months-value))
	ax.scatter(value,pos_PRmedLess[i],marker='|',color=ltpurp,zorder=4,s=20)

n_PRmedLess_times = n_PRmedLess[nor_times10]

ax.text(-18,-0.18, "PR major", color=ltpurp,fontweight='bold', horizontalalignment="left"
	,verticalalignment="center",fontsize=10,fontname=font,clip_on=False)

for i,(nor_times_i,n_PRmedLesst_i) in enumerate(zip(nor_times,n_PRmedLess_times)):
  ax.text(nor_times_i,-0.18, int(n_PRmedLesst_i), color=ltpurp, horizontalalignment="center"
  	,verticalalignment="center",fontsize=10,fontname=font,clip_on=False)

# graph KM curve PRCR > -57 median
m_PRmedGreat,pos_PRmedGreat,c_PFS_PRmedGreat,t,n_PRmedGreat,dur_PRmedGreat,censor1 = KaplanMeier(PRmedGreat_OS,PRmedGreat_OScensor)
ax.plot(m_PRmedGreat, pos_PRmedGreat,color=purp,linewidth=1.5,zorder=4)
		
for value in c_PFS_PRmedGreat:
	i = np.argmin(np.abs(months-value))
	ax.scatter(value,pos_PRmedGreat[i],marker='|',color=purp,zorder=4,s=20)

n_PRmedGreat_times = n_PRmedGreat[nor_times10]

ax.text(-18,-0.22, "PR minor", color=purp,fontweight='bold', horizontalalignment="left"
	,verticalalignment="center",fontsize=10,fontname=font,clip_on=False)

for i,(nor_times_i,n_PRmedGreatt_i) in enumerate(zip(nor_times,n_PRmedGreat_times)):
  ax.text(nor_times_i,-0.22, int(n_PRmedGreatt_i), color=purp, horizontalalignment="center"
  	,verticalalignment="center",fontsize=10,fontname=font,clip_on=False)


#plt.savefig('SD2_F2_KMcompare_draft11_OS.png',dpi=300)
plt.savefig('SD2_F2_KMcompare_draft11_OS.eps')