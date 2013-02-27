import numpy as np
from scipy.ndimage import label
import matplotlib.pyplot as plt

def autolabel(rects):
    for rect in rects:
        height = rect.get_height()
        ax.text(rect.get_x()+rect.get_width()/2., 1.05*height, '%d'%int(height),
                ha='center', va='bottom')

subject_list = ['lap002','lap004','lap005','lap006']

data = {}
data['lap002'] = [(['lap002_run1_cb1_summary.txt'],'vsa_sequence'),(['lap002_run2_cb1_summary.txt'],'vsa_perturb')]
data['lap004'] = [(['lap004_run2_cb2_summary.txt'],'vsa_sequence'),(['lap004_run1_cb2_summary.txt'],'vsa_perturb')]
data['lap005'] = [(['lap005_run1_cb1_summary.txt'],'vsa_sequence'),(['lap005_run2_cb1_summary.txt'],'vsa_perturb')]
data['lap006'] = [(['lap006_run2_cb2_summary.txt'],'vsa_sequence'),(['lap006_run1_cb2_summary.txt'],'vsa_perturb')]

condition={'RD':1,'SD':2,'RI':3}

sd = None
sd_time_to_target=[]
sd_endpoint_error=[]
sd_accuracy=[]

rd = None
rd_time_to_target=[]
rd_endpoint_error=[]
rd_accuracy=[]

for subject in subject_list:
    print subject
    sequence_data_file = data[subject][0][0][0]
    sequence_data = np.genfromtxt(sequence_data_file,skip_header=True)
    sequence_data = sequence_data[sequence_data[:,1]!=0]
    conditions = sequence_data[:,1]

    sd = [conditions == condition['SD']]
    plt.scatter(sequence_data[:,2][sd],sequence_data[:,4][sd],c='b')
    sd_time_to_target.append(np.median(sequence_data[:,4][sd]))
    sd_endpoint_error.append(np.median(np.abs(sequence_data[:,5][sd])))
    sd_accuracy.append(np.mean(sequence_data[:,10][sd]))

    rd = [conditions == condition['RD']]
    plt.scatter(sequence_data[:,2][rd],sequence_data[:,4][rd],c='k')
    rd_time_to_target.append(np.median(sequence_data[:,4][rd]))
    rd_endpoint_error.append(np.median(np.abs(sequence_data[:,5][rd])))
    rd_accuracy.append(np.mean(sequence_data[:,10][rd]))

plt.show()

fig = plt.figure()
ax = fig.add_subplot(111)
ind = np.arange(2)
width = 0.35
rects = ax.bar(ind, 
     [np.mean(sd_time_to_target),np.mean(rd_time_to_target)], 
     width, 
     color = 'b',
     yerr = [np.std(sd_time_to_target)/np.sqrt(len(subject_list)), np.std(rd_time_to_target)/np.sqrt(len(subject_list))],
     ecolor = 'k',
     label='SD')

ax.set_ylabel('Time (s)')
ax.set_title('Mean Time to Target')
ax.set_xticks(ind+width)
ax.set_xticklabels( ('SD','RD') )
#autolabel(rects)
plt.show()


fig = plt.figure()
ax = fig.add_subplot(111)
ind = np.arange(2)
width = 0.35
rects = ax.bar(ind, 
     [np.mean(np.abs(sd_endpoint_error)),np.mean(np.abs(rd_endpoint_error))], 
     width, 
     color = 'b',
     yerr = [np.std(np.abs(sd_endpoint_error))/np.sqrt(len(subject_list)), np.std(np.abs(rd_endpoint_error))/np.sqrt(len(subject_list))],
     ecolor = 'k',
     label='SD')

ax.set_ylabel('Radians')
ax.set_title('Mean Endpoint Error')
ax.set_xticks(ind+width)
ax.set_xticklabels( ('SD','RD') )
#autolabel(rects)
plt.show()
     
fig = plt.figure()
ax = fig.add_subplot(111)
ind = np.arange(2)
width = 0.35
rects = ax.bar(ind, 
     [np.mean(sd_accuracy),np.mean(rd_accuracy)], 
     width, 
     color = 'b',
     yerr = [np.std(sd_accuracy)/np.sqrt(len(subject_list)), np.std(rd_accuracy)/np.sqrt(len(subject_list))],
     ecolor = 'k',
     label='SD')

ax.set_ylabel('% within error bounds')
ax.set_title('Mean Accuracy')
ax.set_xticks(ind+width)
ax.set_xticklabels( ('SD','RD') )
#autolabel(rects)
plt.show()

ri = None
ri_time_to_target=[]
ri_endpoint_error=[]
ri_accuracy=[]

rd = None
rd_time_to_target=[]
rd_endpoint_error=[]
rd_accuracy=[]

for subject in subject_list:
    print subject
    perturb_data_file = data[subject][1][0][0]
    perturb_data = np.genfromtxt(perturb_data_file,skip_header=True)
    perturb_data = perturb_data[perturb_data[:,1]!=0]
    conditions = perturb_data[:,1]

    ri = [conditions == condition['RI']]
    plt.scatter(perturb_data[:,4][ri],perturb_data[:,4][ri],c='r')
    ri_time_to_target.append(np.median(perturb_data[:,4][ri]))
    ri_endpoint_error.append(np.median(np.abs(perturb_data[:,5][ri])))
    ri_accuracy.append(np.mean(perturb_data[:,10][ri]))

    rd = [conditions == condition['RD']]
    plt.scatter(perturb_data[:,2][rd],perturb_data[:,4][rd],c='k')
    rd_time_to_target.append(np.median(perturb_data[:,4][rd]))
    rd_endpoint_error.append(np.median(np.abs(perturb_data[:,5][rd])))
    rd_accuracy.append(np.mean(perturb_data[:,10][rd]))

plt.show()

fig = plt.figure()
ax = fig.add_subplot(111)
ind = np.arange(2)
width = 0.35
rects = ax.bar(ind, 
     [np.mean(ri_time_to_target),np.mean(rd_time_to_target)], 
     width, 
     color = 'r',
     yerr = [np.std(ri_time_to_target)/np.sqrt(len(subject_list)), np.std(rd_time_to_target)/np.sqrt(len(subject_list))],
     ecolor = 'k',
     label='RI')

ax.set_ylabel('Time (s)')
ax.set_title('Mean Time to Target')
ax.set_xticks(ind+width)
ax.set_xticklabels( ('RI','RD') )
#autolabel(rects)
plt.show()


fig = plt.figure()
ax = fig.add_subplot(111)
ind = np.arange(2)
width = 0.35
rects = ax.bar(ind, 
     [np.mean(np.abs(ri_endpoint_error)),np.mean(np.abs(rd_endpoint_error))], 
     width, 
     color = 'r',
     yerr = [np.std(np.abs(ri_endpoint_error))/np.sqrt(len(subject_list)), np.std(np.abs(rd_endpoint_error))/np.sqrt(len(subject_list))],
     ecolor = 'k',
     label='RI')

ax.set_ylabel('Radians')
ax.set_title('Mean Endpoint Error')
ax.set_xticks(ind+width)
ax.set_xticklabels( ('RI','RD') )
#autolabel(rects)
plt.show()
     
fig = plt.figure()
ax = fig.add_subplot(111)
ind = np.arange(2)
width = 0.35
rects = ax.bar(ind, 
     [np.mean(ri_accuracy),np.mean(rd_accuracy)], 
     width, 
     color = 'r',
     yerr = [np.std(ri_accuracy)/np.sqrt(len(subject_list)), np.std(rd_accuracy)/np.sqrt(len(subject_list))],
     ecolor = 'k',
     label='RI')

ax.set_ylabel('% within error bounds')
ax.set_title('Mean Accuracy')
ax.set_xticks(ind+width)
ax.set_xticklabels( ('RI','RD') )
#autolabel(rects)
plt.show()

    
    #repeat for perturb
