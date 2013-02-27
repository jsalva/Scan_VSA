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

condition={'RD':1,'SD':2,'RI':3, 1:'RD',2:'SD',3:'RI'}

time_to_target=[]
endpoint_error=[]
accuracy=[]

chunk_size = 8


for subject in subject_list:

    sequence_labels = []
    sequence_data_file = data[subject][0][0][0]
    sequence_data = np.genfromtxt(sequence_data_file,skip_header=True)
    sequence_data = sequence_data[sequence_data[:,1]!=0]
    conditions = sequence_data[:,1]

    tmp_time_to_target = []
    tmp_endpoint_error = []
    tmp_accuracy = []

    for chunk in range(len(sequence_data)/chunk_size):
        

        chunk_data = sequence_data[chunk*chunk_size:(chunk+1)*chunk_size]
        chunk_cond = condition[np.mean(chunk_data[:,1])]
        chunk_correct = chunk_data[:,10] == 1
        chunk_time_to_target = chunk_data[chunk_correct,4]
        chunk_endpoint_error = chunk_data[chunk_correct,5]
        chunk_accuracy = chunk_data[:,5]
        
        sequence_labels.append(chunk_cond)

        tmp_time_to_target.append(np.median(chunk_time_to_target))
        tmp_endpoint_error.append(np.median(chunk_endpoint_error))
        tmp_accuracy.append(np.mean(chunk_accuracy))

    try:
        time_to_target = np.vstack([time_to_target,tmp_time_to_target])
    except:
        time_to_target = tmp_time_to_target
    try:
        endpoint_error = np.vstack([endpoint_error,tmp_endpoint_error])
    except:
        endpoint_error = tmp_endpoint_error
    try:
        accuracy = np.vstack([accuracy,tmp_accuracy])
    except:
        accuracy = tmp_accuracy

fig = plt.figure()
ax = fig.add_subplot(111)
ind = np.arange(len(sequence_labels))
width = 0
linewidth = 0
rects = ax.bar(ind, 
     np.mean(time_to_target,axis=0), 
     width=width,
     linewidth=linewidth, 
     color = 'b',
     yerr = np.std(time_to_target)/len(subject_list),
     ecolor = 'k')

for idx,mean in enumerate(np.mean(time_to_target,axis=0)):
    color = 'blue' if sequence_labels[idx] == 'SD' else 'black'
    ax.plot(idx,mean,'o',linestyle='None',color = color)


#fix subsequently
for idx,val in enumerate(sequence_labels):
    try:
        if curr_val == val:
                print "y"
                medians.append(np.median(time_to_target,axis=0)[idx])
        else:
                print "n"
                curr_val = val
                end_idx = idx
                ax.plot(np.arange(start_idx,end_idx),np.repeat(np.median(medians),end_idx - start_idx))
                medians = []
                medians.append(np.median(time_to_target,axis=0)[idx])
                start_idx = idx
         
    except:
        medians = []
        curr_val = val
        start_idx = idx
        medians.append(np.median(time_to_target,axis=0)[idx])


ax.set_ylabel('Time (s)')
ax.set_title('Mean Time to Target')
ax.set_xticks(ind+width)
ax.set_xticklabels(sequence_labels)
#autolabel(rects)
plt.show()


