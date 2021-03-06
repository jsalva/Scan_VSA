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
    fig = plt.figure()
    ax = fig.add_subplot(111)

    sequence_labels = []
    sequence_data_file = data[subject][0][0][0]
    sequence_data = np.genfromtxt(sequence_data_file,skip_header=True)
    sequence_data = sequence_data[sequence_data[:,1]!=0]
    conditions = sequence_data[:,1]
    tmp_time_to_target = []

    for chunk in range(len(sequence_data)/chunk_size):
            

        chunk_data = sequence_data[chunk*chunk_size:(chunk+1)*chunk_size]
        chunk_cond = condition[np.mean(chunk_data[:,1])]
        chunk_correct = chunk_data[:,10] == 1
        chunk_time_to_target = chunk_data[chunk_correct,4]
        ax.plot(np.arange(chunk*chunk_size,(chunk+1)*chunk_size)/chunk_size,chunk_data[:,4],'o',color='black',alpha='0.2',linestyle='None')
            
        sequence_labels.append(chunk_cond)

        time_to_target.append(chunk_time_to_target)

    median_time_to_target = [np.median(times) for times in time_to_target]
    ind = np.arange(len(sequence_labels))
    width = 0
    linewidth = 0
    rects = ax.bar(ind, 
         median_time_to_target, 
         width=width,
         linewidth=linewidth, 
         color = 'b',
         yerr = [np.std(times)/np.sqrt(len(times)) for times in time_to_target],
         ecolor = 'k')

    for idx,median in enumerate(median_time_to_target):
        color = 'blue' if sequence_labels[idx] == 'SD' else 'red'
        ax.plot(idx,median,'o',linestyle='None',color = color)


    #fix subsequently
    medians = []
    start_idx = None
    end_idx = None
    curr_val = None
    for idx,val in enumerate(sequence_labels):
        try:
            print idx
            if curr_val == val:
                print "y"
                medians.append(median_time_to_target[idx])
                if idx == len(sequence_labels) - 1:
                    print "DONE!"
                    curr_val = val
                    end_idx = idx+1
                    print start_idx, end_idx
                    ax.plot(np.arange(start_idx,end_idx),np.repeat(np.mean(medians),end_idx - start_idx),linewidth=4,color='black')

                    medians = []
                    medians.append(median_time_to_target[idx])
                    start_idx = idx
                
            
            else:
                    print "n"
                    curr_val = val
                    end_idx = idx
                    print start_idx, end_idx
                    ax.plot(np.arange(start_idx,end_idx),np.repeat(np.mean(medians),end_idx - start_idx),linewidth=4,color='black')
                    medians = []
                    medians.append(median_time_to_target[idx])
                    start_idx = idx
             
        except:
            medians = []
            curr_val = val
            start_idx = idx
            print "start idx!: %d"%(idx)
            medians.append(median_time_to_target[idx])

    skill_scores = []
    indexes = []
    for i in [0,1,2,3]:
        rd_block_1 = sequence_data[(48+16)*i:(48+16)*i+16,4]
        sd = sequence_data[(48+16)*i+16:(48+16)*i+48,4]
        rd_block_2 = sequence_data[(48+16)*i+48:(48+16)*i+48+16,4]
        idx_i = (48+16)*i+48
        indexes.append(idx_i)
        rd = np.hstack([rd_block_1,rd_block_2])
        skill_score = (np.median(rd)-np.median(sd))/np.median(rd)
        skill_scores.append(skill_score)
        ax.text(idx_i/8 - 1,1.05+skill_score,'%0.4f'%(skill_score))
        ax.scatter(idx_i/8 - 1,1.0+skill_score,s=100,facecolor='gray')
    print skill_scores

    ax.set_ylabel('Time (s)')
    ax.set_title('%s Sequence Time to Target'%(subject))
    ax.set_xticks(ind+width)
    ax.set_xticklabels(sequence_labels)
    ax.set_ylim([.25,1.25])
    #autolabel(rects)
    plt.show()


