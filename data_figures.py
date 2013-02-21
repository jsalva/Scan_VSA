import numpy as np
import matplotlib.pyplot as plt
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
    sequence_data_file = data[subject][0][0][0]
    sequence_data = np.genfromtxt(sequence_data_file,skip_header=True)
    sequence_data = sequence_data[sequence_data[:,1]!=0]
    conditions = sequence_data[:,1]

    sd = [conditions == condition['SD']]
    sd_time_to_target.extend(sequence_data[:,4][sd])
    sd_endpoint_error.extend(sequence_data[:,5][sd])
    sd_accuracy.extend(sequence_data[:,10][sd])

    rd = [conditions == condition['RD']]
    rd_time_to_target.extend(sequence_data[:,4][rd])
    rd_endpoint_error.extend(sequence_data[:,5][rd])
    rd_accuracy.extend(sequence_data[:,10][rd])


ri = None
ri_time_to_target=[]
ri_endpoint_error=[]
ri_accuracy=[]

rd = None
rd_time_to_target=[]
rd_endpoint_error=[]
rd_accuracy=[]

for subject in subject_list:
    perturb_data_file = data[subject][1][0][0]
    perturb_data = np.genfromtxt(perturb_data_file,skip_header=True)
    perturb_data = perturb_data[perturb_data[:,1]!=0]
    conditions = perturb_data[:,1]

    ri = [conditions == condition['RI']]
    ri_time_to_target.extend(sequence_data[:,4][ri])
    ri_endpoint_error.extend(sequence_data[:,5][ri])
    ri_accuracy.extend(sequence_data[:,10][ri])

    rd = [conditions == condition['RD']]
    rd_time_to_target.extend(sequence_data[:,4][rd])
    rd_endpoint_error.extend(sequence_data[:,5][rd])
    rd_accuracy.extend(sequence_data[:,10][rd])

    
    #repeat for perturb
