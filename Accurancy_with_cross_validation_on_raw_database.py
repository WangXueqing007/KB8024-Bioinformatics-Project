###Instruction: This script takes the whole dataset and will return the accurancy as well as standard deviation in window size 9-27. Cross validation has also been done. It has to be put in the same folder as file "0309homo70_with_labels.txt", which is generated by doing homology reduction on the raw database with the helo of CD HIT. The sequence identity cut-off is 0.7. After running this script, two mediate files will be created on the same folder.

f=open('0309homo70_with_labels.txt','r')

###CREATE A DICTIONARY TO MATCH EACH AMINO ACID WITH UNIQUE BINARY CODE
n=0
feature=set()
for line in f:
    if n%3==1:
        for i in line.strip():
            feature.add(i)
    n=n+1

f.close()

feature_list=list(feature)
raw_list=[0]*38
raw_list.insert(19,1)


dict_of_aa=dict()

i=0
for amino_acids in feature_list:
    dict_of_aa[amino_acids]=raw_list[i:i+20]
    i=i+1
dict_of_aa['B']=[0]*20


dict_of_labels={'H':0,'S':1,'C':2,'NO':3}

###FOR DIFFERENT WINDOW SIZE, EXTRACT DATA FROM THE 3LINE RAW FILE
for window_size in range(9,29,2):
	features=list()
	labels=list()
	f=open('0309homo70_with_labels.txt','r')
	m=0
	for line in f:
		if m%3==0:
			features.extend(['B']*(int(window_size/2)))
			labels.extend(['NO']*(int(window_size/2)))
		if m%3==1:
			for i in line.strip():
				features.append(i)
		if m%3==2:
			for i in line.strip():
				labels.append(i)
		m=m+1
	features.extend(['B']*(int(window_size/2)))
	labels.extend(['NO']*(int(window_size/2)))

	f.close()

	f=open('0309homo70_with_labels.txt')
	o=open('mediate1.txt','w')
	for i in range(int(window_size/2),int(len(labels)-window_size/2)):
		o.write(labels[i]+' ')
		for window_member in range(i-int(window_size/2),i+int(window_size/2)+1):
			o.write(features[window_member])
		o.write('\n')
	f.close()
	o.close()

	f=open('mediate1.txt')
	p=open('mediate2.txt','w')
	for line in f:
		list_every_line=line.strip().split(' ')
		for key,value in dict_of_labels.items():
			if list_every_line[0]==key:
				p.write(str(value)+' ')
		list_every_line_aa=list(list_every_line[1])
		for num in range(0,window_size-1,1):
			for key,value in dict_of_aa.items():
				if list_every_line_aa[num]==key:
					for num_binary in range(0,20,1):
						p.write(str(value[num_binary]))
		for key, value in dict_of_aa.items():
			if list_every_line_aa[window_size-1]==key:
				for num_last_binary in range(0,20,1):
					p.write(str(value[num_last_binary]))
		p.write('\n')
	f.close()
	p.close()

	f=open('mediate2.txt')
	p=open('input_sklearn.txt','w')
	for line in f:
		if line[0]!='3':
			p.write(line)
	f.close()
	p.close()

	p=open('input_sklearn.txt')
	list_of_labels=list()
	list_of_aa=list()
	for line in p:
		list_of_line=line.strip().split(' ')
		list_of_labels.append(int(list_of_line[0]))
		list_of_aa.append([float(x) for x in list_of_line[1]])
	p.close()

###CROSS VALIDATION(5 SETS)
	import numpy as np
	from sklearn import svm
	from sklearn import datasets
	from sklearn.model_selection import cross_val_score
	clf = svm.SVC(kernel='sigmoid', C=1.0, decision_function_shape='ovr')
	scores = cross_val_score(clf, list_of_aa, list_of_labels, cv=5, verbose=40, n_jobs=-1)
	print(scores)
	print("Accuracy: %0.6f (+/- %0.6f)" % (scores.mean(), scores.std() * 2))
	print ('number %s has been done'%window_size)

