###Instruction: This script takes all the pssm files in the folder "pssm_files" and will return the accurancy as well as standard deviation after cross validation. The window size is set to be in a loop from 7 to 23. It has to be put in the folder "pssm_files" as well, which also contains the raw dataset "0309homo70_with_labels.txt". The sequence identity cut-off is 0.7. After running this script no mediate file will be created.





feature_string='VYWTSPFMKLIHGEQCDNRA'
raw_list=[0]*38
raw_list.insert(19,1)

dict_of_labels={'H':0,'S':1,'C':2,'NO':3}
dict_of_aa=dict()

i=0
for amino_acids in feature_string:
    dict_of_aa[amino_acids]=raw_list[i:i+20]
    i=i+1


for window_size in range(7,25,2):
	half_ws=int(window_size/2)
	ws_list=[]
	key_pssm=0
	pssm_list=[]
	import glob
	for filename in glob.glob('*.pssm'):
		pssm_list.append(filename.strip('.fa.pssm'))
		f=open(filename)
		length_seq_list=[]
		i=half_ws
		for k in range(0,i,1):
			length_seq_list.append([0.0]*20)
		crazy=0
		for line in f:
			crazy=crazy+1
			if crazy>=4:
				if len(line)>=50:
					frequency=line.split()	
					frequency_list=frequency[22:42]
					freq_num_list=[(float (x))/100 for x in frequency_list]
					length_seq_list.append(freq_num_list)
					i=i+1
		for m in range(i,i+half_ws,1):
			length_seq_list.append([0.0]*20)
		f.close()

 ##add window size
		for why in range(half_ws,len(length_seq_list)-half_ws,1):
			list_in_list=[]
			for i in range(why-half_ws,why+half_ws+1,1):
				list_in_list.extend(length_seq_list[i])
			ws_list.append(list_in_list)
		f.close()
	

 ##create label list
	dict_homo70={}
	label_list=[]
	k=open('0309homo70_with_labels.txt')
	n=0
	for line in k:
		if n%3==0:
			a=line
			dict_homo70[a]=0
		if n%3==2:
			dict_homo70[a]=line
		n=n+1

	for genename in pssm_list:
		for key,value in dict_homo70.items():
			if key.strip()==genename.strip():
				label_list.extend(list(value.strip()))
	k.close()

	label_list2=[]
	for label in label_list:
		for key,value in dict_of_labels.items():
			if key==label:
				label_list2.append(value)

	from sklearn import svm
	from sklearn.model_selection import cross_val_score
	clf = svm.LinearSVC(C=1.0)
	scores = cross_val_score(clf, ws_list, label_list2, cv=5, verbose=40, n_jobs=-1)
	print(scores)
	print("Accuracy: %0.6f (+/- %0.6f)" % (scores.mean(), scores.std() * 2))
	print ('number %s has been done'%window_size)
