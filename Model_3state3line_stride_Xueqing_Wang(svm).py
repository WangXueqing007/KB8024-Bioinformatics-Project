###Instruction: This script is my final model. To use it, one has to put this file into the "pssm_files" folder. It takes a query sequence and returns its secondary structure(In the form of a string consisting of three labels, H, S and C. Each label matches to a residue in the query sequence.) The window size is set to be 17 by default, but users can set it to other odd numbers.


X_test_string=input ('Welcome to use our 3-state(stride) secondary structure predictor. Please input your query sequence here(It should be a string of amino acid labels, but we only support 20 essential amino acids): ')
X_test_list1=list(X_test_string)


window_size=input ('Would you like to use 17 as the window size?(y/n) ')
if window_size=='y':
	window_size=17
if window_size=='n':
	window_size=int(input('Input an odd number as window size: '))

print ('Your sequence has been submitted. Processing...')

half_ws=int(window_size/2)
feature_string='VYWTSPFMKLIHGEQCDNRA'
raw_list=[0]*38
raw_list.insert(19,1)

dict_of_labels={'H':0,'S':1,'C':2,'NO':3}
dict_of_aa=dict()

X_test=[]
i=0
for amino_acids in feature_string:
    dict_of_aa[amino_acids]=raw_list[i:i+20]
    i=i+1

X_test_list2=[]
for k in range(0,half_ws,1):
	X_test_list2.append([0.0]*20)
for amino_acid in X_test_list1:
	for key,value in dict_of_aa.items():
		if key==amino_acid:
			X_test_list2.append(value)
for k in range(0,half_ws,1):
	X_test_list2.append([0.0]*20)

for why in range(half_ws,len(X_test_list2)-half_ws,1):
	list_in_list=[]
	for i in range(why-half_ws,why+half_ws+1,1):
		list_in_list.extend(X_test_list2[i])
	X_test.append(list_in_list)



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
X_train=ws_list	

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
y_train=label_list2

from sklearn import svm

clf = svm.LinearSVC(C=1.0)
clf.fit(X_train,y_train)
label_list1=clf.predict(X_test)
label_list2=[]
for element in label_list1:
	for key, value in dict_of_labels.items():
		if int(value)==int(element):
			label_list2.append(key)
label_string=''.join(label_list2)
print ('The secondary structure of your query sequence is as follows:')
print (label_string)
print ('Your work has been done. Thank you for using!')




