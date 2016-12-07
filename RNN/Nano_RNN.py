# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import tensorflow as tf
from DataReader import DataReader

#############################super parameter
TRANNING_READS = 1000000  #Total reads used, if None, use all.
TRANING_EVENTS_TOTAL = 2000000 # If is None, read all the events in the training data
TRAIN_WEIGHT = 0.4  #Proportion of reads used to train
TEST_WEIGHT = 0.4   #Proportion of reads used to test
VALID_WEIGHT = 0.2  #Proportion of reads used to validate

#Structure
EVENT_LENGTH = 40 #Length of each sentence
HIDDEN_UNIT_NUM = 24 #Length of the hidden state of each hidden layer
CellLayer = 3 #Number of the hidden layers

#Training
STEP_RATE = 0.1
BATCH_SIZE = 20
EPOCH = 500
#############################

#############################Path setting
Checkpoint_file = "Model_save.chk"
log_file = "Record.log"
#############################


###############################Read the data
data = DataReader(TRAIN_WEIGHT,TEST_WEIGHT,VALID_WEIGHT,EVENT_LENGTH,TRANNING_READS,event_total = TRANING_EVENTS_TOTAL,file_list = '/home/haotianteng/UQ/BINF7000/Nanopore/GN_003/event_pass.dat')
train_event = data.train_event
train_label = data.train_base
test_event = data.test_event
test_label = data.test_base
###############################
print "training size: ",len(train_event)
print "test size: ",len(test_event)


###############################Construct the graph

x = tf.placeholder(tf.float32,[None,EVENT_LENGTH,3])###each event has three properties, e.g. [mean,dev,length]
y = tf.placeholder(tf.float32,[None,EVENT_LENGTH,4])
cell = tf.nn.rnn_cell.LSTMCell(HIDDEN_UNIT_NUM,state_is_tuple=True)
cell = tf.nn.rnn_cell.MultiRNNCell([cell]*CellLayer,state_is_tuple=True)
lasth,state = tf.nn.dynamic_rnn(cell,x,dtype = tf.float32)
lasth = tf.transpose(lasth,[1,0,2])

###initialize the parameter
w = tf.Variable(tf.truncated_normal([EVENT_LENGTH,HIDDEN_UNIT_NUM,int(y.get_shape()[2])]))
bias = tf.Variable(tf.constant(0.1,shape = [y.get_shape()[1],y.get_shape()[2]]))


###evulation
print lasth.get_shape()
print w.get_shape()
y_ = tf.nn.softmax(tf.transpose(tf.batch_matmul(lasth,w,adj_x=False,adj_y=False),[1,0,2])+bias)
cross_entropy = -tf.reduce_mean(y*tf.log(y_))

###Choose the optimization algorithm
#optimizer = tf.train.AdamOptimizer()
train_step = tf.train.GradientDescentOptimizer(STEP_RATE).minimize(cross_entropy)
###using the stochastic optimizer as the DeepNano suggested  arXiv:1603.09195

##############################

###Saver of the model
saver = tf.train.Saver()
##############################


##############################Execute the graph and train

loss_val_record = []
test_accuracy_record = []
init = tf.initialize_all_variables()
sess = tf.Session()

###Check if there is any save point to be reload
try:
    ckpt = tf.train.import_meta_graph(Checkpoint_file+'.meta')
    # Restores from checkpoint
    ckpt.restore(sess,Checkpoint_file)
    print "Model loaded"
except IOError:
    print "No checkpoint file found"
    sess.run(init)
###

batches_num = int(len(train_event)/BATCH_SIZE)
for i in range(EPOCH):
    ptr = 0
    for j in range(batches_num):
        batch_x, batch_y = train_event[ptr:ptr+BATCH_SIZE], train_label[ptr:ptr+BATCH_SIZE]
        ptr+=BATCH_SIZE
        _,loss_val = sess.run([train_step,cross_entropy],feed_dict={x: batch_x, y: batch_y})
    loss_val_record.append(loss_val)
    print "Epoch - ",str(i),'loss = ',loss_val
    
    
############################Evulating the model
    correct_prediction = tf.equal(tf.argmax(y,2), tf.argmax(y_,2))
    accuracy = tf.reduce_mean(tf.cast(correct_prediction, tf.float32))
    correct = sess.run(accuracy,feed_dict={x: test_event, y: test_label})
    test_accuracy_record.append(correct)
    print('Epoch {:2d} error {:3.3f}%'.format(i + 1, 100 * correct))

############################Save the model and record
saver.save(sess,Checkpoint_file)
log_f = open(log_file,'w+')
log_f.write("#Loss_Value\n")
log_f.write(",".join(str(x) for x in loss_val_record))
log_f.write("\n")
log_f.write("#Accuracy_Test\n")
log_f.write(",".join(str(x) for x in test_accuracy_record))
log_f.write("\n")
log_f.close()
sess.close()





#lstm = tf.nn.rnn_cell.BasicLSTMCell(20);
#print DataReader.testingChars