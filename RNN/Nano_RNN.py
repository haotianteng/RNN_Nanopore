# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import tensorflow as tf
from DataReader import DataReader

#############################super parameter
TRANNING_READS = 10  #Total reads used
TRAIN_WEIGHT = 0.4  #Proportion of reads used to train
TEST_WEIGHT = 0.3   #Proportion of reads used to test
VALID_WEIGHT = 0.3  #Proportion of reads used to validate

EVENT_LENGTH = 15
HIDDEN_UNIT_NUM = 100#24
STEP_RATE = 2
BATCH_SIZE = 200
EPOCH = 5000
#############################

###############################Read the data
data = DataReader(TRAIN_WEIGHT,TEST_WEIGHT,VALID_WEIGHT,EVENT_LENGTH,TRANNING_READS)
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

##############################Execute the graph and train

loss_val_record = []
init = tf.initialize_all_variables()
sess = tf.Session()
sess.run(init)
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
    correct_prediction = tf.equal(tf.argmax(y,1), tf.argmax(y_,1))
    accuracy = tf.reduce_mean(tf.cast(correct_prediction, tf.float32))
    incorrect = sess.run(accuracy,feed_dict={x: test_event, y: test_label})
    print('Epoch {:2d} error {:3.1f}%'.format(i + 1, 100 * incorrect))
sess.close()





#lstm = tf.nn.rnn_cell.BasicLSTMCell(20);
#print DataReader.testingChars