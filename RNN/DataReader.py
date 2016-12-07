import h5py as h5, numpy as np, re, sys, urllib2, random,os
from nanoraw import correct_raw as cr
from tensorflow.contrib.learn.python.learn.datasets import base
from tensorflow.python.platform import gfile
URL = "http://data.genomicsresearch.org/internal/Nanopore/GN_003_R9_2D_20160928/reads/downloads/pass/"
DATA_PATH= "/home/haotianteng/UQ/BINF7000/Nanopore/GN_003_R9_2D_20160928/"
GMAP_DIRECTORY = '/home/haotianteng/UQ/BINF7000/Nanopore/graphmap/bin/Linux-x64/graphmap'
REFERENCE_GENOME = 'pacbio_ref.fa'
WORK_PATH = "/home/haotianteng/UQ/BINF7000/Nanopore/data/"
class DataReader:
	# Three numpy arrays that get output to the user
    train_event = None
    valid_event = None
    test_event = None
    
    train_base = None
    valid_base = None
    test_base = None

    def __init__(self, train_weight=0.3, test_weight=0.3, valid_weight=0.4 ,windowsize = 100,reads_num= None,event_total = None,file_list=None,reference_genome =  REFERENCE_GENOME):
        
        total_weight = float(train_weight+test_weight+valid_weight)
        train_weight = train_weight/total_weight
        test_weight = test_weight/total_weight
        valid_weight = valid_weight/total_weight
        
        if file_list is None:
            file_list = self.download_reads(reads_num)
            data_size = len(file_list)
            train_size = int(data_size * train_weight)
            test_size = int(data_size*test_weight)
            valid_size = data_size - train_size - test_size
            [train_event,train_label] = self.fast5_to_event(file_list[0:train_size],windowsize,DATA_PATH+REFERENCE_GENOME)
            [test_event,test_label] = self.fast5_to_event(file_list[train_size:train_size+valid_size],windowsize,DATA_PATH+REFERENCE_GENOME)
            [valid_event,valid_label] = self.fast5_to_event(file_list[train_size+valid_size:],windowsize,DATA_PATH+REFERENCE_GENOME)
        else:
            event = []
            label = []
            if not isinstance(file_list,list):
                file_list = [file_list]
            for file_path in file_list:
                [file_name,filext] = os.path.splitext(file_path)
                if filext == '.fast5':
                    [tmp_event,tmp_label] = self.fast5_to_event(file_path,windowsize,DATA_PATH+REFERENCE_GENOME)
                elif filext  == '.dat':
                    [tmp_event,tmp_label] = self.dat_to_event(file_path,event_total)
                event.append(tmp_event)
                label.append(tmp_label)
            event = self._window_slide(event[0],windowsize)
            label = self._window_slide(label[0],windowsize)
            pack_for_shuffle = zip(event,label)
            random.shuffle(pack_for_shuffle)
            event = [pack[0] for pack in pack_for_shuffle]
            label = [pack[1] for pack in pack_for_shuffle]    
            if reads_num is not None:
                event = event[0:reads_num]
                label = label[0:reads_num]
            data_size = len(label)
            train_size = int(data_size * train_weight)
            test_size = int(data_size*test_weight)
            valid_size = data_size - train_size - test_size
            train_event = event[0:train_size]
            train_label = label[0:train_size]
            test_event = event[train_size:train_size+valid_size]
            test_label = label[train_size:train_size+valid_size]
            valid_event = event[train_size+valid_size:]
            valid_label = label[train_size+valid_size:]
        self.train_event = np.array(train_event)
        self.test_event = np.array(test_event)
        self.valid_event = np.array(valid_event)
        self.train_base = np.array(train_label)
        self.test_base =  np.array(test_label)
        self.valid_base = np.array(valid_label)
#    def read(self,file_list)

        
    def mapping_extract(self, item,reference) :
        print item
        # Dependant on the file path of the user's computer currently defaulting to your current working directory
        b_group = "Basecall_1D_000"
        c_group = "RawGenomeCorrected_000"
        # output the numpy array perhaps this should be returned and added to an array or object so we have a 3D Array
        # With each layer(file) containting an inner array with the numpy data
        try:
            array = cr.correct_raw_data(item, reference, GMAP_DIRECTORY, b_group,c_group, in_place = False)
            event_array = []
            base_array = []
            # Format this array, grabbing the first three
            print 'Array length: ',len(array)
            for elem in array:
                event_array.append([float(elem[0]), float(elem[1]), float(elem[3])])
                # Convert to a char, hopefully it takes up a bit less memory
                base_array.append(self._base2one_hot_vector(elem[4][0]))
                #
            print "Successfully read %s to array" %(item)
            return event_array, base_array;
        except RuntimeError:
            print "Runtime error on file: %s" %(item);
            return None;
    def _base2one_hot_vector(self,base):
        Alphabeta = ['A','G','C','T']
        alphabeta = ['a','g','c','t']
        label_array = [0,0,0,0]
        if ord(base)<97:
            label_array[Alphabeta.index(base)]=1
        else:
            label_array[alphabeta.index(base)]=1
        return label_array
    def download_reads(self,reads_num):
        #Load in the URL of the file
        # Establish a connection to the server
        print "Contact the server."
        page = urllib2.urlopen(URL).read()
        print "Get the file list."
        links = re.findall('<a href=(.*?)>.*?</a>', page)
        fast5Links = []
        fileList = []
        for link in links:
            #Determine if fast 5 file
            if "fast5" in link:
                fast5Links.append(link)
        print "Totally ",len(fast5Links),' fast5 files in the list.'
        # Download files
        if not gfile.Exists(DATA_PATH+REFERENCE_GENOME):          
            print 'Error, reference genome ', REFERENCE_GENOME, 'missing!'
            return           
        CurrentReadIndex=1
        for item in fast5Links[0: reads_num]:
            fileName = item[1:-1]
            try:
                output_file_path = base.maybe_download(fileName,DATA_PATH,URL+fileName)
                fileList.append(output_file_path)
            except:
                print "Error:",sys.exc_info()[0]
                continue
            print CurrentReadIndex,'/',reads_num
            CurrentReadIndex +=1
            ###Try to download the file
        return fileList
        
    def fast5_to_event(self,fileList,chopsize,reference,save_to_file = None):
        event = []
        label = []
        chopsize = int(chopsize)
        for output_file_path in fileList:
            try:
                output = self.mapping_extract(output_file_path,reference)
                if output is not None:
                    print "Mapping completed, chop into short reads and add into the dataset."
                    event+=self._window_slide(output[0],chopsize)
                    #label+=[i[chopsize/2] for i in self.window_slide(output[1],chopsize)]
                    label+=self._window_slide(output[1],chopsize)
            except IOError:
                print output_file_path,'file broken, skipped.'
            except:
                print "Unexpected Error on",output_file_path
        ###Stop extract if we have enough reads,
        pack_for_shuffle = zip(event,label)
        random.shuffle(pack_for_shuffle)
        event = [pack[0] for pack in pack_for_shuffle]
        label = [pack[1] for pack in pack_for_shuffle]
        return event,label
    def dat_to_event(self,file_path,event_num = None,save_to_file = WORK_PATH+'rnnp_', align_index = 2):
        ###################read the .dat file
        #align_index is the index of the 
        output_event = []
        output_label = []
        file_h = open(file_path,'r')
        self._mkdir_p(save_to_file)
        event_h = open(save_to_file+'event.txt','w+')
        label_h = open(save_to_file+'label.txt','w+')
        count = 0
        for line in file_h:
            if line[0]=="#":
            #Line begin with a # is not the data line.
                continue;
            else:
                data = line.split()
                label = data[0][align_index]
                event = data[1:4]
                label_h.write(label+'\n')
                event_h.write(' '.join(event)+'\n')
                count +=1
                event = [float(item) for item in event]
                label = self._base2one_hot_vector(label)         
                output_event.append(event)
                output_label.append(label)
                if count%1000 == 0:
                    if event_num is not None:
                        sys.stdout.write("%d/%d lines read.   \r" % (count,event_num) )
                    else:
                        sys.stdout.write("%d lines read.   \r" % (count) )
                    sys.stdout.flush()
            if event_num is not None:
                if count >= event_num:
                    break                    
        return output_event,output_label
    def _window_slide(self,L,win):
        assert len(L)>=win
        out = [L[i:(i+win)] for i in range(len(L)-win+1)]
        return out               
    def _mkdir_p(self,path):
        directory = os.path.dirname(path)
        if not os.path.exists(directory):
            os.makedirs(directory)
            print path," not exist, directory auto-create."
    def _Printer(self,data):
        """Print things to stdout on one line dynamically"""
        sys.stdout.write("\r\x1b[K"+data.__str__())
        sys.stdout.flush()

#if __name__ == '__main__':
#	dataRead = DataReader(5, 50, 25 ,25)