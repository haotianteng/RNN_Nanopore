import h5py as h5, numpy as np, re, sys, urllib2, random
from nanoraw import correct_raw as cr
from tensorflow.contrib.learn.python.learn.datasets import base
from tensorflow.python.platform import gfile
URL = "http://data.genomicsresearch.org/internal/Nanopore/GN_003_R9_2D_20160928/reads/downloads/pass/"
DATA_PATH= "/home/haotianteng/UQ/BINF7000/Nanopore/GN_003_R9_2D_20160928/"
GMAP_DIRECTORY = '/home/haotianteng/UQ/BINF7000/Nanopore/graphmap/bin/Linux-x64/graphmap'
REFERENCE_GENOME = 'pacbio_ref.fa'
class DataReader:
	# Three numpy arrays that get output to the user
    train_event = None
    valid_event = None
    test_event = None
    
    train_base = None
    valid_base = None
    test_base = None

    def __init__(self, train_weight=0.3, test_weight=0.3, valid_weight=0.4 ,windowsize = 100,reads_num= 300):
        
        total_weight = float(train_weight+test_weight+valid_weight)
        train_weight = train_weight/total_weight
        test_weight = test_weight/total_weight
        valid_weight = valid_weight/total_weight
        file_list = self.download_reads(reads_num)
        data_size = len(file_list)
        train_size = int(data_size * train_weight)
        test_size = int(data_size*test_weight)
        valid_size = data_size - train_size - test_size
        [train_event,train_label] = self.prepare_data_for_training(file_list[0:train_size],windowsize)
        [test_event,test_label] = self.prepare_data_for_training(file_list[train_size:train_size+valid_size],windowsize)
        [valid_event,valid_label] = self.prepare_data_for_training(file_list[train_size+valid_size:],windowsize)
        self.train_event = np.array(train_event)
        self.test_event = np.array(test_event)
        self.valid_event = np.array(valid_event)
        self.train_base = np.array(train_label)
        self.test_base =  np.array(test_label)
        self.valid_base = np.array(valid_label)
        
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
        
    def prepare_data_for_training(self,fileList,chopsize,save_to_file = None):
        event = []
        label = []
        chopsize = int(chopsize)
        for output_file_path in fileList:
            try:
                output = self.mapping_extract(output_file_path,DATA_PATH+REFERENCE_GENOME)
                if output is not None:
                    print "Mapping completed, chop into short reads and add into the dataset."
                    event+=self.window_slide(output[0],chopsize)
                    #label+=[i[chopsize/2] for i in self.window_slide(output[1],chopsize)]
                    label+=self.window_slide(output[1],chopsize)
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
    def window_slide(self,L,win):
        assert len(L)>=win
        out = [L[i:(i+win)] for i in range(len(L)-win+1)]
        return out               
        

#if __name__ == '__main__':
#	dataRead = DataReader(5, 50, 25 ,25)