import numpy


fp_features = open("features.txt", "w")
fp_train = open("training.txt", "w")
fp_test = open("test.txt", "w")

fp_features.write("instance, x\n")

i = 0
for x in numpy.arange(0, 1.01 ,0.01):
    inst = "x1_%.2f" %(x)
    fp_features.write("%s, %.2f\n" %(inst, x))
    
    if i % 2 == 0:
        fp_train.write("%s\n" %(inst))
    else:
        fp_test.write("%s\n" %(inst))
        
    i += 1
        
fp_features.flush()
fp_features.close()
fp_train.flush()
fp_train.close()
fp_test.flush()
fp_test.close()            
    
    