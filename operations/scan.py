
def load_scan(inputfile):
    import numpy as np
    import re
    #read&match
    with open(inputfile,'r') as file_object:
        line = file_object.readlines()
            
    r_array = []
    a_array = []
    d_array = []
        
    for i in range(len(line)):
        if len(line[i].split()) == 5:
            temp_r = re.findall('R\s\d+\s\d+\s\d+\.\d+\s\d+',line[i])
            if len(temp_r) > 0 :
                temp_r = temp_r[0].split()
                for j in range(len(temp_r)-1):
                    r_array = np.append(r_array,float(temp_r[j+1]))
        
        elif len(line[i].split()) == 6:            
            temp_a = re.findall('A\s\d+\s\d+\s\d+\s\d+\.\d+\s\d+',line[i])
            if len(temp_a) > 0 :
                temp_a = temp_a[0].split()
                for j in range(len(temp_a)-1):
                    a_array = np.append(a_array,float(temp_a[j+1]))
        
        elif len(line[i].split()) == 7:  
            temp_d = re.findall('D\s\d+\s\d+\s\d+\s\d+\s\d+\.\d+\s\d+',line[i])
            if len(temp_d) > 0 :
                temp_d = temp_d[0].split()
                for j in range(len(temp_d)-1):
                    d_array = np.append(d_array,float(temp_d[j+1]))            
        
    if (len(r_array) > 0) :
        r_array = r_array.reshape(len(r_array)//4,4)
    if (len(a_array) > 0) :
        a_array = a_array.reshape(len(a_array)//5,5)
    if (len(d_array) > 0) :
        d_array = d_array.reshape(len(d_array)//6,6)

    return r_array, a_array, d_array