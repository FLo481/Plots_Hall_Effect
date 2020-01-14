import numpy as np
import math
import matplotlib.pyplot as plt

d = 0.0001 #thickness of sample in m
Pi = math.pi

data = np.genfromtxt(r"C:\Users\Flo\Desktop\LabCourse\Hall Effect\Test data\test3.dat", skip_header=2, dtype=None, delimiter='\t')
#print(data)
temp = data[0]
current = data[1]
error_current = data[2]

#data[3] = V_43
#data[5] = V_23
#data[7] = V_24
#data[9] = V_24 RF
#data[11] = V_31
#data[13] = V_31 RF


#x = [1, 2, 3, 5]   
#y = [4, 3, 6, 2]   
#errors = [0.1, 0.5, .2, 0.25] 


#plt.errorbar(x, y, yerr=y_errors, fmt='x', capsize=5) 
#plt.grid()
#plt.title("Dat is eine Gülle")

#plt.show()
#plt.clf()