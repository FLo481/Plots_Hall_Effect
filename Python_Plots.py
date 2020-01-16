import numpy as np
import math
import re
import matplotlib.pyplot as plt
import os

def uncertpropmult(num, num_sigma, denom, denom_sigma):
    value = num/denom
    result = math.sqrt((num_sigma/num)**2+(denom_sigma/denom)**2)*value
    return result

def function_f(x):
    result = 1/math.cosh(math.log(x)/2.403)
    return result

def getmagneticfield(i, pathtoFile):

    if i < 2:

        pattern = '....KGs'

        with open (pathtoFile, 'rt') as myfile:
            for myline in myfile:   
                break
        magneticfield = re.findall(pattern, myline)
        result = magneticfield[i].replace('KGs', '')
    else:
        result = 0
        print("getmagneticfield(i) can only take values 0 and 1")
   
    return result

def calc_resistivity ():

    #read in voltages, currents and magnetic field strengths

    curr_numb = [12, 14, 21, 41]
    folders = []
    files = []
    values = [[]]


    dirName = r"C:\Users\Flo\Desktop\LabCourse\Hall Effect\resistivity data 15.01.20\Flo"

    for root, directory, folder in os.walk(dirName):
        for folder in directory:
            folders.append(os.path.join(root, folder))

    for folder in folders:
        #print(folder)
        for root, directory, files in os.walk(folder):
            for file in files:
                if '.txt' in file:
                    files.append(os.path.join(root, file))

        for file in files:
            #print(file)
            if str(curr_numb[0]) in file:
                #print(file)
                data = np.genfromtxt(folder + "\\" + file, delimiter="\t", skip_header = 2, dtype=None)
                #print(data)
                R_12_43 = data[3]/data[1]
                R_12_43_sigma = uncertpropmult(data[3], data[4], data[1], data[2])
                print("R_12_43 =", R_12_43, "+/-", R_12_43_sigma,"Ohm" ,"@ ", getmagneticfield(0, folder + "\\" + file), "kG")

   


    #calculate transresistance

    #R_12_43 = V_43/I
    #R_12_43_sigma = uncertpropmult(V_43, V_43_sigma, I, I_sigma)

    #R_14_23 = V_23/I
    #R_14_23_sigma = uncertpropmult(V_23, V_23_sigma, I, I_sigma)

    #R_21_34 = V_34/I
    #R_21_34_sigma = uncertpropmult(V_34, V_34_sigma, I, I_sigma)

    #R_41_32 = V_32/I
    #R_41_32_sigma = uncertpropmult(V_32, V_32_sigma, I, I_sigma)

    #calculate resistivity

    #rho = Pi*d/(2*math.log(2))*(R_12_43+R_14_23)*function_f(R_12_43/R_14_23)
    #print(rho)
    #rho_sigma = 
    #print(rho_sigma)

    return

d = 350E-6 #thickness of sample in m
d_sigma = 25E-6
Pi = math.pi
#path = r"C:\Users\Flo\Desktop\LabCourse\Hall Effect\resistivity data 15.01.20\Flo\0.8 KG\I_12 and V_43.dat"

calc_resistivity()

#plot graph(s)


#x = [1, 2, 3, 5]   
#y = [4, 3, 6, 2]   
#y_errors = [0.1, 0.5, .2, 0.25] 


#plt.errorbar(x, y, yerr=y_errors, fmt='x', capsize=5) 
#plt.grid()
#plt.title("Dat is eine Gülle")

#plt.show()
#plt.clf()