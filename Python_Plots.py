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
    rho = {}
    rho_reverse = {}
    Pi = math.pi
    ln = math.log
    i = 0
    R_B_12_43 = {}
    R_B_14_23 = {}
    R_B_21_34 = {}
    R_B_41_32 = {}
    d = 350E-6 #thickness of sample in m
    d_sigma = 25E-6 #thickness uncertainty of sample in m
    x_plt = []
    y_plt = []
    y_err_plt = []


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
        
        #calculate transresistance

        for file in files:
            #print(file)
            if str(curr_numb[0]) in file:
                
                data = np.genfromtxt(folder + "\\" + file, delimiter="\t", skip_header = 2, dtype=None)
                R_B_12_43[i, 0] = data[3]/data[1]
                R_B_12_43[i, 1] = uncertpropmult(data[3], data[4], data[1], data[2])
                R_B_12_43[i, 2] = getmagneticfield(0, folder + "\\" + file)
                R_B_12_43[i, 3] = getmagneticfield(1, folder + "\\" + file)
                #print("R_12_43 =", R_B_12_43[i, 0], "+/-", R_B_12_43[i, 1],"Ohm" ,"@ B_12_43", R_B_12_43[i, 2], "+/-", R_B_12_43[i, 3], "kG")

            if str(curr_numb[1]) in file:
                
                data = np.genfromtxt(folder + "\\" + file, delimiter="\t", skip_header = 2, dtype=None)
                R_B_14_23[i, 0] = data[3]/data[1]
                R_B_14_23[i, 1] = uncertpropmult(data[3], data[4], data[1], data[2])
                R_B_14_23[i, 2] = getmagneticfield(0, folder + "\\" + file)
                R_B_14_23[i, 3] = getmagneticfield(1, folder + "\\" + file)
                #print("R_14_23 =", R_14_23, "+/-", R_14_23_sigma,"Ohm" ,"@ B_14_23 =", B_14_23, "+/-", B_14_23_sigma, "kG")


            if str(curr_numb[2]) in file:
                
                data = np.genfromtxt(folder + "\\" + file, delimiter="\t", skip_header = 2, dtype=None)
                R_B_21_34[i, 0] = data[3]/data[1]
                R_B_21_34[i, 1] = uncertpropmult(data[3], data[4], data[1], data[2])
                R_B_21_34[i, 2] = getmagneticfield(0, folder + "\\" + file)
                R_B_21_34[i, 3] = getmagneticfield(1, folder + "\\" + file)
                #print("R_21_34 =", R_21_34, "+/-", R_21_34_sigma,"Ohm" ,"@ B_21_34 =", B_21_34, "+/-", B_21_34_sigma, "kG")

            if str(curr_numb[3]) in file:
                
                data = np.genfromtxt(folder + "\\" + file, delimiter="\t", skip_header = 2, dtype=None)
                R_B_41_32[i, 0] = data[3]/data[1]
                R_B_41_32[i, 1] = uncertpropmult(data[3], data[4], data[1], data[2])
                R_B_41_32[i, 2] = getmagneticfield(0, folder + "\\" + file)
                R_B_41_32[i, 3] = getmagneticfield(1, folder + "\\" + file)
                #print("R_41_32 =", R_41_32, "+/-", R_41_32_sigma,"Ohm" ,"@ B_41_32 =", B_41_32 , "+/-", B_41_32_sigma , "kG")
        
        #print("\n")
        i+=1

    #TODO
    #MAGNETIC FIELD STRENGTH NEEDS TO BE AVERAGED IN THE RESISTIVITY CALCULATION

    #calculate resistivity

    for j in range(len(folders)):
        rho[j, 0] = Pi*d/(2*ln(2))*(R_B_12_43[j, 0]+R_B_14_23[j, 0])*function_f(R_B_12_43[j, 0]/R_B_14_23[j, 0])
        #rho_1_sigma = 
        #print("rho from R_12_43 and R_14_23 : rho =", rho[j, 0], "@ B_12_43", R_B_12_43[j, 2], "+/-", R_B_12_43[j, 3], "kG")
    
    #print("\n")

    for j in range(len(folders)):
        rho_reverse[j, 0] = Pi*d/(2*ln(2))*(R_B_21_34[j, 0]+R_B_41_32[j, 0])*function_f(R_B_21_34[j, 0]/R_B_41_32[j, 0])
        #rho_reverse_sigma = 
        #print("rho from R_21_34 and R_41_32 : rho =", rho_reverse[j, 0], "@ B_21_34", R_B_21_34[j, 2], "+/-", R_B_21_34[j, 3], "kG")

    for l in range(len(folders)):
        x_plt.append(0.25*(float(R_B_21_34[l, 2]) + float(R_B_41_32[l, 2]) + float(R_B_12_43[l, 2]) + float(R_B_14_23[l, 2])))
        y_plt.append(0.5*(rho[l, 0] + rho_reverse[l, 0]))
   
    

    #plt.errorbar(x, y, yerr=y_errors, fmt='x', capsize=5) 
    plt.errorbar(x_plt, y_plt, fmt='x', capsize=5) 
    plt.grid()
    plt.title("Resistivity in dependence of the magnetic field strength")

    plt.show()
    plt.clf()

    return

def main():
        
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

if __name__ == "__main__" :
    main()