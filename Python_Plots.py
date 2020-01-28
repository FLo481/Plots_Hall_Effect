import numpy as np
import math
import re
import matplotlib.pyplot as plt
import os
import sympy as sym
import csv

def uncertpropdiv(num, num_sigma, denom, denom_sigma):
    value = num/denom
    result = math.sqrt((num_sigma/num)**2+(denom_sigma/denom)**2)*value
    return result

def uncertpropmult(a, a_sigma, b, b_sigma):
    value = a*b
    result = math.sqrt((a_sigma/a)**2+(b_sigma/b)**2)*value
    return result

def uncertpropadd(a_sigma, b_sigma, c_sigma, d_sigma):

    return math.sqrt(a_sigma**2 + b_sigma**2 + c_sigma**2 + d_sigma**2)

def function_f(x):
    result = 1/math.cosh(math.log(x)/2.403)
    return result

def calc_V_H(V_antiparallel, V_parallel):

    return 0.5*(V_antiparallel - V_parallel)

def getmagneticfield(i, pathtoFile):

    if i < 2:

        pattern = '....KGs'

        with open (pathtoFile, 'rt') as myfile:
            for myline in myfile:   
                break
        magneticfield = re.findall(pattern, myline)
        result = float(magneticfield[i].replace('KGs', ''))
    else:
        result = 0
        print("getmagneticfield(i) can only take values 0 and 1")
   
    return result

def rho_for_deriv(R_12_43, R_14_23, d):

    return math.pi*d/(2*math.log(2))*(R_12_43 + R_14_23)*1/(sym.cosh(sym.log(R_12_43/R_14_23)/2.403))

def deriv_of_rho(derive_by, point1, point2, point3):

    R_12_43, R_14_23, d = sym.symbols('R_12_43 R_14_23 d')    

    if derive_by == R_12_43 or derive_by == R_14_23 or derive_by == d:
        del_rho = sym.diff(rho_for_deriv(R_12_43, R_14_23, d), derive_by)
        result = del_rho.evalf(subs={R_12_43: point1, R_14_23: point2, d: point3})

    else:
        print("Derivative has to be with respect to R_12_43 or R_14_23 or d!")
        result = None

    return result

def readinvalues_Hall_coeff(dirName):

    data = {}
    i = 0 #counts the files
    j = 0 #for skipping the first two lines
    line = 0 #counts the lines in a file
    pattern = '....KGs' #pattern for magnetic field and magnetic field error search

   

    for root, directory, files in os.walk(dirName):
         for file in files:
             if '.txt' in file:
                 files.append(os.path.join(root, file))
     

    for file in files:
        with open(dirName + "\\" + file) as f:
            reader = csv.reader(f, delimiter = '\t')
            for row in reader:
                if j > 1:
                    if float(row[0]) != 0.0:
                        data[line, 0] = magneticfield
                        data[line, 1] = magneticfield_error
                        for k in range (2, 29):
                            data[line, k] = float(row[k - 2])
                            #print(data[line, k])
                        line += 1
  
                elif j == 0:
                    s1 = re.findall(pattern, row[0])
                    magneticfield = float(s1[0].replace('KGs', ''))
                    magneticfield_error = float(s1[1].replace('KGs', ''))
                    #print(data[line, 0], "+/-", data[line, 1])
                
                
                j += 1
                
            #    print("New line")
            #print("New file")
        j = 0    
        i += 1

    #replace all 0.0 error with 0.001
    for o in range(line - 1):
        for u in range(3,29):
            if u % 2 == 0:
                if data[o, u] == 0.0:
                    data[o, u] = 0.001

    #for o in range(6):
    #    print("Line %d" % o)
    #    for u in range(29):
    #        print(data[o , u])
         
 
    return data, i, line

def calc_resistivity ():

    #read in voltages, currents and magnetic field strengths

    d_sample = 350E-6 #thickness of sample in m
    d_sample_sigma = 25E-6 #thickness uncertainty of sample in m
    s = math.sqrt(10) #to get the standard error instead of the standard deviation
    curr_numb = [12, 14, 21, 41]
    folders = []
    files = []
    Pi = math.pi
    ln = math.log
    i = 0
    R_B_12_43 = {}
    R_B_14_23 = {}
    R_B_21_34 = {}
    R_B_41_32 = {}
    rho = {}
    rho_reverse = {}
    x_plt = []
    y_plt = []
    y_err_plt = []
    R_12_43, R_14_23, d = sym.symbols('R_12_43 R_14_23 d')

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
        print("Reading in from ", folder)
        for file in files:
            print("Reading in ",file)
            if str(curr_numb[0]) in file:
                
                data = np.genfromtxt(folder + "\\" + file, delimiter="\t", skip_header = 2, dtype=None)
                R_B_12_43[i, 0] = data[3]/data[1]
                R_B_12_43[i, 1] = uncertpropdiv(data[3], data[4]/s, data[1], data[2]/s)
                R_B_12_43[i, 2] = getmagneticfield(0, folder + "\\" + file)
                R_B_12_43[i, 3] = getmagneticfield(1, folder + "\\" + file)/s
                #print("R_12_43 =", R_B_12_43[i, 0], "+/-", R_B_12_43[i, 1],"Ohm" ,"@ B_12_43" , R_B_12_43[i, 2], "+/-", R_B_12_43[i, 3], "kG")

            if str(curr_numb[1]) in file:
                
                data = np.genfromtxt(folder + "\\" + file, delimiter="\t", skip_header = 2, dtype=None)
                R_B_14_23[i, 0] = data[3]/data[1]
                R_B_14_23[i, 1] = uncertpropdiv(data[3], data[4]/s, data[1], data[2]/s)
                R_B_14_23[i, 2] = getmagneticfield(0, folder + "\\" + file)
                R_B_14_23[i, 3] = getmagneticfield(1, folder + "\\" + file)/s
                #print("R_14_23 =", R_B_14_23[i, 0], "+/-", R_B_14_23[i, 1],"Ohm" ,"@ B_14_23 =" , R_B_14_23[i, 2], "+/-", R_B_14_23[i, 3], "kG")


            if str(curr_numb[2]) in file:
                
                data = np.genfromtxt(folder + "\\" + file, delimiter="\t", skip_header = 2, dtype=None)
                R_B_21_34[i, 0] = data[3]/data[1]
                R_B_21_34[i, 1] = uncertpropdiv(data[3], data[4]/s, data[1], data[2]/s)
                R_B_21_34[i, 2] = getmagneticfield(0, folder + "\\" + file)
                R_B_21_34[i, 3] = getmagneticfield(1, folder + "\\" + file)/s
                #print("R_21_34 =", R_21_34, "+/-", R_21_34_sigma,"Ohm" ,"@ B_21_34 =", B_21_34, "+/-", B_21_34_sigma, "kG")

            if str(curr_numb[3]) in file:
                
                data = np.genfromtxt(folder + "\\" + file, delimiter="\t", skip_header = 2, dtype=None)
                R_B_41_32[i, 0] = data[3]/data[1]
                R_B_41_32[i, 1] = uncertpropdiv(data[3], data[4]/s, data[1], data[2]/s)
                R_B_41_32[i, 2] = getmagneticfield(0, folder + "\\" + file)
                R_B_41_32[i, 3] = getmagneticfield(1, folder + "\\" + file)/s
                #print("R_41_32 =", R_41_32, "+/-", R_41_32_sigma,"Ohm" ,"@ B_41_32 =", B_41_32 , "+/-", B_41_32_sigma , "kG")
        
        #print("\n")
        i += 1

    #calculate resistivity

    for j in range(len(folders)):
        rho[j, 0] = Pi*d_sample/(2*ln(2))*(R_B_12_43[j, 0]+R_B_14_23[j, 0])*function_f(R_B_12_43[j, 0]/R_B_14_23[j, 0])
        rho[j, 1] = math.sqrt((deriv_of_rho(R_12_43, R_B_12_43[j, 0], R_B_14_23[j, 0], d_sample)*R_B_12_43[j, 1])**2 + (deriv_of_rho(R_14_23, R_B_12_43[j, 0], R_B_14_23[j, 0], d_sample)*R_B_14_23[j, 1])**2 + (deriv_of_rho(d, R_B_12_43[j, 0], R_B_14_23[j, 0], d_sample)*d_sample_sigma)**2)
        #print("rho from R_12_43 and R_14_23 : rho =", rho[j, 0]*1E6,"+/-", rho[j, 1]*1E6, "µmOhm", "@ B_12_43", R_B_12_43[j, 2], "+/-", R_B_12_43[j, 3], "kG")
    
    #print("\n")

    for j in range(len(folders)):
        rho_reverse[j, 0] = Pi*d_sample/(2*ln(2))*(R_B_21_34[j, 0]+R_B_41_32[j, 0])*function_f(R_B_21_34[j, 0]/R_B_41_32[j, 0])
        rho_reverse[j, 1] = math.sqrt((deriv_of_rho(R_12_43, R_B_21_34[j, 0], R_B_41_32[j, 0], d_sample)*R_B_21_34[j, 1])**2 + (deriv_of_rho(R_14_23, R_B_21_34[j, 0], R_B_41_32[j, 0], d_sample)*R_B_41_32[j, 1])**2 + (deriv_of_rho(R_14_23, R_B_21_34[j, 0], R_B_41_32[j, 0], d_sample)*d_sample_sigma)**2)
        #print("rho from R_21_34 and R_41_32 : rho =", rho_reverse[j, 0]*1E6, "+/-", rho_reverse[j, 1]*1E6, "µmOhm", "@ B_21_34", R_B_21_34[j, 2], "+/-", R_B_21_34[j, 3], "kG")

    for l in range(len(folders)):
        x_plt.append(0.25*1E-1*(float(R_B_21_34[l, 2]) + float(R_B_41_32[l, 2]) + float(R_B_12_43[l, 2]) + float(R_B_14_23[l, 2])))
        y_plt.append(0.5*1E6*(rho[l, 0] + rho_reverse[l, 0]))
        #y_err_plt.append(0.5*1E6*math.sqrt(rho[l, 1]**2 + rho_reverse[l, 1]**2))
        y_err_plt.append(0.5*1E6*uncertpropadd(rho[l, 1], rho_reverse[l, 1], 0, 0))
   
    
    plt.errorbar(x_plt, y_plt, y_err_plt, fmt='x', capsize=5) 
    plt.grid()
    plt.title("Resistivity in dependence of the magnetic field strength")
    plt.xlabel("B[T]")
    plt.ylabel(r"$\rho$ $[\mu m$ $\Omega]$")

    plt.show()
    plt.clf()

    return

def calc_Hall_coefficient():

    #electric charge
    e = 16.022E-19
    d = 350E-6
    d_sigma = 25E-6
    s = math.sqrt(10) #to get the standard error instead of the standard deviation

    dirName = r"C:\Users\Flo\Desktop\LabCourse\Hall Effect\Hall coefficient data 17.01.20\Flo"
    
    data = {}
    x_plt = []
    y_plt = []
    y_err_plt = []
    denom = {}
    num = {}


    data, n, j = readinvalues_Hall_coeff(dirName)

       
    #prints the content of every file in the directory dirName
    #for i in range(n):
    #    print("File ", i + 1)
    #    for j in range(0,29):
    #        print(data[i, j])
    
    #calculate Hall coefficient

    for l in range(n):
        x_plt.append(data[l, 0]*1E-1)
        denom[l, 0] = data[l, 3] * data[l, 0] * 1E-1
        denom[l, 1] = uncertpropmult(data[l, 0], data[l, 1]/s, data[l, 3], data[l, 4]/s)
        num[l, 0] = 0.25*(calc_V_H(data[l, 9], data[l, 11])+calc_V_H(data[l, 13], data[l, 15])+calc_V_H(data[l, 21], data[l, 23])+calc_V_H(data[l, 25], data[l, 27]))*d
        num[l, 1] = uncertpropmult(num[l, 0], uncertpropadd(0.5*uncertpropadd(data[l, 10]/s, data[l, 12]/s, 0, 0), 0.5*uncertpropadd(data[l, 14]/s, data[l, 16]/s, 0, 0), 0.5*uncertpropadd(data[l, 22]/s, data[l, 24]/s, 0, 0), 0.5*uncertpropadd(data[l, 26]/s, data[l, 28]/s, 0, 0)) , d, d_sigma)
        y_plt.append(num[l, 0]/denom[l, 0])
        y_err_plt.append(uncertpropdiv(num[l, 0], num[l, 1], denom[l, 0], denom[l, 1]))
        #print("R_H =", y_plt[l], "+/-", y_err_plt[l])

    #removes the first value, because it has a huge error bar
    del x_plt[0]
    del y_plt[0]
    del y_err_plt[0]
    plt.errorbar(x_plt, y_plt, y_err_plt, fmt='x', capsize=5) 
    plt.grid()
    plt.title("Hall coefficient in dependence of the magnetic field")
    plt.xlabel("B[T]")
    plt.ylabel("R_H[?]")

    plt.show()
    plt.clf()         
   
    return

def calc_tempdep_Hall_coefficient():

    dirName = r"C:\Users\Flo\Desktop\LabCourse\Hall Effect\Temperature dependent data 22.01.20\Konny"

    data = {}
    x_plt = []
    y_plt = []
    y_err_plt = []

    data, i, j = readinvalues_Hall_coeff(dirName)

    print("Number of lines", j)



    for l in range (j):
        x_plt.append(data[l, 2])
        


    return

def main():
        
    #calc_resistivity()
    calc_Hall_coefficient()
    #calc_tempdep_Hall_coefficient()

    
if __name__ == "__main__" :
    main()