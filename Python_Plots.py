import numpy as np
import math
import re
import matplotlib.pyplot as plt
import os
import sympy as sym

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
        result = magneticfield[i].replace('KGs', '')
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
        print("Derivative has to be with respect to R_12_43 or R_14_23!")
        result = None

    return result

def readinvalues_Hallcoefficient(dirName):

    data = {}

    i = 0

    for root, directory, files in os.walk(dirName):
        for file in files:
           if '.txt' in file:
                files.append(os.path.join(root, file))
     

    for file in files:
        #print(file)

        if str(i) in file:
            print("Reading in", file)
            readindata = np.genfromtxt(dirName + "\\" + file, delimiter="\t", skip_header = 2, dtype=None)

            data[i, 0] = float(getmagneticfield(0, dirName + "\\" + file)) #B
            data[i, 1] = float(getmagneticfield(1, dirName + "\\" + file)) #sigma_B

            for j in range(2,29):
                data[i, j] = float(readindata[j - 2])
                
        else:

            print("No files found.")
     
        i += 1


    return data, i

def calc_resistivity ():

    #read in voltages, currents and magnetic field strengths

    d_sample = 350E-6 #thickness of sample in m
    d_sample_sigma = 25E-6 #thickness uncertainty of sample in m
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

        for file in files:
            print("Reading in ",file)
            if str(curr_numb[0]) in file:
                
                data = np.genfromtxt(folder + "\\" + file, delimiter="\t", skip_header = 2, dtype=None)
                R_B_12_43[i, 0] = data[3]/data[1]
                R_B_12_43[i, 1] = uncertpropdiv(data[3], data[4], data[1], data[2])
                R_B_12_43[i, 2] = getmagneticfield(0, folder + "\\" + file)
                R_B_12_43[i, 3] = getmagneticfield(1, folder + "\\" + file)
                #print("R_12_43 =", R_B_12_43[i, 0], "+/-", R_B_12_43[i, 1],"Ohm" ,"@ B_12_43" , R_B_12_43[i, 2], "+/-", R_B_12_43[i, 3], "kG")

            if str(curr_numb[1]) in file:
                
                data = np.genfromtxt(folder + "\\" + file, delimiter="\t", skip_header = 2, dtype=None)
                R_B_14_23[i, 0] = data[3]/data[1]
                R_B_14_23[i, 1] = uncertpropdiv(data[3], data[4], data[1], data[2])
                R_B_14_23[i, 2] = getmagneticfield(0, folder + "\\" + file)
                R_B_14_23[i, 3] = getmagneticfield(1, folder + "\\" + file)
                #print("R_14_23 =", R_B_14_23[i, 0], "+/-", R_B_14_23[i, 1],"Ohm" ,"@ B_14_23 =" , R_B_14_23[i, 2], "+/-", R_B_14_23[i, 3], "kG")


            if str(curr_numb[2]) in file:
                
                data = np.genfromtxt(folder + "\\" + file, delimiter="\t", skip_header = 2, dtype=None)
                R_B_21_34[i, 0] = data[3]/data[1]
                R_B_21_34[i, 1] = uncertpropdiv(data[3], data[4], data[1], data[2])
                R_B_21_34[i, 2] = getmagneticfield(0, folder + "\\" + file)
                R_B_21_34[i, 3] = getmagneticfield(1, folder + "\\" + file)
                #print("R_21_34 =", R_21_34, "+/-", R_21_34_sigma,"Ohm" ,"@ B_21_34 =", B_21_34, "+/-", B_21_34_sigma, "kG")

            if str(curr_numb[3]) in file:
                
                data = np.genfromtxt(folder + "\\" + file, delimiter="\t", skip_header = 2, dtype=None)
                R_B_41_32[i, 0] = data[3]/data[1]
                R_B_41_32[i, 1] = uncertpropdiv(data[3], data[4], data[1], data[2])
                R_B_41_32[i, 2] = getmagneticfield(0, folder + "\\" + file)
                R_B_41_32[i, 3] = getmagneticfield(1, folder + "\\" + file)
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

    plt.show()
    plt.clf()

    return

def calc_Hall_coefficient():

    #electric charge
    e = 16.022E-19
    d = 350E-6
    d_sigma = 25E-6

    dirName = r"C:\Users\Flo\Desktop\LabCourse\Hall Effect\Hall coefficient data 17.01.20\Flo"
    
    data = {}
    x_plt = []
    y_plt = []
    y_err_plt = []
    denom = {}
    num = {}

    data, n = readinvalues_Hallcoefficient(dirName)

       
    #prints the content of every file in the directory dirName
    #for i in range(n):
    #    print("File ", i + 1)
    #    for j in range(0,29):
    #        print(data[i, j])
    
    #calculate Hall coefficient

    for l in range(n):
        x_plt.append(data[l, 0]*1E-1)
        denom[l, 0] = data[l, 3] * data[l, 0] * 1E-1
        denom[l, 1] = uncertpropmult(data[l, 0], data[l, 1], data[l, 3], data[l, 4])
        num[l, 0] = 0.25*(calc_V_H(data[l, 9], data[l, 11])+calc_V_H(data[l, 13], data[l, 15])+calc_V_H(data[l, 21], data[l, 23])+calc_V_H(data[l, 25], data[l, 27]))*d
        num[l, 1] = uncertpropmult(num[l, 0], uncertpropadd(0.5*uncertpropadd(data[l, 10], data[l, 12], 0, 0), 0.5*uncertpropadd(data[l, 14], data[l, 16], 0, 0), 0.5*uncertpropadd(data[l, 22], data[l, 24], 0, 0), 0.5*uncertpropadd(data[l, 26], data[l, 28], 0, 0)) , d, d_sigma)
        y_plt.append(num[l, 0]/denom[l, 0])
        y_err_plt.append(uncertpropdiv(num[l, 0], num[l, 1], denom[l, 0], denom[l, 1]))
        print("R_H =", y_plt[l], "+/-", y_err_plt[l])

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

def main():
        
    #calc_resistivity()
    calc_Hall_coefficient()

    
if __name__ == "__main__" :
    main()