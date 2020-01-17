import numpy as np
import math
import re
import matplotlib.pyplot as plt
import os
import sympy as sym

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
            #print(file)
            if str(curr_numb[0]) in file:
                
                data = np.genfromtxt(folder + "\\" + file, delimiter="\t", skip_header = 2, dtype=None)
                R_B_12_43[i, 0] = data[3]/data[1]
                R_B_12_43[i, 1] = uncertpropmult(data[3], data[4], data[1], data[2])
                R_B_12_43[i, 2] = getmagneticfield(0, folder + "\\" + file)
                R_B_12_43[i, 3] = getmagneticfield(1, folder + "\\" + file)
                #print("R_12_43 =", R_B_12_43[i, 0], "+/-", R_B_12_43[i, 1],"Ohm" ,"@ B_12_43" , R_B_12_43[i, 2], "+/-", R_B_12_43[i, 3], "kG")

            if str(curr_numb[1]) in file:
                
                data = np.genfromtxt(folder + "\\" + file, delimiter="\t", skip_header = 2, dtype=None)
                R_B_14_23[i, 0] = data[3]/data[1]
                R_B_14_23[i, 1] = uncertpropmult(data[3], data[4], data[1], data[2])
                R_B_14_23[i, 2] = getmagneticfield(0, folder + "\\" + file)
                R_B_14_23[i, 3] = getmagneticfield(1, folder + "\\" + file)
                #print("R_14_23 =", R_B_14_23[i, 0], "+/-", R_B_14_23[i, 1],"Ohm" ,"@ B_14_23 =" , R_B_14_23[i, 2], "+/-", R_B_14_23[i, 3], "kG")


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
        y_err_plt.append(0.5*1E6*math.sqrt(rho[l, 1]**2 + rho_reverse[l, 1]**2))
   
    
    plt.errorbar(x_plt, y_plt, y_err_plt, fmt='x', capsize=5) 
    plt.grid()
    plt.title("Resistivity in dependence of the magnetic field strength")

    plt.show()
    plt.clf()

    return

def calc_Hall_coefficient():

    print("Hello World")
    
    return

def main():
        
    calc_resistivity()
    #calc_Hall_coefficient()

    
if __name__ == "__main__" :
    main()