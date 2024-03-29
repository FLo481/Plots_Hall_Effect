import numpy as np
import math
import re
import matplotlib.pyplot as plt
import os
import sympy as sym
import csv
import scipy.optimize
#from scipy.stats import chisquare

def fit_func1(x, a):

   return a*x

def fit_func2(x, a, b, c):

    return a*(x-b)**(3/2)+c

def fit_func3(x, a):

    return x/x*a

def fit_func4(x, a, b):
    e = 1.6E-19

    return -1/e*1/((a*x)**(3)*np.exp(-b/x))

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
    
    return 1/math.cosh(math.log(x)/2.403)

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
    s = np.sqrt(10)

   

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
                        data[line, 0] = magneticfield * 1E-1
                        data[line, 1] = magneticfield_error * 1E-1/s
                        for k in range (2, 29):
                            if k > 2:
                                data[line, k] = float(row[k - 2]) * 1E-3
                            else:
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

    #replace 0.0 with the device error


    for o in range(line):
        if data[o, 1] == 0.0:
            #magnetic field strength error replacement
            data[o, 1] == 0.01 * 1E-1/s
        for u in range(3,29):
            if u % 2 == 0:
                if data[o, u] == 0.0:
                    if k < 5:
                        #current error replacement
                        data[o, u] = (0.0001 * data[o, u - 1] + 0.00004 * 100) * math.sqrt(10) * 1E-3/s
                    elif k > 5:
                        #voltage error replacement
                        data[o, u] = (0.00003 * data[o, u - 1] + 0.000035 * 100) * 1E-3/s
                else:
                    data[o, u] = data[o, u]/s

   
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
    x_temp = []
    y_temp = []
    y_err_temp = []
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

    #expludes the last two points, due to a wrong magnetic field readout of the labview program
    n_excl = 2

    #calculate resistivity

    for j in range(len(folders)-n_excl):
        rho[j, 0] = Pi*d_sample/(2*ln(2))*(R_B_12_43[j, 0]+R_B_14_23[j, 0])*function_f(R_B_12_43[j, 0]/R_B_14_23[j, 0])
        #rho[j, 1] = math.sqrt((deriv_of_rho(R_12_43, R_B_12_43[j, 0], R_B_14_23[j, 0], d_sample)*R_B_12_43[j, 1])**2 + (deriv_of_rho(R_14_23, R_B_12_43[j, 0], R_B_14_23[j, 0], d_sample)*R_B_14_23[j, 1])**2 + (deriv_of_rho(d, R_B_12_43[j, 0], R_B_14_23[j, 0], d_sample)*d_sample_sigma)**2)
        rho[j, 1] = math.sqrt((deriv_of_rho(R_12_43, R_B_12_43[j, 0], R_B_14_23[j, 0], d_sample)*R_B_12_43[j, 1])**2 + (deriv_of_rho(R_14_23, R_B_12_43[j, 0], R_B_14_23[j, 0], d_sample)*R_B_14_23[j, 1])**2)
    
    #print("\n")

    for j in range(len(folders)-n_excl):
        rho_reverse[j, 0] = Pi*d_sample/(2*ln(2))*(R_B_21_34[j, 0]+R_B_41_32[j, 0])*function_f(R_B_21_34[j, 0]/R_B_41_32[j, 0])
        #rho_reverse[j, 1] = math.sqrt((deriv_of_rho(R_12_43, R_B_21_34[j, 0], R_B_41_32[j, 0], d_sample)*R_B_21_34[j, 1])**2 + (deriv_of_rho(R_14_23, R_B_21_34[j, 0], R_B_41_32[j, 0], d_sample)*R_B_41_32[j, 1])**2 + (deriv_of_rho(R_14_23, R_B_21_34[j, 0], R_B_41_32[j, 0], d_sample)*d_sample_sigma)**2)
        rho_reverse[j, 1] = math.sqrt((deriv_of_rho(R_12_43, R_B_21_34[j, 0], R_B_41_32[j, 0], d_sample)*R_B_21_34[j, 1])**2 + (deriv_of_rho(R_14_23, R_B_21_34[j, 0], R_B_41_32[j, 0], d_sample)*R_B_41_32[j, 1])**2)
        #print("rho from R_21_34 and R_41_32 : rho =", rho_reverse[j, 0]*1E6, "+/-", rho_reverse[j, 1]*1E6, "µmOhm", "@ B_21_34", R_B_21_34[j, 2], "+/-", R_B_21_34[j, 3], "kG")

    for l in range(len(folders)-n_excl):
        x_temp.append(0.25*(float(R_B_21_34[l, 2]) + float(R_B_41_32[l, 2]) + float(R_B_12_43[l, 2]) + float(R_B_14_23[l, 2])))
        y_temp.append(0.5*(rho[l, 0] + rho_reverse[l, 0]))
        #y_err_plt.append(0.5*1E6*math.sqrt(rho[l, 1]**2 + rho_reverse[l, 1]**2))
        y_err_temp.append(0.5*uncertpropadd(rho[l, 1], rho_reverse[l, 1], 0, 0))
   
    x_plt = np.empty(i-n_excl, dtype = float)
    y_plt = np.empty(i-n_excl, dtype = float)
    y_err_plt = np.empty(i-n_excl, dtype = float)

    x_plt[:] = x_temp
    y_plt[:] = y_temp
    y_err_plt[:] = y_err_temp
    
    plt.errorbar(x_plt * 1E-1, y_plt * 1E6, y_err_plt * 1E6, fmt='x', capsize=5, markersize=11, label="Resistivity at constant temperature") 
    params, params_cov = scipy.optimize.curve_fit(fit_func3, x_plt * 1E-1, y_plt* 1E6, sigma = y_err_plt * 1E6, absolute_sigma = True)
    plt.plot(x_plt * 1E-1, fit_func3(x_plt * 1E-1, params[0]), label= r"constant fit $y=a$")
    #plt.axhline(y=params[0], xmin=0.0, xmax=1.0, color='r')
    plt.grid()
    #plt.title("Resistivity in dependence of the magnetic field strength")
    plt.xlabel("Magnetic field [T]")
    plt.ylabel(r"Resistivity [$\mu$m $\Omega$]")
    plt.ylim([213.5, 214.2])
    #plt.legend()

    perr = np.sqrt(np.diag(params_cov))/np.sqrt(len(x_plt))
    print("Resistivity =", params[0], "+/-", perr[0], r"$\mu$m $\Omega$")

    plt.show()
    plt.clf()

    return

def calc_Hall_coefficient():

    d = 350E-6
    d_sigma = 25E-6

    dirName = r"C:\Users\Flo\Desktop\LabCourse\Hall Effect\Hall coefficient data 17.01.20\Flo"
    
    data = {}
    x_plt_1 = []
    y_plt_1 = []
    y_err_plt_1 = []
    x_temp = []
    y_temp_1 = []
    y_err_temp_1 = []
    y_temp_2 = []
    y_err_temp_2 = []
    denom = {}
    num = {}
    I_x = []
    I_x_sigma = []
    I_x_value = 0
    I_x_value_sigma = 0

    excl_datap = 2

    data, n, j = readinvalues_Hall_coeff(dirName)

    #n-2, because we are neglecting the last two points
    x_plt = np.empty(n-excl_datap, dtype = float)
    y_plt_1 = np.empty(n-excl_datap, dtype = float)
    y_err_plt_1 = np.empty(n-excl_datap, dtype = float)
    y_plt_2 = np.empty(n-excl_datap, dtype = float)
    y_err_plt_2 = np.empty(n-excl_datap, dtype = float)
    fit = np.empty(n-excl_datap, dtype = float)

       
    #prints the content of every file in the directory dirName
    #for i in range(n):
    #    print("File ", i + 1)
    #    for j in range(0,29):
    #        print(data[i, j])
    
    #calculate Hall coefficient (removes the first value, because it has a huge error bar)

    for l in range(n-excl_datap):
        #x_plt_1.append(data[l, 0])
        x_temp.append(data[l, 0])
        denom[l, 0] = data[l, 3] * data[l, 0] 
        denom[l, 1] = uncertpropmult(data[l, 0], data[l, 1], data[l, 3], data[l, 4])
        num[l, 0] = 0.25*(calc_V_H(data[l, 9], data[l, 11]) + calc_V_H(data[l, 13], data[l, 15]) + calc_V_H(data[l, 21], data[l, 23]) + calc_V_H(data[l, 25], data[l, 27]))*d
        num[l, 1] = uncertpropmult(num[l, 0]/d, 0.25*uncertpropadd(0.5*uncertpropadd(data[l, 10], data[l, 12], 0, 0), 0.5*uncertpropadd(data[l, 14], data[l, 16], 0, 0), 0.5*uncertpropadd(data[l, 22], data[l, 24], 0, 0), 0.5*uncertpropadd(data[l, 26], data[l, 28], 0, 0)) , d, d_sigma)
        y_temp_1.append(num[l, 0]/denom[l, 0])
        y_err_temp_1.append(uncertpropdiv(num[l, 0], num[l, 1], denom[l, 0], denom[l, 1]))
        y_temp_2.append(0.25*(calc_V_H(data[l, 9], data[l, 11]) + calc_V_H(data[l, 13], data[l, 15]) + calc_V_H(data[l, 21], data[l, 23]) + calc_V_H(data[l, 25], data[l, 27])))
        y_err_temp_2.append(0.25*uncertpropadd(0.5*uncertpropadd(data[l, 10], data[l, 12], 0, 0), 0.5*uncertpropadd(data[l, 14], data[l, 16], 0, 0), 0.5*uncertpropadd(data[l, 22], data[l, 24], 0, 0), 0.5*uncertpropadd(data[l, 26], data[l, 28], 0,0)))
        #print("R_H =", y_plt[l], "+/-", y_err_plt[l])

    x_plt[:] = x_temp
    y_plt_1[:] = y_temp_1
    y_err_plt_1[:] = y_err_temp_1
    y_plt_2[:] = y_temp_2
    y_err_plt_2[:] = y_err_temp_2

    #plots the hall coefficient
    #plt.errorbar(x_plt, y_plt_1 * 1E6, y_err_plt_1 * 1E6, fmt='x', capsize=5, label="Hall coefficient")
    #plt.grid()
    #plt.title("Hall effect in a doped GaAs semiconductor")
    #plt.xlabel("Magnetic field [T]")
    #plt.ylabel(r"Hall coefficient [$cm^3$ $C^{-1}$]")
    #plt.legend()

    ##Plot for paper
    ##neglecting the errorbars, since they are pretty small
    plt.errorbar(x_plt, y_plt_2 * 1E3, fmt='x', label="Hall voltage at constant temperature", markersize=11)
    params, params_cov = scipy.optimize.curve_fit(fit_func1, x_plt, y_plt_2* 1E3, sigma = y_err_plt_2* 1E3, absolute_sigma = True)
    plt.plot(x_plt, fit_func1(x_plt, params[0]), label= r"linear fit $y=ax$")
    plt.grid()
    #plt.title("Hall effect in a doped GaAs semiconductor")
    plt.xlabel("Magnetic field [T]")
    plt.ylabel("Hall voltage [mV]")
    plt.legend()


    fig, (ax1, ax2) = plt.subplots(2)

    #Plots the Hall voltage
    #neglecting the errorbars, since they are pretty small
    ax1.errorbar(x_plt, y_plt_2 * 1E3, fmt='x', markersize=11, capsize=5, label="Measured Hall voltage")
    params, params_cov = scipy.optimize.curve_fit(fit_func1, x_plt, y_plt_2 * 1E3, sigma = y_err_plt_2 * 1E3, absolute_sigma = True)
    ax1.plot(x_plt, fit_func1(x_plt, params[0]), label= r"linear fit $y=ax$")
    ax1.grid()
    #ax1.title.set_text("Hall voltage")
    ax1.set_xlabel("Magnetic field [T]")
    ax1.set_ylabel("Hall voltage [mV]")
    ax1.legend()

    ax2.errorbar(x_plt, y_plt_2 * 1E3 - fit_func1(x_plt, params[0]), y_err_plt_2 * 1E3, fmt='x', capsize=5, markersize=11, label = "Residuals")
    ax2.axhline(y=0.0, xmin=0.0, xmax=1.0, color='r')
    #ax2.title.set_text(r"Residuals of $V_H$ plot")
    ax2.set_xlabel("Magnetic field [T]")
    ax2.set_ylabel("data - fit [mV]")
    ax2.legend()

    #standard deviation = sqrt(covariance matrix elements)
    perr = np.sqrt(np.diag(params_cov))/np.sqrt(len(x_plt))
    print("a =", params[0], "+/-", perr[0], "mV/T")

    #calculation of R_H with the slope in the V_H(B) diagram

    for j in range(0, len(x_plt)):
        I_x.append(data[j, 3])
        I_x_sigma.append(data[j, 4])
        I_x_value += I_x[j]
   

    I_x_value = I_x_value/len(I_x)
    #EXPLUDED THOSE POINTS. INSERT IF INCLUDED
    #+I_x[7]**2+I_x[8]**2
    I_x_value_sigma = 1/(len(I_x))*math.sqrt(I_x[0]**2+I_x[1]**2+I_x[2]**2+I_x[3]**2+I_x[4]**2+I_x[5]**2+I_x[6]**2)
    numerator_sigma = 0
    numerator_sigma = uncertpropmult(params[0], perr[0], d, d_sigma)
    
    print("R_H =", params[0]*d/I_x_value * 1E3, "+/-", uncertpropdiv(params[0]*d, numerator_sigma, I_x_value, I_x_value_sigma) * 1E3, "cm^3/C")

    #Chi squared test

    chi_squared = 0
    fit[:] = fit_func1(x_plt, params[0])

    for l in range(0, len(y_plt_2)):
        chi_squared += (y_plt_2[l] * 1E3 - fit[l])**2/(y_err_plt_2[l] * 1E3)**2 

    print("Chi^2_red =" , chi_squared/(len(y_plt_2) - 1))

    #plt.errorbar(x_plt, y_plt, y_err_plt, fmt='x', capsize=5) 
    #plt.grid()
    #plt.title("Hall coefficient in dependence of the magnetic field")
    #plt.xlabel("B[T]")
    #plt.ylabel("R_H[?]")

    fig.tight_layout()

    plt.show()
    plt.clf()         
   
    return

def calc_tempdep_Hall_coefficient():

    d = 350E-6
    d_sigma = 25E-6
    T_0 = 273.15 #0° Celsius in Kelvin

    dirName = r"C:\Users\Flo\Desktop\LabCourse\Hall Effect\Temperature dependent data 22.01.20\Flo"

    data = {}
    x_temp = []
    y_temp = []
    y_err_temp = []
    denom = {}
    num = {}

    data, i, j = readinvalues_Hall_coeff(dirName)

    x_plt = np.empty(j, dtype = float)
    y_plt = np.empty(j, dtype = float)
    y_err_plt = np.empty(j, dtype = float)
    fit = np.empty(j, dtype = float)

    #print("Number of lines", j)


    for l in range (0, j):
        x_temp.append(data[l, 2] + T_0)
        denom[l, 0] = data[l, 3] * data[l, 0] 
        denom[l, 1] = uncertpropmult(data[l, 0], data[l, 1], data[l, 3], data[l, 4])
        num[l, 0] = 0.25*(calc_V_H(data[l, 9], data[l, 11]) + calc_V_H(data[l, 13], data[l, 15]) + calc_V_H(data[l, 21], data[l, 23]) + calc_V_H(data[l, 25], data[l, 27]))*d
        num[l, 1] = uncertpropmult(num[l, 0]/d, 0.25*uncertpropadd(0.5*uncertpropadd(data[l, 10], data[l, 12], 0, 0), 0.5*uncertpropadd(data[l, 14], data[l, 16], 0, 0), 0.5*uncertpropadd(data[l, 22], data[l, 24], 0, 0), 0.5*uncertpropadd(data[l, 26], data[l, 28], 0, 0)) , d, d_sigma)
        y_temp.append(num[l, 0]/denom[l, 0])
        y_err_temp.append(uncertpropdiv(num[l, 0], num[l, 1], denom[l, 0], denom[l, 1]))
        

    x_plt[:] = x_temp
    y_plt[:] = y_temp
    y_err_plt[:] = y_err_temp

    
    #1E6 factors are added for "optimal" units
    plt.errorbar(x_plt, y_plt * 1E6, y_err_plt * 1E6, fmt='x', capsize=5, label = r"Hall coefficient at constant magnetic field", markersize=11)
    params, params_cov = scipy.optimize.curve_fit(fit_func4, x_plt, y_plt * 1E6, sigma = y_err_plt * 1E6, absolute_sigma = True)
    plt.plot(x_plt, fit_func4(x_plt, params[0], params[1]), label= r"linear fit $y=ax$")
    plt.grid()
    #plt.title("Doped GaAs sample at different temperatures")
    plt.xlabel("Temperature [K]")
    plt.ylabel(r"Hall coefficient [$cm^3$ $C^{-1}$]")
    #plt.legend()

    perr_a = np.sqrt(np.diag(params_cov[0]))/np.sqrt(len(x_plt))
    perr_b = np.sqrt(np.diag(params_cov[1]))/np.sqrt(len(x_plt))

    print("a =", params[0], "+/-", perr_b[0][0])
    print("b =", params[1], "+/-", perr_a[0][0])

    print("n =", (params[0]*293)**(3)*np.exp(-params[1]/293), "+/-", "1/cm^3")

    chi_squared_value = 0
    fit[:] = fit_func4(x_plt, params[0], params[1])

    for k in range(len(x_plt)):
        chi_squared_value += (y_plt[k] * 1E6 - fit[k])**2/(y_err_plt[k] * 1E6)**2 

    print("red. Chi^2 =", chi_squared_value/(len(x_plt)-2))

    plt.show()
    plt.clf()
      
    return

def calc_tempdep_resistivity():

    d_sample = 350E-6 
    d_sample_sigma = 25E-6
    T_0 = 273.15
    Pi = math.pi
    ln = math.log

    data = {}

    R_12_43, R_14_23, d = sym.symbols('R_12_43 R_14_23 d')

    dirName = r"C:\Users\Flo\Desktop\LabCourse\Hall Effect\Temperature dependent data 22.01.20\Flo"

    data, n, j = readinvalues_Hall_coeff(dirName)

    x_temp = []
    y_temp = []
    y_err_temp = []
    x_plt = np.empty(j, dtype = float)
    y_plt = np.empty(j, dtype = float)
    y_err_plt = np.empty(j, dtype = float)
    fit = np.empty(j, dtype= float)

    R_B_12_43 = {}
    R_B_14_23 = {}
    R_B_21_34 = {}
    R_B_41_32 = {}
    rho = {}
    rho_reverse = {}


    for l in range(j):
        R_B_12_43[l, 0] = data[l, 5]/data[l, 3] 
        R_B_12_43[l, 1] = uncertpropdiv(data[l, 5], data[l, 6], data[l, 3], data[l, 4])
        R_B_14_23[l, 0] = data[l, 7]/data[l, 3] 
        R_B_14_23[l, 1] = uncertpropdiv(data[l, 7], data[l, 8], data[l, 3], data[l, 4])
        R_B_21_34[l, 0] = data[l, 17]/data[l, 3]
        R_B_21_34[l, 1] = uncertpropdiv(data[l, 17], data[l, 18], data[l, 3], data[l, 4])
        R_B_41_32[l, 0] = data[l, 19]/data[l, 3]
        R_B_41_32[l, 1] = uncertpropdiv(data[l, 19], data[l, 20], data[l, 3], data[l, 4])

    #we neglected the d_sigma error, because it is a systematic uncertainty
    for l in range(j):
        x_temp.append(data[l, 2] + T_0)
        rho[l, 0] = Pi*d_sample/(2*ln(2))*(R_B_12_43[l, 0]+R_B_14_23[l, 0])*function_f(R_B_12_43[l, 0]/R_B_14_23[l, 0])
        #rho[l, 1] = math.sqrt((deriv_of_rho(R_12_43, R_B_12_43[l, 0], R_B_14_23[l, 0], d_sample)*R_B_12_43[l, 1])**2 + (deriv_of_rho(R_14_23, R_B_12_43[l, 0], R_B_14_23[l, 0], d_sample)*R_B_14_23[l, 1])**2 + (deriv_of_rho(d, R_B_12_43[l, 0], R_B_14_23[l, 0], d_sample)*d_sample_sigma)**2)
        rho[l, 1] = math.sqrt((deriv_of_rho(R_12_43, R_B_12_43[l, 0], R_B_14_23[l, 0], d_sample)*R_B_12_43[l, 1])**2 + (deriv_of_rho(R_14_23, R_B_12_43[l, 0], R_B_14_23[l, 0], d_sample)*R_B_14_23[l, 1])**2)
        rho_reverse[l, 0] = Pi*d_sample/(2*ln(2))*(R_B_21_34[l, 0]+R_B_41_32[l, 0])*function_f(R_B_21_34[l, 0]/R_B_41_32[l, 0])
        #rho_reverse[l, 1] = math.sqrt((deriv_of_rho(R_12_43, R_B_21_34[l, 0], R_B_41_32[l, 0], d_sample)*R_B_21_34[l, 1])**2 + (deriv_of_rho(R_14_23, R_B_21_34[l, 0], R_B_41_32[l, 0], d_sample)*R_B_41_32[l, 1])**2 + (deriv_of_rho(R_14_23, R_B_21_34[l, 0], R_B_41_32[l, 0], d_sample)*d_sample_sigma)**2)
        rho_reverse[l, 1] = math.sqrt((deriv_of_rho(R_12_43, R_B_21_34[l, 0], R_B_41_32[l, 0], d_sample)*R_B_21_34[l, 1])**2 + (deriv_of_rho(R_14_23, R_B_21_34[l, 0], R_B_41_32[l, 0], d_sample)*R_B_41_32[l, 1])**2)
        y_temp.append(0.5*(rho[l, 0] + rho_reverse[l, 0]))
        y_err_temp.append(0.5*uncertpropadd(rho[l, 1], rho_reverse[l, 1], 0, 0))

    x_plt[:] = x_temp
    y_plt[:] = y_temp
    y_err_plt[:] = y_err_temp

    #we neglect the errorbars, since they are too small
    plt.errorbar(x_plt, y_plt * 1E6, fmt='x', label="Resistivity at constant magnetic field", markersize=11) 
    #plt.errorbar(x_plt, y_plt * 1E6, y_err_plt * 1E6, fmt='x', capsize=5, label="Resistivity") 
    params, params_cov = scipy.optimize.curve_fit(fit_func2, x_plt, y_plt * 1E6, sigma = y_err_plt * 1E6, absolute_sigma = True)
    plt.plot(x_plt, fit_func2(x_plt, params[0], params[1], params[2]), label= r"$\rho = a(T-b)^{3/2}+c$ fit")

    #Chi squared calculation

    chi_squared_value = 0
    fit[:] = fit_func2(x_plt, params[0], params[1], params[2])

    for k in range(len(x_plt)):
        #print("rho errors :", y_err_plt[k])
        #print("rel. rho errors :", y_err_plt[k]/y_plt[k])
        chi_squared_value += (y_plt[k] * 1E6 - fit[k])**2/(y_err_plt[k]* 1E6)**2 

    #in the denominator we have (# of data points) - (# of fitting parameters)
    print("Chi^2_red =", chi_squared_value/(len(x_plt) - 2))
    print("a =", params[0])
    print("b =", params[1])
    print("c =", params[2])

    plt.grid()
    #plt.title("GaAs sample at different temperatures")
    plt.xlabel("Temperature [K]")
    plt.ylabel(r"Resistivity [$\mu m$ $\Omega$]")
    plt.legend()

    plt.show()
    plt.clf()    

    return

def main():
        
    #calc_resistivity()
    #calc_Hall_coefficient()
    #calc_tempdep_resistivity()
    calc_tempdep_Hall_coefficient()
    
if __name__ == "__main__" :
    main()