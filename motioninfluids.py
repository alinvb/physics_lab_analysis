# physics_lab_analysis
finding the terminal velocity of a ball falling through fluids of various viscosities 
#define my constants

n = 0.01 # viscosity of the liquid 
density_glycerine = 1.00 # density of the liquid 
D = 0.095 # the width of the box 
import numpy as np
import matplotlib.pyplot as pl

# this loads in the various files. They key is the diamter of each ball
files = {"Ball1-run1-2.txt":  2.21 ,"Ball1-run2-2.txt": 2.21 , "Ball1-run3-2.txt": 2.21}
# creates the list of terminal velocities for each ball 
v_terminal_list = []
# creates the list of errros for each terminal velocity 
v_terminal_errl = []
# this is the list of reynolds numbers corrosoping to each terminal velocity
Rey_term = []
# this list will contain the errors in the reynolds number
Rey_term_err = []
# this list will contain the radius for each ball 
radius = []
# this dictionary will contain the file name as a key and the terminal velocity for each file
v_terminal_dictionary = {}
# this dictionary will contain the file name as a key and the reynolds number for each file 
reynolds_term_dictionary = {}
# this loop iterates over each file in the dictionary 
for key in files:
        # the time and positon data are loaded in from the file 
        time,position = np.loadtxt(key, unpack = True, skiprows = 1)
        # the position is converted to metres
        position = position/1000.0 

        # an empty array to contain the error in the time measurement is created
        dtime = np.empty(len(time))

        # an empty array to contain the error in the position data is created
        dposition = np.empty(len(position))
        # the time error is filled in 
        dtime.fill(0.0005)
        # the position error is filled in
        dposition.fill(0.0005)
        dposition = dposition
        # the l is the diamter of the ball converted to metres
        l = files[key]/1000.0
        # the error is the error in the measurement converted to metres 
        lerr = 0.00002/ 1000.0
        # the inital counter for the loop below 
        i = 1
        # the velocity list will eventually contain the velocity for each time step 
        velocity = np.empty(len(time))
        velocity.fill(0.0)
        # the inital velocity is the inital position divided by the initial time 
        velocity[0] = position[0]/time[0]
        # this list will eventually contain the error in the velocity
        verr = []
        verr.append(0.00002)
        # since the velocity will go to zero for some time intervals where the position appears not to have changed this list will filter those
        velocity_withoutzero = []
        # the corrospoding times
        time_withoutzero = []
        # the inital velocity
        velocity_withoutzero.append(position[0]/time[0])
        time_withoutzero.append(0.0)
        #the acceleration will be determined
        acceleration = np.empty(len(time))
        acceleration.fill(0.0)
        # this while loop will iterate over each position in the array and find the velocity and acceleration
        # velocity is calculated by determing delta(position)/delta(time)
        # acceleration is determined by doing delta(velocity)/delta(time)
        while i < len(position):
                # velocity is determined
                velocity[i] = (position[i] - position[i - 1])/(time[i] - time[i -1])
                
                # acceleration is determined
                acceleration[i-1] = (velocity[i] - velocity[i-1])/(time[i] - time[i-1])
                
                # this calculated the relative error in the position subtraction
                relerr_position = np.sqrt(2)*0.000005/(position[i]+0.5)
                
                # this calculates the relative error in the time subtraction
                relerr_time = (np.sqrt(2)*0.000005)/(time[i])

                relerr = np.sqrt(relerr_position**2 + relerr_time**2)


                # this calculates the error and adds it to the list 
                verr.append(float(relerr) * float(velocity[i]))
                
                i = i + 1

               
        j = 1
        # this loops filters out all the times where the velocity goes to zero because of a slow frame rate 
        while j < len(position):
                if velocity[j] > 0:
                        # appends the velocity with a value greater than zero 
                        velocity_withoutzero.append(velocity[j])
                        time_withoutzero.append(time[j])
                        
                j = j + 1 

        k = 0
        # terminal velocity is when the acceleration goes to 0
        # to determine when the acceleration goes to zero for a prolonged time the average of several accelerations is taken 
        term_v_list = []
        while k < len(acceleration):
                # calculated the average accelerations to determine if they are zero 
                average = (acceleration[k-5] + acceleration[k - 5] + acceleration[k-3] + acceleration[k-2] + acceleration[k-1] + acceleration[k])/5
                # code will accept all accelerations less than 0.1 and greater than -0.1 
                if average < 0.1:
                        if average > -0.1:
                                term_v_list.append(k)
                k = k + 1


        # calculated the correction in the velocity
        
        # calculates the error in the velocity and turns it into an array 
        verr = np.array(verr)

        # takes the average index of zero acceleration to find when the terminal velocity first appears
        start = sum(term_v_list)/len(term_v_list)

        # the terminal velocity is the average of all the velocities after the start index
        


        v_terminal = sum(velocity[start:])/len(velocity[start:])

        print v_terminal
        
        v_velocity = np.array(v_terminal/((1-(2.104*(l/D)) + (2.089*(l/D))**2)))

        print v_velocity

        velocity = np.array(velocity/((1-(2.104*(l/D)) + (2.089*(l/D))**2)))

        # calculates the error in the terminal velocities
        v_terminal_err = np.sqrt(sum((v_terminal - velocity[start:])**2)/(len(velocity[start:]) -1))

        # adds to the dictionary of terminal velocities
        v_terminal_dictionary[key] = [v_terminal, v_terminal_err]

        # determines the reynolds number for all velocities
        reynolds = np.array((density_glycerine * l * velocity)/n)

        # determines the reynolds number for the terminal velocity
        Reynolds_term = (density_glycerine * (files[key]/10) * (v_terminal * 100))/n

        print Reynolds_term 

        # calculates the error in the terminal velocity reynolds number 
        err = np.sqrt(((v_terminal_err/v_terminal)**2) + ((lerr/l)**2))*Reynolds_term

        print err 
        # calculates the error in the reynolds numbers
        err_reynolds = np.sqrt(((verr**2) + ((lerr/(l))**2))) * reynolds

        # adds to the reynolds number dictionarys
        reynolds_term_dictionary[key] = [Reynolds_term, err]

        # creates an array of terminal velocities 
        v_terminal_list.append(v_terminal)
        # creates an array of terminal velocity errors
        v_terminal_errl.append(v_terminal_err)

        
        Rey_term.append(Reynolds_term)
        Rey_term_err.append(err)

        # creates a list with all the radii
        radius.append(files[key])

        # this plots the position, acceleration, velocity, reynolds number vs time for all the data files
        pl.subplot(4,1,1)

        pl.title(files[key])

        
        pl.errorbar(time,position, xerr = dtime, yerr = dposition,linestyle = 'none', marker = '.', color = 'r')

        pl.ylabel("Position(m)")
        
        pl.subplot(4,1,2)
        
        pl.errorbar(time_withoutzero,velocity_withoutzero, linestyle = 'none', marker = '.', color = 'b')

        pl.ylabel("Velocity(m/s)")
        
        pl.subplot(4,1,3)
        
        pl.plot(time,acceleration, linestyle = 'none', marker = '.', color = 'k')

        pl.ylabel("Acceleration(m/s^2)")
        
        pl.subplot(4,1,4)
        
        pl.errorbar(time,reynolds, xerr= dtime, yerr = err_reynolds, linestyle = 'none', marker = '.', color = 'g')

        pl.xlabel("Time(s)")
        pl.ylabel("Reynolds Number")
        
        pl.show()

print v_terminal_dictionary
print reynolds_term_dictionary 


v_terminal_errl = np.array(v_terminal_errl)

v_terminal_list = np.array(v_terminal_list)

radius_1 = np.sqrt((np.array(radius)))

# imports curve fit from scipy to fit the terminal velocity vs radius data 
from scipy.optimize import curve_fit

# creates a linear fit 
def linear(x,b,m):
        return b + x*m

# fits the data 
fitpar, covmat = curve_fit(linear, radius_1, v_terminal_list, sigma = v_terminal_errl)

# the erorr in the slope of the data
sigma_slope = (covmat[1,1])

# creates the fit residuals
fit_residuals = (v_terminal_list - linear(radius_1, *fitpar))

# determines the chi squared value
chi_squared = sum(((fit_residuals)**2)/(v_terminal_errl**2))

dof = len(v_terminal_list) - len(fitpar)

# prints the desired information
print fitpar[1]

print chi_squared

print chi_squared/dof

print dof

print sigma_slope

# plots all the data 
pl.errorbar(radius_1, linear(radius_1, *fitpar), color = 'b')

pl.errorbar(radius_1, fitpar[0] + radius_1 * (fitpar[1] + sigma_slope), linestyle = 'dotted')

pl.errorbar(radius_1, fitpar[0] + radius_1 * (fitpar[1] - sigma_slope), linestyle = 'dotted' ) 

pl.errorbar(radius_1, v_terminal_list, xerr = 0.002, yerr = v_terminal_errl, linestyle = 'none', marker = '.', color = 'r')

pl.title("Velocity vs  square root Radius for Water ")
pl.xlabel("square root Radius(m)")
pl.ylabel("Velocity (mm/s)")
pl.show()

