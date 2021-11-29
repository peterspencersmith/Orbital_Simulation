#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on %(date)s

@author: %(username)s
"""
import numpy as np
import matplotlib.pyplot as plt

G=6.67408E-11 #Gravitational constant
Me=5.972E24 #Mass of the Earth


MyInput = '0'
while MyInput != 'q':
    MyInput = input( 'For a full explanation of the figures produced by the following simulation\
 please see the file "Double body orbits.pdf" located in the GitHub repository located at https://github.com/peterspencersmith/Project-1/blob/main/Orbital%20simulation.pdf \
 \n \nPlease enter one of the following choices, "a", "b" or "q" to quit.\
 Entering "a" runs a simulation of a single body orbit, entering "b" runs a simulation of a double body orbit: ')
    print('\n\nUser has entered the choice: ',MyInput)

    if MyInput == 'a':

        input_co= input('Please enter a value for the initial velocity coeffecient for the rocketship. \
A value of around 1.3 produces an eliptical orbit and a value of 1 a circular orbit: \
')
        co=float(input_co)

        m=1000 #mass of the orbital body
        h=0.5 #time step
        
        #establishes the functions to be called within the while loop
        def f1(vx):
            return vx
        
        def f2(vy):
            return vy
        
        def f3(x,y):
            a=x**2+y**2
            return (-G * Me * x)/(a * np.sqrt(a))
        
        def f4(x,y):
            a=x**2 + y**2
            return (-G * Me * y)/(a * np.sqrt(a))
        
        #The functions which calculate the stable orbital velocity
        def stable_v_t(G, Me, r):
            return np.sqrt(G * Me / np.sqrt(r**2))
        
        # We assume a clockwise orbit 
        def stable_v(G, Me, x, y):
            r = np.sqrt(x**2 + y**2)
            v_t = stable_v_t(G, Me, r)
            return v_t * y / r, v_t * x / r
        
        #Establishes the empty lists to be filled and plotted
        x_values = []
        y_values = []
        ke_values = []
        pe_values = []
        u_values = []
        t_values = []
        
#        Establish the initial conditions..        
        x=0
        y=-7000000
        vx = co*stable_v_t(G, Me, y)
        vy = 0
        t=0
        
        while (t<50000):
            
            k1x = f1(vx)
            k1y = f2(vy)
            k1vx = f3(x,y)
            k1vy = f4(x,y)
            
            k2x = f1(vx+((h*k1vx)/2))
            k2y = f2(vy+((h*k1vy)/2))
            k2vx = f3(x+((h*k1x)/2),y+((h*k1y)/2))
            k2vy = f4(x+((h*k1x)/2),y+((h*k1y)/2))
            
            k3x = f1(vx+((h*k2vx)/2))
            k3y = f2(vy+((h*k2vy)/2))
            k3vx = f3(x+((h*k2x)/2),y+((h*k2y)/2))
            k3vy = f4(x+((h*k2x)/2),y+((h*k2y)/2))
            
            k4x = f1(vx+(h*k3vx))
            k4y = f2(vy+(h*k3vy))
            k4vx = f3(x+(h*k3x),y+(h*k3y))
            k4vy = f4(x+(h*k3x),y+(h*k3y))
            
            x = x + (h/6)*(k1x + 2*k2x + 2*k3x + k4x)
            y = y + (h/6)*(k1y + 2*k2y + 2*k3y + k4y)
            vx = vx + (h/6)*(k1vx + 2*k2vx + 2*k3vx + k4vx)
            vy = vy + (h/6)*(k1vy + 2*k2vy + 2*k3vy + k4vy)
            t = t + h
            
            x_values.append(x)
            y_values.append(y)
            
        #    appends the energy values onto the empty lists
            ke_values.append(0.5*m*(vx**2 + vy**2))
            pe_values.append(-G*Me*m/(np.sqrt(x**2 + y**2)))
            u_values.append(0.5*m*(vx**2 + vy**2) + -G*Me*m/(np.sqrt(x**2 + y**2)))
            t_values.append(t)
        
        #position plot
        plt.figure()
        plt.plot(x_values,y_values)
        plt.title('Plot of position for the single bodied problem')
        plt.ylabel('x(m)')
        plt.xlabel('y(m)')
        plt.show()
        
        #energy plots
        plt.figure()
        plt.plot(t_values,ke_values, label='Kinetic')
        plt.plot(t_values,pe_values, label='Potential')
        plt.plot(t_values,u_values, label='Total')
        plt.legend()
        plt.title('Plot of Energy as a function of time for the single bodied problem')
        plt.ylabel('Energy')
        plt.xlabel('Time (s)')
        plt.show()

    elif MyInput == 'b':

        Mm=7.34767309E22 #Mass of the moon
        y_m=384400E3 #Y-axis position of the moon
        vx0=7545.822473395462 #Stable Earth orbital velocity 
        x=0
        y=-7000000 #Initial position of rocket
        vx = vx0*1.3995 #Initial velocity of rocket
        vy = 0
        t=0
        h=10 #Time step
        
        def f1(vx):
            return vx
        
        def f2(vy):
            return vy
        
        def f3(x,y):
            a=x**2 + y**2
            b=np.abs(x**2+(y-y_m)**2)
            return -(G * Me * x)/(a * np.sqrt(a)) -(G * Mm * x)/(b * np.sqrt(b))
        
        def f4(x,y):
            a=x**2 + y**2
            b=np.abs(x**2+(y-y_m)**2)
            return -(G * Me * y)/(a * np.sqrt(a)) -(G * Mm)*(y-y_m)/(b * np.sqrt(b))
        
        x_values = []
        y_values = []
                
        while (t<845000):
            
            k1x = f1(vx)
            k1y = f2(vy)
            k1vx = f3(x,y)
            k1vy = f4(x,y)
            
            k2x = f1(vx+((h*k1vx)/2))
            k2y = f2(vy+((h*k1vy)/2))
            k2vx = f3(x+((h*k1x)/2),y+((h*k1y)/2))
            k2vy = f4(x+((h*k1x)/2),y+((h*k1y)/2))
            
            k3x = f1(vx+((h*k2vx)/2))
            k3y = f2(vy+((h*k2vy)/2))
            k3vx = f3(x+((h*k2x)/2),y+((h*k2y)/2))
            k3vy = f4(x+((h*k2x)/2),y+((h*k2y)/2))
            
            k4x = f1(vx+(h*k3vx))
            k4y = f2(vy+(h*k3vy))
            k4vx = f3(x+(h*k3x),y+(h*k3y))
            k4vy = f4(x+(h*k3x),y+(h*k3y))
            
            x = x + (h/6)*(k1x + 2*k2x + 2*k3x + k4x)
            y = y + (h/6)*(k1y + 2*k2y + 2*k3y + k4y)
            vx = vx + (h/6)*(k1vx + 2*k2vx + 2*k3vx + k4vx)
            vy = vy + (h/6)*(k1vy + 2*k2vy + 2*k3vy + k4vy)
            t = t + h
            
            x_values.append(x)
            y_values.append(y)
            
        print('\n\nflight time (s) = ',t)
        
        plt.plot(x_values,y_values)
        plt.title('Plot of position for the doubble bodied problem')
        plt.ylabel('x(m)')
        plt.xlabel('y(m)')
        plt.scatter(0,y_m)
        plt.scatter(0,0)
        plt.show()
        
    elif MyInput != 'q':
        print('Sorry, this is not a recognised choice')
        
print('\n\nGoodbye Rocketman..')
    




    