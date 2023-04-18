import paho.mqtt.client as mqtt
import math
import numpy as np
import time
import pandas as pd
from teprob_OOP import teprob

#C               Tennessee Eastman Process Control Test Problem
#C
#C                    James J. Downs and Ernest F. Vogel
#C
#C                  Process and Control Systems Engineering
#C                        Tennessee Eastman Company
#C                              P.O. Box 511
#C                          Kingsport, TN  37662
#C
#C
#C  Reference:
#C    "A Plant-Wide Industrial Process Control Problem"
#C    Presented at the AIChE 1990 Annual Meeting
#C    Industrial Challenge Problems in Process Control, Paper #24a
#C    Chicago, Illinois, November 14, 1990
#C
#C    "A Plant-Wide Industrial Process Control Problem"
#C    Computers and Chemical Engineering, Vol. 17, No. 3, pp. 245-255
#C    (1993).
#C    
#C
#C  Main program for demonstrating application of the Tennessee Eastman
#C  Process Control Test Problem
#C
#C
#C=============================================================================

#C=============================================================================

def on_connect(client, userdata, flags, rc):
    if rc == 0:
        client.connected_flag = True
        print("connected OK")
    else:
        print("Bad connection Returned code=". rc)
        client.bad_connection_flag = True
        
def on_message(client, userdata, message):
    global start,stop,restart,xmv,idv
    global gain,tau_I
    
    temp = float(message.payload.decode("utf-8"))
    if str(message.topic) == "START_PAUSE":
        start = temp
    elif str(message.topic) == "STOP":
        stop = temp
    elif str(message.topic) == "RESTART":
        restart = temp
    for j in range(12):
        if str(message.topic) == ("XMV("+str(j+1)+")"):
            xmv[j] = temp
    for index in range(20):
        if str(message.topic) == ("IDV("+str(index+1)+")"):
            idv[index] = temp
    for m in range(19):
        a = m+1 if m<=10 else m+2
        if str(message.topic) == ("gain_"+str(a)):
            gain[m] = temp
        elif str(message.topic) == ("tau_I_"+str(a)):
            tau_I[m] = temp

def contrl(xmeas, xmv, setpt, gain_, delta_t, tau_I_, errold_):
    #C  Proportional-Integral Controller (Velocity Form)
    #C  gain_ = Controller gain_
    #C  tau_I_ = Reset Time (min)
    err = setpt - xmeas[14]
    dxmv = gain_*((err-errold_) + err*delta_t*60.0/tau_I_)
    xmv7 = xmv[7] - dxmv
    errold_a = err
    return xmv7, errold_a

def contrl_1(xmeas, xmv, setpt, gain,delta_t, tau_I, errold):
    #global xmv, errold_1
    err1 = (setpt[0] - xmeas[1]) * 100. / 5811.
    dxmv = gain[1-1] * ((err1 - errold[1-1]) + (err1*delta_t*3./tau_I[1-1] if tau_I[1-1] != 0 else 0))
    xmv0 = xmv[0] + dxmv
    errold_1a = err1
    return xmv0, errold_1a

def contrl_2(xmeas, xmv, setpt, gain,delta_t, tau_I,errold):
    #global xmv, errold_2
    err2 = (setpt[1] - xmeas[2]) * 100. / 8354.
    dxmv = gain[2-1] * (( err2 - errold[2-1] ) + (err2*delta_t*3./tau_I[2-1] if tau_I[2-1] != 0 else 0))
    xmv1 = xmv[1] + dxmv    
    errold_2a = err2
    return xmv1, errold_2a

def contrl_3(xmeas, xmv, setpt, gain,delta_t, tau_I,errold):
    err3 = (setpt[2] - xmeas[0]) * 100. / 1.017
    dxmv = gain[3-1] * (( err3 - errold[3-1] ) + (err3*delta_t*3./tau_I[3-1] if tau_I[3-1] != 0 else 0))
    xmv2 = xmv[2] + dxmv
    errold_3a = err3
    return xmv2, errold_3a

def contrl_4(xmeas, xmv, setpt, gain,delta_t, tau_I,errold):
    err4 = (setpt[3] - xmeas[3]) * 100. / 15.25
    dxmv = gain[4-1] * (( err4 - errold[4-1] ) + (err4*delta_t*3./tau_I[4-1] if tau_I[4-1] != 0 else 0))
    xmv3 = xmv[3] + dxmv
    errold_4a = err4
    return xmv3, errold_4a

def contrl_5(xmeas, xmv, setpt, gain, delta_t, tau_I,errold):    
    err5 = (setpt[4] - xmeas[4])  * 100. / 53.
    dxmv = gain[5-1]*((err5 - errold[5-1]) + err5*delta_t*3./tau_I[5-1])
    xmv4 = xmv[4] + dxmv 
    errold_5a = err5
    return xmv4, errold_5a

def contrl_6(xmeas, xmv, setpt, gain, delta_t, tau_I, errold):
    #global flag
    flag = int(0)
    if (xmeas[12] >= 2950.0):
        flag = 1 
        xmv[5] = 100.0
    elif (flag==1 and xmeas[12]>=2633.7):
        xmv[5] = 100.0
    elif (flag==1 and xmeas[12]<=2633.7):
        xmv[5]=  40.060
        setpt[5] = 0.33712
        errold_6 = 0.0
        flag = 0
    elif (xmeas[12] <= 2300.):
        xmv[5] = 0.0
        flag = 2
    elif (flag==2 and xmeas[12]<=2633.7):
        xmv[5]=0.0
    elif (flag==2 and xmeas[12]>=2633.7):
        xmv[5] = 40.060
        setpt[5] = 0.33712
        errold_6 = 0.0
        flag = 0
    else:
        flag=0
    err6 = (setpt[5] - xmeas[9]) * 100. /1.
    dxmv = gain[6-1] * ((err6 - errold[6-1]) + (err6*delta_t*3./tau_I[6-1] if tau_I[6-1] != 0 else 0))
    xmv5 = xmv[5] + dxmv
    errold_6a = err6
    return xmv5, errold_6a

def contrl_7(xmeas, xmv, setpt, gain, delta_t, tau_I, errold):
    err7 = (setpt[6] - xmeas[11]) * 100./70.
    dxmv = gain[7-1]*((err7 - errold[7-1]) + (err7*delta_t*3./tau_I[7-1] if tau_I[7-1] != 0 else 0))
    xmv6 = xmv[6] + dxmv
    errold_7a = err7
    return xmv6, errold_7a

def contrl_8(xmeas, xmv, setpt, gain, delta_t, tau_I,errold):
    err8 = (setpt[7] - xmeas[14])*100./70.
    dxmv =  gain[8-1]*((err8 - errold[8-1]) + (err8*delta_t*3./tau_I[8-1] if tau_I[8-1] != 0 else 0))
    xmv7 = xmv[7] + dxmv
    errold_8a = err8
    return xmv7, errold_8a

def contrl_9(xmeas, xmv, setpt, gain, delta_t, tau_I,errold):
    err9 = (setpt[8] - xmeas[18])*100./460. 
    dxmv = gain[9-1]*((err9 - errold[9-1]) + (err9*delta_t*3./tau_I[9-1] if tau_I[9-1] != 0 else 0))
    xmv8 = xmv[8] + dxmv
    errold_9a = err9
    return xmv8, errold_9a

def contrl_10(xmeas, xmv, setpt, gain,delta_t, tau_I,errold):
    err10 = (setpt[9] - xmeas[20])*100./150.
    dxmv = gain[10-1]*((err10 - errold[10-1]) + err10*delta_t*3./tau_I[10-1])
    xmv9 = xmv[9] + dxmv
    errold_10a = err10
    return xmv9, errold_10a

def contrl_11(xmeas, xmv, setpt, gain,delta_t,tau_I,errold):
    err11 = (setpt[10] - xmeas[16]) * 100. / 46.
    dxmv = gain[11-1]*((err11 - errold[11-1]) + err11*delta_t*3./tau_I[11-1])
    xmv10 = xmv[10] + dxmv  
    errold_11a = err11
    return xmv10, errold_11a

def contrl_13(xmeas, xmv, setpt, gain,delta_t, tau_I,errold):
    err13 = (setpt[12] - xmeas[22])*100./100.
    dxmv = gain[13-2] * ((err13 - errold[13-2])+err13*delta_t*360./tau_I[13-2])
    setpt2 = setpt[2] + dxmv*1.017/100. 
    errold_13a = err13
    return setpt2, errold_13a

def contrl_14(xmeas, xmv, setpt, gain,delta_t, tau_I,errold):
    err14 = (setpt[13] - xmeas[25])*100./100.
    dxmv = gain[14-2]*((err14 - errold[14-2]) + err14*delta_t*360./tau_I[14-2])
    setpt0 = setpt[0] + dxmv*5811./100.
    errold_14a = err14
    return setpt0, errold_14a

def contrl_15(xmeas, xmv, setpt, gain,delta_t, tau_I,errold):
    err15 = (setpt[14] - xmeas[26]) * 100./100.
    dxmv = gain[15-2]*((err15 - errold[15-2]) + err15*delta_t*360./tau_I[15-2])
    setpt1 = setpt[1] + dxmv*8354./100.
    errold_15a = err15
    return setpt1, errold_15a

def contrl_16(xmeas, xmv, setpt, gain,delta_t, tau_I,errold):
    err16 = (setpt[15] - xmeas[17])*100./130.
    dxmv = gain[16-2] * ((err16 - errold[16-2]) + err16*delta_t*3./tau_I[16-2])
    setpt8 = setpt[8] + dxmv*460./100.
    errold_16a = err16
    return setpt8, errold_16a

def contrl_17(xmeas, xmv, setpt, gain,delta_t, tau_I,errold):
    err17 = (setpt[16] - xmeas[7])*100./50.
    dxmv = gain[17-2]*((err17 - errold[17-2]) + err17*delta_t*3./tau_I[17-2])
    setpt3 = setpt[3] + dxmv*15.25/100.
    errold_17a = err17
    return setpt3, errold_17a

def contrl_18(xmeas, xmv, setpt, gain,delta_t, tau_I,errold):
    err18 = (setpt[17] - xmeas[8])*100./150. 
    dxmv = gain[18-2]*((err18 - errold[18-2]) + err18*delta_t*3./tau_I[18-2])
    setpt9 = setpt[9] + dxmv*150./100.
    errold_18a = err18
    return setpt9, errold_18a

def contrl_19(xmeas, xmv, setpt, gain,delta_t, tau_I,errold):
    err19 = (setpt[18] - xmeas[29])*100./26.
    dxmv = gain[19-2]*((err19 - errold[19-2]) + err19*delta_t*360./tau_I[19-2])
    setpt5 = setpt[5] + dxmv * 1. / 100.
    errold_19a = err19
    return setpt5, errold_19a

def contrl_20(xmeas, xmv, setpt, gain,delta_t, tau_I,errold):
    err20 = (setpt[19] - xmeas[37])*100./1.6
    dxmv = gain[20-2]*((err20 - errold[20-2]) + err20*delta_t*900./tau_I[20-2])
    setpt15 = setpt[15] + dxmv  * 130./100.
    errold_20a = err20
    return setpt15, errold_20a

def contrl_22(xmeas, xmv, setpt, gain,delta_t, tau_I,errold):
    err22 = setpt[11] - xmeas[12]
    dxmv = gain[22-3]*((err22 - errold[22-3]) + err22*delta_t*3./tau_I[22-3])
    xmv5 = xmv[5] + dxmv
    errold_22a = err22
    return xmv5, errold_22a


def OUTPUT(i,xmeas):
    #print("Reac Temp = {0}, Stripper Lev = {1}, Sripper Underflow = {2}".format(xmeas[8], xmeas[14], xmv[7]))
    #print("Reac Press = {0}, Reac Lev = {1}, Reac Temp = {2}".format(xmeas[6], xmeas[7], xmeas[8]))
    #print("xmeas(37) :",xmeas[38], ", xmeas(24): ",xmeas[23])
    print("Time = {0}, A Feed Flow (Stream 1) = {1}, D Feed Flow (Stream 2) = {2} ".format(i, xmeas[0], xmeas[1]))
    return 0

#C=============================================================================

def intgtr(tep, Time,nn, yy, yp, xmv, idv, delta_t):
    K1 = np.zeros(nn)
    K2 = np.zeros(nn)
    K3 = np.zeros(nn)
    K4 = np.zeros(nn)
    K5 = np.zeros(nn)
    K6 = np.zeros(nn)
    K7 = np.zeros(nn)

#C  Forward-Euler Integration Algorithm
    yp = tep.tefunc(nn, Time, yy, yp, xmv , idv)
    Time = Time + delta_t
    for i in range(1-1, nn):
        yy[i] = yy[i] + yp[i]*delta_t

#   ------------------------------------

#   Improved-Euler Integration Algorithm
#    yp= tep.tefunc(Time, yy, yp, xmv , idv)
#    yy_0 = yy
#    Time0 = Time
#    K1 = yp
#    Time = Time + delta_t
#    for i in range(0,nn):
#        yy[i] += K1[i]*delta_t
#    yp= tep.tefunc(Time, yy, yp, xmv , idv)
#    K2 = yp
#    Time = Time0 + delta_t
#    for i in range(0,nn):
#        yy[i] = yy_0[i] + (K1[i] + K2[i]*delta_t)*delta_t
        
#   ------------------------------------

#   ------------------------------------

#   Range-Kutta Integration Algorithm
#    yp = tep.tefunc(Time, yy, yp, xmv , idv)
#    yy_0 = yy
#    Time0 = Time
#    K1 = yp
#    Time = Time0 + delta_t/2
#    for i in range(0,nn):
#        yy[i] = yy_0[i] + K1[i]*delta_t/2
#    yp = tep.tefunc(Time, yy, yp, xmv , idv)
#    K2 = yp
#    Time = Time0 + delta_t/2
#    for i in range(0,nn):
#        yy[i] = yy_0[i] + K2[i]*delta_t/2
#    yp = tep.tefunc(Time, yy, yp, xmv , idv)
#    K3 = yp
#    Time = Time0 + delta_t
#    for i in range(0,nn):
#        yy[i] = yy_0[i] + K3[i]*delta_t
#    yp = tep.tefunc(Time, yy, yp, xmv , idv)
#    K4 = yp
#    Time = Time0 + delta_t
#    for i in range(0,nn):
#        yy[i] = yy_0[i] + (K1[i] + 2*K2[i] + 2*K3[i] + K4[i])*delta_t/6
    
#   ------------------------------------

#   Bogacki-Shampine (Approximation) method
#    yp = tep.tefunc(Time, yy, yp, xmv , idv)
#    yy_0 = yy
#    Time0 = Time
#    K1 = yp
#    Time = Time0 + delta_t/2
#    for i in range(0,nn):
#        yy[i] = yy_0[i] + K1[i]*delta_t/2
#    yp = tep.tefunc(Time, yy, yp, xmv , idv)
#    K2 = yp
#    Time = Time0 + 3*delta_t/4
#    for i in range(0,nn):
#        yy[i] = yy_0[i] + 3*K2[i]*delta_t/4
#    yp = tep.tefunc(Time, yy, yp, xmv , idv)
#    K3 = yp
#    Time = Time0 + delta_t
#    for i in range(0,nn):
#        yy[i] = yy_0[i] + 2*K1[i]*delta_t/9 + K2[i]*delta_t/3 + 4*K3[i]*delta_t/9
#    yp = tep.tefunc(Time, yy, yp, xmv , idv)
#    K4 = yp
#    Time = Time0 + delta_t
#    for i in range(0,nn):
 #       yy[i] = yy_0[i] + (7*K1[i]/24 + K2[i]/3 + K3[i]/3 + K4[i]/8)*delta_t
    
#   ------------------------------------
    
#   Dormand-Prince method
#    yy_0 = yy
#    Time0 = Time
#    yp = tep.tefunc(Time, yy, yp, xmv , idv)
#    K1 = yp
#    Time = Time0 + delta_t/5
#    for i in range(0,nn):
#        yy[i] = yy_0[i] + K1[i]*delta_t/5
#    yp = tep.tefunc(Time, yy, yp, xmv , idv)
#    K2 = yp
#    Time = Time0 + 3*delta_t/10
#    for i in range(0,nn):
#        yy[i] = yy_0[i] + (3*K1[i]/40 + 9*K2[i]/40)*delta_t
#    yp = tep.tefunc(Time, yy, yp, xmv , idv)
#    K3 = yp
#    Time = Time0 + 4*delta_t/5
#    for i in range(0,nn):
#        yy[i] = yy_0[i] + (44*K1[i]/45 + 56*K2[i]/15 + 32*K3[i]/9)*delta_t/9
#    yp = tep.tefunc(Time, yy, yp, xmv , idv)
#    K4 = yp
#    Time = Time0 + 8*delta_t/9
#    for i in range(0,nn):
#        yy[i] = yy_0[i] + (19372*K1[i]/6561 + 25360*K2[i]/2187 + 64448*K3[i]/6561 + \
#                212*K4[i]/729)*delta_t
#    yp = tep.tefunc(Time, yy, yp, xmv , idv)
#    K5 = yp
#    Time = Time0 + delta_t
#    for i in range(0,nn):
#        yy[i] = yy_0[i] + (9017*K1[i]/3168 + 355*K2[i]/33 + 46732*K3[i]/5247 + \
#                49*K4[i]/176 + 5103*K5[i]/18656)*delta_t
#    yp = tep.tefunc(Time, yy, yp, xmv , idv)
#    K6 = yp
#    Time = Time0 + delta_t
#    for i in range(0,nn):
#        yy[i] = yy_0[i] + (35*K1[i]/354 + 500*K3[i]/1113 + 125*K4[i]/192 + \
#                2187*K5[i]/6784 + 11*K6[i]/84)*delta_t
#    yp = tep.tefunc(Time, yy, yp, xmv , idv)
#    K7 = yp
#    Time = Time0 + delta_t
#    for i in range(0,nn):
#        yy[i] = yy_0[i] + (5179*K1[i]/57600 + 7571*K3[i]/16695 + 393*K4[i]/640 + \
#                92097*K5[i]/339200 + 187*K6[i]/2100 + K7[i]/40)*delta_t      
    
    return Time, yy, yp
#   ------------------------------------ 

    
def conshand(xmv_):
    for i in range(0, 11):
        if (xmv_[i] <= 0.0):
            xmv_[i]=0.0
        if (xmv_[i] >= 100.0):
            xmv_[i]= 100.0
    return xmv_


def main():
    global start,stop,restart, xmv, idv
    global gain, tau_I
    
    # Deklarasi simulation
    start = 0
    stop = 0
    restart = 0
    topik1 = ["XMEAS(1)","XMEAS(2)","XMEAS(3)","XMEAS(4)","XMEAS(5)","XMEAS(6)",
            "XMEAS(7)","XMEAS(8)","XMEAS(9)","XMEAS(10)","XMEAS(11)","XMEAS(12)",
            "XMEAS(13)","XMEAS(14)","XMEAS(15)","XMEAS(16)","XMEAS(17)","XMEAS(18)",
            "XMEAS(19)","XMEAS(20)","XMEAS(21)","XMEAS(22)","XMEAS(23)","XMEAS(24)",
            "XMEAS(25)","XMEAS(26)","XMEAS(27)","XMEAS(28)","XMEAS(29)","XMEAS(30)",
            "XMEAS(31)","XMEAS(32)","XMEAS(33)","XMEAS(34)","XMEAS(35)","XMEAS(36)",
            "XMEAS(37)","XMEAS(38)","XMEAS(39)","XMEAS(40)","XMEAS(41)"]
    topik2 = ["XMV(1)","XMV(2)","XMV(3)","XMV(4)","XMV(5)","XMV(6)","XMV(7)","XMV(8)",
            "XMV(9)","XMV(10)","XMV(11)","XMV(12)"]
    topik3 = ["IDV(1)","IDV(2)","IDV(3)","IDV(4)","IDV(5)","IDV(6)","IDV(7)","IDV(8)",
            "IDV(9)","IDV(10)","IDV(11)","IDV(12)","IDV(13)","IDV(14)","IDV(15)","IDV(16)",
            "IDV(17)","IDV(18)","IDV(19)","IDV(20)"]
    
    
    # Number of state simulation
    #C  Set the number of differential equations (states).  The process has 50
    #C  states.  If the user wishes to integrate additional states, nn must be
    #C  increased by the number of additional differential equations.
    nn = 50

    Time = float()
    yy = np.zeros(nn)
    yp = np.zeros(nn)
    setpt = np.zeros(20)

    #C  Set the number of points to simulate
    npts = 34300

    #C  Set the number of pints to simulate in steady state operation
    sspts = 3600 * 8

    #C  Integrator Step Size:  1 Second Converted to Hours
    delta_t = 1.0/3600.0
    
    #C  Initialize Process
    #C  (Sets Time to zero)
    tep = teprob()
    yy, yp = tep.teinit(nn, Time, yy, yp)
    
    xmeas, xmv, idv = tep.getOutput()
    
    #C  Set Controller Parameters
    #C  Make a Stripper Level Set Point Change of +15%
    #setpt = xmeas[14] + 15.0
    #gain_ = 2.0
    #tau_I_ = 5.0
    #errold_ = 0.0
    #C  Example Disturbance:
    #C  Change Reactor Cooling
    #xmv[9] = 38.0

    setpt[0] = 3664.0  
    gain_1 = 1.0
    setpt[1] = 4509.3
    gain_2 = 1.0
    setpt[2] = .25052
    gain_3 = 1.
    setpt[3] = 9.3477
    gain_4 = 1.
    setpt[4] = 26.902
    gain_5 = -0.083        
    tau_I_5 = 1./3600.   
    setpt[5] = 0.33712  
    gain_6 = 1.22  
    setpt[6] = 50.0
    gain_7 = -2.06
    setpt[7] = 50.0
    gain_8 = -1.62
    setpt[8] = 230.31
    gain_9 = 0.41    
    setpt[9] = 94.599
    gain_10 = -0.156 * 10.
    tau_I_10 = 1452./3600. 
    setpt[10] = 22.949    
    gain_11 = 1.09    
    tau_I_11 = 2600./3600.
    setpt[12] = 32.188
    gain_13 = 18.  
    tau_I_13 = 3168./3600.   
    setpt[13] = 6.8820
    gain_14 = 8.3     
    tau_I_14 = 3168.0/3600.
    setpt[14] = 18.776          
    gain_15 = 2.37      
    tau_I_15 = 5069./3600.    
    setpt[15] = 65.731
    gain_16 = 1.69 / 10.
    tau_I_16 = 236./3600.
    setpt[16] = 75.000
    gain_17 = 11.1 / 10.
    tau_I_17 = 3168./3600.  
    setpt[17] = 120.40
    gain_18 = 2.83 * 10.
    tau_I_18 = 982./3600.
    setpt[18] = 13.823
    gain_19 = -83.2 / 5. /3.  
    tau_I_19 = 6336./3600. 
    setpt[19] = 0.83570  
    gain_20 = -16.3 / 5.       
    tau_I_20 = 12408./3600.  
    setpt[11] = 2633.7
    gain_22 = -1.0  * 5.       
    tau_I_22 = 1000./3600.
    
    gain = [gain_1,gain_2,gain_3,gain_4,gain_5,gain_6,gain_7,gain_8,
            gain_9,gain_10,gain_11,gain_13,gain_14,gain_15,gain_16,
            gain_17,gain_18,gain_19,gain_20]
    tau_I = [0, 0, 0, 0, tau_I_5, 0, 0, 0,
            0,tau_I_10, tau_I_11, tau_I_13, tau_I_14, tau_I_15,tau_I_16,
            tau_I_17, tau_I_18, tau_I_19,tau_I_20]
    errold = np.zeros(20)
    
    flag = int()

    #C    Example Disturbance:
    #C    Change Manipulated Variable
    xmv[0] = 63.053 + 0.
    xmv[1] = 53.980 + 0.
    xmv[2] = 24.644 + 0.    
    xmv[3] = 61.302 + 0.
    xmv[4] = 22.210 + 0.
    xmv[5] = 40.064 + 0.
    xmv[6] = 38.100 + 0.
    xmv[7] = 46.534 + 0.
    xmv[8] = 47.446 + 0.
    xmv[9] = 41.106 + 0.
    xmv[10]= 18.114 + 0.

    #C  Set all Disturbance Flags to OFF
    for i in range(1-1, 20):
        idv[i] = 0  
        
    a = int(0)
    mv1 = np.array([])
    mv2 = np.array([])
    mv3 = np.array([])
    me01 = np.array([])
    me02 = np.array([])
    me03 = np.array([])
    me04 = np.array([])
    me05 = np.array([])
    me06 = np.array([])
    me07 = np.array([])
    me08 = np.array([])
    me09 = np.array([])
    me10 = np.array([])
    me11 = np.array([])
    
    broker = "127.0.0.1"      # Local Host
    client = mqtt.Client('Tennessee Eastmann Process')
    client.on_connect = on_connect  # bind call back function
    client.on_message = on_message
    #client.on_message = on_message  #attach function to callback
    print("Connecting to broker", broker)
    client.connect(broker,1883,60)
    client.loop_start() 
    
    k = 1
    #C  Simulation Loop
    while k!=0:
        client.subscribe("START_PAUSE")
        if start == 1:
            xmeas, xmv, idv = tep.getOutput()
            # active Disturbance
            if (k >= sspts):
                idv[19]= 0

            # Control aktif 3 detik sekali    
            #test = np.remainder(i,3)
            test = k%3
            if (test == 0):
                xmv[0], errold[1-1] = contrl_1(xmeas, xmv, setpt, gain, delta_t, tau_I, errold)
                xmv[1], errold[2-1] = contrl_2(xmeas, xmv, setpt, gain, delta_t, tau_I,errold)
                xmv[2], errold[3-1] = contrl_3(xmeas, xmv, setpt, gain, delta_t, tau_I,errold)
                xmv[3], errold[4-1] = contrl_4(xmeas, xmv, setpt, gain, delta_t, tau_I,errold)
                xmv[4], errold[5-1] = contrl_5(xmeas, xmv, setpt, gain, delta_t, tau_I,errold)
                xmv[5], errold[6-1] = contrl_6(xmeas, xmv, setpt, gain, delta_t, tau_I, errold)
                xmv[6], errold[7-1] = contrl_7(xmeas, xmv, setpt, gain, delta_t, tau_I, errold)
                xmv[7], errold[8-1] = contrl_8(xmeas, xmv, setpt, gain, delta_t, tau_I,errold)
                xmv[8], errold[9-1] = contrl_9(xmeas, xmv, setpt, gain, delta_t, tau_I,errold)
                xmv[9], errold[10-1] = contrl_10(xmeas, xmv, setpt, gain,delta_t, tau_I,errold)
                xmv[10], errold[11-1] = contrl_11(xmeas, xmv, setpt, gain,delta_t, tau_I,errold)
                setpt[8], errold[16-2] = contrl_16(xmeas, xmv, setpt, gain,delta_t, tau_I,errold)
                setpt[3], errold[17-2] = contrl_17(xmeas, xmv, setpt, gain,delta_t, tau_I,errold)
                setpt[9], errold[18-2] = contrl_18(xmeas, xmv, setpt, gain,delta_t, tau_I,errold)

            # Control aktif 360 detik sekali
            test_1 = k%360
            #test_1 = np.remainder(i,360)
            if (test_1 == 0):
                setpt[2], errold[13-2] = contrl_13(xmeas, xmv, setpt, gain,delta_t, tau_I,errold)
                setpt[0], errold[14-2] = contrl_14(xmeas, xmv, setpt, gain,delta_t, tau_I,errold)
                setpt[1], errold[15-2] = contrl_15(xmeas, xmv, setpt, gain,delta_t, tau_I,errold)     
                setpt[5], errold[19-2] = contrl_19(xmeas, xmv, setpt, gain,delta_t, tau_I,errold)

            # Control aktif 900 detik sekali    
            test_2 = k%900
            #test_2 = np.remainder(i,900)
            if (test_2 == 0):
                setpt[15], errold[20-2] = contrl_20(xmeas, xmv, setpt, gain,delta_t, tau_I,errold)

            #TEST3 = k%5000
            #TEST3 = np.remainder(i,5000)       
            #if(TEST3 == 0):
                #print('Simulation Time (in seconds) = ', k)

            # Sampling data pada 1 detik
            #test_4 = np.remainder(i,1)
            test_4 = k%1
            if test_4 == 0: 
                for d in range(41):
                    xmeas[d] = np.format_float_scientific(xmeas[d],precision=4)
                    #xmeas[d] = np.around(xmeas[d], decimals=4)
                for c in range(12):
                    xmv[c] = np.format_float_scientific(xmv[c],precision=4)
                    #xmv[c] = np.around(xmv[c], decimals=4)
                mv1 = np.append(mv1,[xmv[0],xmv[1],xmv[2],xmv[3]])
                mv2 = np.append(mv2,[xmv[4],xmv[5],xmv[6],xmv[7]])
                mv3 = np.append(mv3,[xmv[8],xmv[9],xmv[10],xmv[11]])
                me01 = np.append(me01,[xmeas[0],xmeas[1],xmeas[2],xmeas[3]])
                me02 = np.append(me02,[xmeas[4],xmeas[5],xmeas[6],xmeas[7]])
                me03 = np.append(me03,[xmeas[8],xmeas[9],xmeas[10],xmeas[11]])
                me04 = np.append(me04,[xmeas[12],xmeas[13],xmeas[14],xmeas[15]])
                me05 = np.append(me05,[xmeas[16],xmeas[17],xmeas[18],xmeas[19]])
                me06 = np.append(me06,[xmeas[20],xmeas[21],xmeas[22],xmeas[23]])
                me07 = np.append(me07,[xmeas[24],xmeas[25],xmeas[26],xmeas[27]])
                me08 = np.append(me08,[xmeas[28],xmeas[29],xmeas[30],xmeas[31]])
                me09 = np.append(me09,[xmeas[32],xmeas[33],xmeas[34],xmeas[35]])
                me10 = np.append(me10,[xmeas[36],xmeas[37],xmeas[38],xmeas[39]])
                me11 = np.append(me11, xmeas[40])
                a += 1
                    
                OUTPUT(k,xmeas)
                
                for i in range(41):
                    client.publish(topik1[i],xmeas[i])  # publish XMEAS

                for j in range(12):
                    client.subscribe(topik2[j])    # subscribe XMV
                    client.publish(topik2[j],xmv[j])    # publish XMV

                for m in range(20):
                    client.subscribe(topik3[m])    # subscribe IDV
                    client.publish(topik3[m],float(idv[m]))    # publish IDV
                
                for n in range(19):
                    client.subscribe("gain_"+str(n+1 if n<=10 else n+2))
                    client.publish("gain_"+str(n+1 if n<=10 else n+2), float(gain[n]))
                    
                    client.subscribe("tau_I_"+str(n+1 if n<=10 else n+2))
                    client.publish("tau_I_"+str(n+1 if n<=10 else n+2), float(tau_I[n]))

            Time, yy, yp = intgtr(tep, Time,nn, yy, yp,xmv, idv, delta_t)
            xmv = conshand(xmv)
            k += 1
        
        client.subscribe("RESTART")
        if restart == 1:
            yy, yp = tep.teinit(nn, Time, yy, yp)
            xmeas, xmv, idv = tep.getOutput()
            for j in range(12):
                client.publish(topik2[j],xmv[j])    # publish XMV
            for m in range(20):
                client.publish(topik3[m],float(idv[m]))    # publish IDV
            k = 1
            restart = 0
        
        client.subscribe("STOP")
        if stop == 1.:
            break
        
        time.sleep(1)
        
    print("Simulation is done")
    client.loop_stop()
    client.disconnect()
    
    mv1.resize(a,4)
    mv2.resize(a,4)
    mv3.resize(a,4)
    me01.resize(a,4)
    me02.resize(a,4)
    me03.resize(a,4)
    me04.resize(a,4)
    me05.resize(a,4)
    me06.resize(a,4)
    me07.resize(a,4)
    me08.resize(a,4)
    me09.resize(a,4)
    me10.resize(a,4)
    me11.resize(a,1)

    df_mv1 = pd.DataFrame(mv1)
    df_mv1.to_csv(r'C:/Users/Lenovo/Documents/TUGAS AKHIR/Dataset/Data 12 Jam/Tuning TC 10, Eksternal IDV 11, dan IDV 4/mv1.csv', header=False, index=False)
    df_mv2 = pd.DataFrame(mv2)
    df_mv2.to_csv(r'C:/Users/Lenovo/Documents/TUGAS AKHIR/Dataset/Data 12 Jam/Tuning TC 10, Eksternal IDV 11, dan IDV 4/mv2.csv', header=False, index=False)
    df_mv3 = pd.DataFrame(mv3)
    df_mv3.to_csv(r'C:/Users/Lenovo/Documents/TUGAS AKHIR/Dataset/Data 12 Jam/Tuning TC 10, Eksternal IDV 11, dan IDV 4/mv3.csv', header=False, index=False)
    df_me01 = pd.DataFrame(me01)
    df_me01.to_csv(r'C:/Users/Lenovo/Documents/TUGAS AKHIR/Dataset/Data 12 Jam/Tuning TC 10, Eksternal IDV 11, dan IDV 4/me01.csv', header=False, index=False)
    df_me02 = pd.DataFrame(me02)
    df_me02.to_csv(r'C:/Users/Lenovo/Documents/TUGAS AKHIR/Dataset/Data 12 Jam/Tuning TC 10, Eksternal IDV 11, dan IDV 4/me02.csv', header=False, index=False)
    df_me03 = pd.DataFrame(me03)
    df_me03.to_csv(r'C:/Users/Lenovo/Documents/TUGAS AKHIR/Dataset/Data 12 Jam/Tuning TC 10, Eksternal IDV 11, dan IDV 4/me03.csv', header=False, index=False)
    df_me04 = pd.DataFrame(me04)
    df_me04.to_csv(r'C:/Users/Lenovo/Documents/TUGAS AKHIR/Dataset/Data 12 Jam/Tuning TC 10, Eksternal IDV 11, dan IDV 4/me04.csv', header=False, index=False)
    df_me05 = pd.DataFrame(me05)
    df_me05.to_csv(r'C:/Users/Lenovo/Documents/TUGAS AKHIR/Dataset/Data 12 Jam/Tuning TC 10, Eksternal IDV 11, dan IDV 4/me05.csv', header=False, index=False)
    df_me06 = pd.DataFrame(me06)
    df_me06.to_csv(r'C:/Users/Lenovo/Documents/TUGAS AKHIR/Dataset/Data 12 Jam/Tuning TC 10, Eksternal IDV 11, dan IDV 4/me06.csv', header=False, index=False)
    df_me07 = pd.DataFrame(me07)
    df_me07.to_csv(r'C:/Users/Lenovo/Documents/TUGAS AKHIR/Dataset/Data 12 Jam/Tuning TC 10, Eksternal IDV 11, dan IDV 4/me07.csv', header=False, index=False)
    df_me08 = pd.DataFrame(me08)
    df_me08.to_csv(r'C:/Users/Lenovo/Documents/TUGAS AKHIR/Dataset/Data 12 Jam/Tuning TC 10, Eksternal IDV 11, dan IDV 4/me08.csv', header=False, index=False)
    df_me09 = pd.DataFrame(me09)
    df_me09.to_csv(r'C:/Users/Lenovo/Documents/TUGAS AKHIR/Dataset/Data 12 Jam/Tuning TC 10, Eksternal IDV 11, dan IDV 4/me09.csv', header=False, index=False)
    df_me10 = pd.DataFrame(me10)
    df_me10.to_csv(r'C:/Users/Lenovo/Documents/TUGAS AKHIR/Dataset/Data 12 Jam/Tuning TC 10, Eksternal IDV 11, dan IDV 4/me10.csv', header=False, index=False)
    df_me11 = pd.DataFrame(me11)
    df_me11.to_csv(r'C:/Users/Lenovo/Documents/TUGAS AKHIR/Dataset/Data 12 Jam/Tuning TC 10, Eksternal IDV 11, dan IDV 4/me11.csv', header=False, index=False)

if __name__=='__main__': main()








