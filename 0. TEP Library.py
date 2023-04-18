
import math
import numpy as np

class teprob:
    def __init__(self):
        self.cpflmx = self.cpprmx = float(0)
        self.vtr = self.vts = self.vtc = self.vtv = float(0)
        self.hwr = self.hws = float(0)
        
        self.G = float()

        self.xmeas = np.zeros(41)
        self.xmv = np.zeros(12)
        self.idv = np.zeros((20,), dtype=int)
        self.uclr = np.zeros(8)
        self.ucvr = np.zeros(8)
        self.xlr = np.zeros(8)
        self.xvr = np.zeros(8)
        self.ppr = np.zeros(8)
        self.crxr = np.zeros(8)
        self.ucls = np.zeros(8)
        self.ucvs = np.zeros(8)
        self.xls = np.zeros(8)
        self.xvs = np.zeros(8)
        self.pps = np.zeros(8)
        self.uclc = np.zeros(8)
        self.xlc = np.zeros(8)
        self.ucvv = np.zeros(8)
        self.xvv = np.zeros(8)
        self.sfr = np.zeros(8)
        self.rr = np.zeros(4)
        self.vcv = np.zeros(12)      #
        self.vrng = np.zeros(12)     #
        self.vtau = np.zeros(12)     #
        self.ftm = np.zeros(13)
        self.fcm = np.zeros((8,13))
        self.xst = np.zeros((8,13))
        self.xmws = np.zeros(13)
        self.hst = np.zeros(13)
        self.tst = np.zeros(13)
        self.htr = np.zeros(2)
        self.xdel = np.zeros(41)
        self.xns = np.zeros(41)
        self.vst = np.zeros(12)
        self.ivst = np.zeros((12,), dtype=int)
        #yy = np.zeros(50)

        self.adist = np.zeros(12)
        self.bdist = np.zeros(12)
        self.cdist = np.zeros(12)
        self.ddist = np.zeros(12)
        self.tlast = np.zeros(12)
        self.tnext = np.zeros(12)
        self.hspan = np.zeros(12)
        self.hzero = np.zeros(12)
        self.szero = np.zeros(12)
        self.sspan = np.zeros(12)
        self.spspan = np.zeros(12)
        self.idvwlk = np.zeros((12,), dtype=int)

        # konstanta
        self.avp = np.zeros(8)
        self.bvp = np.zeros(8)
        self.cvp = np.zeros(8)
        self.ah = np.zeros(8)
        self.bh = np.zeros(8)
        self.ch = np.zeros(8)
        self.ag = np.zeros(8)
        self.bg = np.zeros(8)
        self.cg = np.zeros(8)
        self.av = np.zeros(8)
        self.ad = np.zeros(8)
        self.bd = np.zeros(8)
        self.cd = np.zeros(8)
        self.xmw = np.zeros(8)

        self.fin = np.zeros(8)
        self.vpos = np.zeros(12)
        self.xcmp = np.zeros(41)


    def tefunc(self,nn, Time, yy, yp, xmv , idv):
        tcr = tcs = tcc = tcv = float(0)
        dlr = dls = dlc = float(0)
        xmns = tgas = tprod = float(0)

        # Limiting of disturbace activations
        for j in range(0,20):
            if (idv[j] > 0):
                idv[j]=1
            else:
                idv[j]=0
        
        # Assignment of disturbance activations
        self.idvwlk[0] = idv[7]
        self.idvwlk[1] = idv[7]
        self.idvwlk[2] = idv[8]
        self.idvwlk[3] = idv[9]
        self.idvwlk[4] = idv[10]
        self.idvwlk[5] = idv[11]
        self.idvwlk[6] = idv[12]
        self.idvwlk[7] = idv[12]
        self.idvwlk[8] = idv[15]
        self.idvwlk[9] = idv[16]
        self.idvwlk[10] = idv[17]
        self.idvwlk[11] = idv[19]
        
        # Recalcultion of disturbace process parameters - determination of
        # process to be updated (1-9/ 10-12)
        for j in range(0,9):
            if (Time >= self.tnext[j]):
                hwlk = self.tnext[j] - self.tlast[j]
                swlk = self.adist[j] + hwlk*(self.bdist[j] + hwlk*(self.cdist[j] + hwlk*self.ddist[j]))
                spwlk = self.bdist[j] + hwlk*(2.0*self.cdist[j] + 3.0*hwlk*self.ddist[j])
                self.tlast[j] = self.tnext[j]
                self.adist[j], self.bdist[j], self.cdist[j], self.ddist[j], self.tnext[j] = self.tesub5(swlk,spwlk,self.adist[j],self.bdist[j],self.cdist[j],
                          self.ddist[j],self.tlast[j],self.tnext[j],self.hspan[j],self.hzero[j],
                          self.sspan[j],self.szero[j],self.spspan[j],self.idvwlk[j])
                          
        for j in range (9,12):  #
            if(Time >= self.tnext[j]):
                hwlk = self.tnext[j] - self.tlast[j]
                swlk = self.adist[j] + hwlk*(self.bdist[j] + hwlk*(self.cdist[j] + hwlk*self.ddist[j]))
                spwlk = self.bdist[j] + hwlk*(2.0*self.cdist[j] + 3.0*hwlk*self.ddist[j])
                self.tlast[j] = self.tnext[j]
                if(swlk > 0.10):
                    self.adist[j]=swlk
                    self.bdist[j]=spwlk
                    self.cdist[j]= -(3.0*swlk + 0.20*spwlk)/0.010
                    self.ddist[j]=(2.0*swlk + 0.10*spwlk)/0.0010
                    self.tnext[j]=self.tlast[j]+0.10
                else:
                    isd= -1
                    hwlk = self.hspan[j]*self.tesub7(isd) + self.hzero[j]
                    self.adist[j]=0.0
                    self.bdist[j]=0.0
                    self.cdist[j]=np.double(self.idvwlk[j])/(hwlk**2)
                    self.ddist[j]=0.0
                    self.tnext[j]=self.tlast[j] + hwlk
        
        # Initialization of disturbace processes parameters
        if(Time == 0.0):
            for j in range(1-1,12):
                self.adist[j]=self.szero[j]
                self.bdist[j]=0.0
                self.cdist[j]=0.0
                self.ddist[j]=0.0
                self.tlast[j]=0.00
                self.tnext[j]=0.10
        
        # Determination of disturbed values
        self.xst[0,3] = self.tesub8(1,Time) - idv[0]*0.030 - idv[1]*2.43719e-3    #
        self.xst[1,3] = self.tesub8(2,Time) + idv[1]*0.0050                       #
        self.xst[2,3] = 1.0- self.xst[0,3] - self.xst[1,3]
        self.tst[0] = self.tesub8(3,Time) + idv[2]*5.0
        self.tst[3] = self.tesub8(4,Time)
        tcwr = self.tesub8(5,Time) + idv[3]*5.0
        tcws = self.tesub8(6,Time) + idv[4]*5.0
        r1f = self.tesub8(7,Time)
        r2f = self.tesub8(8,Time)
        
        #Retrieving of current states
        for j in range(0,3):
            self.ucvr[j] = yy[j]
            self.ucvs[j] = yy[j+9]
            self.uclr[j] = 0.0
            self.ucls[j] = 0.0

        for j in range(3,8):
            self.uclr[j] = yy[j]
            self.ucls[j] = yy[j+9]

        for j in range(0,8):
            self.uclc[j] = yy[j+18]
            self.ucvv[j] = yy[j+27]

        etr = yy[9-1]
        ets = yy[18-1]
        etc = yy[27-1]
        etv = yy[36-1]
        twr = yy[37-1]
        tws = yy[38-1]
        for j in range(0,12):
            self.vpos[j] = yy[j+38]  
        
        #Calculation of collective holdup
        utlr = 0.0
        utls = 0.0
        utlc = 0.0
        utvv = 0.0
        for j in range(0,8):
            utlr = utlr + self.uclr[j]
            utls = utls + self.ucls[j]
            utlc = utlc + self.uclc[j]
            utvv = utvv + self.ucvv[j]
        
        # Calculation of component concentration
        for j in range(0,8):
            self.xlr[j] = self.uclr[j]/utlr
            self.xls[j] = self.ucls[j]/utls
            self.xlc[j] = self.uclc[j]/utlc
            self.xvv[j] = self.ucvv[j]/utvv
        
        #Calculation of specific internal energy
        esr = etr/utlr
        ess = ets/utls
        esc = etc/utlc
        esv = etv/utvv
        
        #Calculation of tempretures
        a = tcr
        tcr = self.tesub2(self.xlr,a,esr,0)
        tkr = tcr + 273.15

        b = tcs
        tcs = self.tesub2(self.xls,b ,ess,0)
        tks = tcs + 273.15
        
        c = tcc
        tcc = self.tesub2(self.xlc,c,esc,0)
        
        d = tcv
        tcv = self.tesub2(self.xvv,d,esv,2)
        tkv = tcv + 273.15
        
        # Calculation of densities
        dlr = self.tesub4(self.xlr,tcr,dlr)
        dls = self.tesub4(self.xls,tcs,dls)
        dlc = self.tesub4(self.xlc,tcc,dlc)
        
        #Calculation of volume of liquid and vapor phase
        vlr = utlr/dlr
        vls = utls/dls
        vlc = utlc/dlc
        vvr = self.vtr - vlr
        vvs = self.vts - vls
        
        # Calculation of pressure
        rg = 998.9
        ptr = 0.0
        pts = 0.0
        for j in range(0,3):
            self.ppr[j] = self.ucvr[j]*rg*tkr/vvr
            ptr = ptr + self.ppr[j]
            self.pps[j] = self.ucvs[j]*rg*tks/vvs
            pts = pts + self.pps[j]

        for j in range(3,8):
            vpr = np.exp(self.avp[j] + self.bvp[j]/(tcr+self.cvp[j]))
            self.ppr[j] = vpr*self.xlr[j]
            ptr = ptr + self.ppr[j]
            vpr = np.exp(self.avp[j] + self.bvp[j]/(tcs + self.cvp[j]))
            self.pps[j] = vpr*self.xls[j]
            pts = pts + self.pps[j]

        ptv = utvv*rg*tkv/self.vtv
        
        #Calculation of component concentration in vapor phase (reactor and
        # separator)
        for j in range(0,8):
            self.xvr[j] = self.ppr[j]/ptr
            self.xvs[j] = self.pps[j]/pts
        
        #Calculation of collective holdup of component in vapor phase (reactor
        # and separator)
        utvr = ptr*vvr/rg/tkr
        utvs = pts*vvs/rg/tks
        
        # Calculation of single holdup of components in vapor phase (reactor
        #  and separator)
        for j in range(3,8):
            self.ucvr[j] = utvr*self.xvr[j]
            self.ucvs[j] = utvs*self.xvs[j]
        
        # Reaction kinetics
        self.rr[0] = np.exp(31.5859536 - 40000.0/1.987 / tkr)*r1f
        self.rr[1] = np.exp(3.00094014-20000.0/1.987 / tkr)*r2f
        self.rr[2] = np.exp(53.4060443-60000.0/1.987 / tkr)
        self.rr[3] = self.rr[2]*0.7674883340

        if (self.ppr[0]>0.0 and self.ppr[2]>0.0):
            r1f = self.ppr[0]**1.1544
            r2f = self.ppr[2]**0.3735
            self.rr[0] = self.rr[0]*r1f*r2f*self.ppr[3]
            self.rr[1] = self.rr[1]*r1f*r2f*self.ppr[4]
        else:
            self.rr[0] = 0.0
            self.rr[1] = 0.0
            
        self.rr[2] = self.rr[2]*self.ppr[0]*self.ppr[4]
        self.rr[3] = self.rr[3]*self.ppr[0]*self.ppr[3]
        
        for j in range(0,4):
            self.rr[j] = self.rr[j]*vvr
        
        #Consumption and creation of components in reactor
        self.crxr[0] = -self.rr[0] - self.rr[1] - self.rr[2]
        self.crxr[2] = -self.rr[0] - self.rr[1]
        self.crxr[3] = -self.rr[0] - 1.50*self.rr[3]
        self.crxr[4] = -self.rr[1] - self.rr[2]
        self.crxr[5] = self.rr[2] + self.rr[3]
        self.crxr[6] = self.rr[0]
        self.crxr[7] = self.rr[1]
        rh = self.rr[0]*self.htr[0] + self.rr[1]*self.htr[1]
        self.xmws[0] = 0.0
        self.xmws[1] = 0.0
        self.xmws[5] = 0.0
        self.xmws[7] = 0.0
        self.xmws[8] = 0.0
        self.xmws[9] = 0.0
        for j in range (0,8):
            self.xst[j,6-1] = self.xvv[j]
            self.xst[j,8-1] = self.xvr[j]
            self.xst[j,9-1] = self.xvs[j]
            self.xst[j,10-1] = self.xvs[j]
            self.xst[j,11-1] = self.xls[j]
            self.xst[j,13-1] = self.xlc[j]
            self.xmws[1-1] = self.xmws[1-1] + self.xst[j,1-1]*self.xmw[j]
            self.xmws[2-1] = self.xmws[2-1] + self.xst[j,2-1]*self.xmw[j]
            self.xmws[6-1] = self.xmws[6-1] + self.xst[j,6-1]*self.xmw[j]
            self.xmws[8-1] = self.xmws[8-1] + self.xst[j,8-1]*self.xmw[j]
            self.xmws[9-1] = self.xmws[9-1] + self.xst[j,9-1]*self.xmw[j]
            self.xmws[10-1] = self.xmws[10-1] +self.xst[j,10-1]*self.xmw[j]
            
        self.tst[5] = tcv
        self.tst[7] = tcr
        self.tst[8] = tcs
        self.tst[9] = tcs
        self.tst[10] = tcs
        self.tst[12] = tcc
        self.hst[0] = self.tesub1(self.xst[:,0],self.tst[0],self.hst[0],1)
        self.hst[1] = self.tesub1(self.xst[:,1],self.tst[1],self.hst[1],1)
        self.hst[2] = self.tesub1(self.xst[:,2],self.tst[2],self.hst[2],1)
        self.hst[3] = self.tesub1(self.xst[:,3],self.tst[3],self.hst[3],1)
        self.hst[5] = self.tesub1(self.xst[:,5],self.tst[5],self.hst[5],1)
        self.hst[7] = self.tesub1(self.xst[:,7],self.tst[7],self.hst[7],1)
        self.hst[8] = self.tesub1(self.xst[:,8],self.tst[8],self.hst[8],1)
        self.hst[9] = self.hst[8]
        self.hst[10] = self.tesub1(self.xst[:,10],self.tst[10],self.hst[10],0)
        self.hst[12] = self.tesub1(self.xst[:,12],self.tst[12],self.hst[12],0)

        self.ftm[0] = self.vpos[0]*self.vrng[0]/100.0
        self.ftm[1] = self.vpos[1]*self.vrng[1]/100.0
        self.ftm[2] = self.vpos[2]*(1.0-idv[5])*self.vrng[2]/100.0
        self.ftm[3] = self.vpos[3]*(1.0 - (idv[6]*0.20))*(self.vrng[3]/100.0) + (1.0e-10)
        self.ftm[10] = self.vpos[6]*self.vrng[6]/100.0
        self.ftm[12] = self.vpos[7]*self.vrng[7]/100.0
        uac = self.vpos[8]*self.vrng[8]*(1.0 + self.tesub8(9,Time))/100.0
        fwr = self.vpos[9]*self.vrng[9]/100.0
        fws = self.vpos[10]*self.vrng[10]/100.0
        agsp = (self.vpos[11]+150.0)/100.0
        dlp = ptv - ptr
        if (dlp < 0.0):
            dlp=0.0
        flms = 1937.60*math.sqrt(dlp)
        self.ftm[5] = flms/self.xmws[5]
        
        dlp = ptr-pts
        if (dlp < 0.0):
            dlp = 0.0
        flms = 4574.210*math.sqrt(dlp)*(1.0 - 0.250*self.tesub8(12,Time))
        self.ftm[7] = flms/self.xmws[7]
        
        dlp = pts - 760.0
        if (dlp < 0.0):
            dlp=0.0
        flms = self.vpos[5]*0.1511690*math.sqrt(dlp)
        self.ftm[9] = flms/self.xmws[9]
        self.pr = ptv/pts
        if (self.pr < 1.0):
            self.pr = 1.0
        if (self.pr > self.cpprmx):
            self.pr = self.cpprmx
        flcoef = self.cpflmx/1.1970
        flms = self.cpflmx + flcoef*(1.0 - self.pr**3)
        cpdh = flms*(tcs + 273.150)*(1.8e-6)*1.98720*(ptv - pts)/(self.xmws[8]*pts)
        dlp = ptv-pts
        if (dlp < 0.0):
            dlp=0.0
        flms = flms-self.vpos[4]*53.3490*math.sqrt(dlp)
        if (flms < 1.0e-3):
            flms = 1.0e-3
        self.ftm[8] = flms/self.xmws[8]
        self.hst[8] = self.hst[8] + cpdh/self.ftm[8]

        for j in range(0,8):
            self.fcm[j,0] = self.xst[j,0]*self.ftm[0]
            self.fcm[j,1] = self.xst[j,1]*self.ftm[1]
            self.fcm[j,2] = self.xst[j,2]*self.ftm[2]
            self.fcm[j,3] = self.xst[j,3]*self.ftm[3]
            self.fcm[j,5] = self.xst[j,5]*self.ftm[5]
            self.fcm[j,7] = self.xst[j,7]*self.ftm[7]
            self.fcm[j,8] = self.xst[j,8]*self.ftm[8]
            self.fcm[j,9] = self.xst[j,9]*self.ftm[9]
            self.fcm[j,10] = self.xst[j,10]*self.ftm[10]
            self.fcm[j,12] = self.xst[j,12]*self.ftm[12]

        if (self.ftm[10] > 0.1):
            if (tcc > 170.0):
                tmpfac = tcc - 120.262
            elif (tcc < 5.292):
                tmpfac=0.1
            else:
                tmpfac = 363.744/(177.0 - tcc) - 2.22579488

            vovrl = self.ftm[3]/self.ftm[10]*tmpfac
            self.sfr[3] = 8.5010*vovrl/(1.0 + 8.5010*vovrl)
            self.sfr[4] = 11.402*vovrl/(1.0 + 11.402*vovrl)
            self.sfr[5] = 11.795*vovrl/(1.0 + 11.795*vovrl)
            self.sfr[6] = 0.0480*vovrl/(1.0 + 0.0480*vovrl)
            self.sfr[7] = 0.0242*vovrl/(1.0 + 0.0242*vovrl)
        else:
            self.sfr[3]=0.9999
            self.sfr[4]=0.999
            self.sfr[5]=0.999
            self.sfr[6]=0.99
            self.sfr[7]=0.98

        for j in range(0,8):
            self.fin[j] = 0.0                    #
            self.fin[j] = self.fin[j] + self.fcm[j,3]
            self.fin[j] = self.fin[j] + self.fcm[j,10]

        self.ftm[4] = 0.0
        self.ftm[11] = 0.0
        for j in range(0,8):
            self.fcm[j,4] = self.sfr[j]*self.fin[j]
            self.fcm[j,11] = self.fin[j] - self.fcm[j, 4]
            self.ftm[4] = self.ftm[4] + self.fcm[j,4]
            self.ftm[11] = self.ftm[11] + self.fcm[j,11] 

        for j in range(0,8):
            self.xst[j,4] = self.fcm[j,4]/self.ftm[4]
            self.xst[j,11] = self.fcm[j,11]/self.ftm[11]

        self.tst[4] = tcc
        self.tst[11] = tcc
        self.hst[4] = self.tesub1(self.xst[:,4],self.tst[4],self.hst[4],1)
        self.hst[11] = self.tesub1(self.xst[:,11],self.tst[11],self.hst[11],0)
        self.ftm[6] = self.ftm[5]
        self.hst[6] = self.hst[5]
        self.tst[6] = self.tst[5]
        for j in range(0,8):
            self.xst[j,6]=self.xst[j,5]
            self.fcm[j,6]=self.fcm[j,5]
        
        #Calculation of heat transfer in reactor
        if (vlr/7.8 > 50.0):
            uarlev = 1.0
        elif (vlr/7.8 < 10.0):
            uarlev = 0.0
        else:
            uarlev=0.025*vlr/7.8 - 0.25
        
        uar = uarlev*(-0.5*(agsp**2) + 2.75*agsp - 2.5)*855490.0e-6
        qur = uar*(twr - tcr)*(1.0 - 0.350*self.tesub8(10,Time))
        
        #Calculation of heat transfer in condenser(separator)
        uas = 0.404655*(1.0 - 1.0/(1.0 + (self.ftm[7]/3528.73)**4))
        qus = uas*(tws-self.tst[7])*(1.0 - 0.250*self.tesub8(11,Time))
        
        #Calculation of heat transfer in stripper
        quc = 0.0
        if (tcc < 100.0):
            quc = uac*(100.0-tcc)
        
        #Settig of measured values
        self.xmeas[0] = self.ftm[2]*0.359/35.3145         #
        self.xmeas[1] = self.ftm[0]*self.xmws[0]*0.454
        self.xmeas[2] = self.ftm[1]*self.xmws[1]*0.454
        self.xmeas[3] = self.ftm[3]*(0.359/35.3145)
        self.xmeas[4] = self.ftm[8]*0.359/35.3145
        self.xmeas[5] = self.ftm[5]*0.359/35.3145
        self.xmeas[6] = (ptr-760.0)/760.0*101.325
        self.xmeas[7] = (vlr-84.6)/666.7*100.0
        self.xmeas[8] = tcr
        self.xmeas[9] = self.ftm[9]*0.359/35.3145
        self.xmeas[10] = tcs
        self.xmeas[11] = (vls-27.5)/290.0*100.0
        self.xmeas[12] = (pts-760.0)/760.0*101.325
        self.xmeas[13] = self.ftm[10]/dls/35.3145
        self.xmeas[14] = (vlc-78.25)/self.vtc*100.0
        self.xmeas[15] = (ptv-760.0)/760.0*101.325
        self.xmeas[16] = self.ftm[12]/dlc/35.3145
        self.xmeas[17] = tcc
        self.xmeas[18] = quc*1.04e3 * 0.454
        self.xmeas[19] = cpdh*0.0003927e6
        self.xmeas[19] = cpdh*0.29307e3
        self.xmeas[20] = twr
        self.xmeas[21] = tws
        
        #Checking of shutdown-constraints
        isd = 0
        if (self.xmeas[6] > 3000.0):
            isd = 1
            print("High Reactor Pressure!!  Shutting down.")
        if (vlr/35.3145 > 24.0):
            isd = 2
            print("High Reactor Liquid Level!!  Shutting down.")
        if (vlr/35.3145 < 2.0):
            isd = 3
            print("Low Reactor Liquid Level!!  Shutting down.")
        if (self.xmeas[8] > 175.0):
            isd=4
            print("High Reactor Temperature!!  Shutting down.")
        if (vls/35.3145 > 12.0):
            isd=5
            print("High Separator Liquid Level!!  Shutting down.")
        if (vls/35.3145 < 1.0):
            isd=6
            print("Low Separator Liquid Level!!  Shutting down.")
        if (vlc/35.3145 > 8.0):
            isd=7
            print("High Stripper Liquid Level!!  Shutting down.")
        if (vlc/35.3145 < 1.0):
            isd=8
            print("Low Stripper Liquid Level!!  Shutting down.")
        if (Time>0.0 and isd==0):
            for j in range(0,22):
                xmns = self.tesub6(self.xns[j], xmns)
                self.xmeas[j] = self.xmeas[j] + xmns
        
        #Analyzer outputs
        self.xcmp[22] = self.xst[0,6]*100.0
        self.xcmp[23] = self.xst[1,6]*100.0
        self.xcmp[24] = self.xst[2,6]*100.0
        self.xcmp[25] = self.xst[3,6]*100.0
        self.xcmp[26] = self.xst[4,6]*100.0
        self.xcmp[27] = self.xst[5,6]*100.0
        self.xcmp[28] = self.xst[0,9]*100.0
        self.xcmp[29] = self.xst[1,9]*100.0
        self.xcmp[30] = self.xst[2,9]*100.0
        self.xcmp[31] = self.xst[3,9]*100.0
        self.xcmp[32] = self.xst[4,9]*100.0
        self.xcmp[33] = self.xst[5,9]*100.0
        self.xcmp[34] = self.xst[6,9]*100.0
        self.xcmp[35] = self.xst[7,9]*100.0
        self.xcmp[36] = self.xst[3,12]*100.0
        self.xcmp[37] = self.xst[4,12]*100.0
        self.xcmp[38] = self.xst[5,12]*100.0
        self.xcmp[39] = self.xst[6,12]*100.0
        self.xcmp[40] = self.xst[7,12]*100.0
        
        ## dari N. L. Ricker
        #self.xmeas[49] = self.xcmp[39]
        #self.xmeas[50] = self.xcmp[40]
        
        if (Time==0.0):
            for j in range(22,41):
                self.xdel[j] = self.xcmp[j]
                self.xmeas[j] = self.xcmp[j]

            tgas=0.1
            tprod=0.25

        if (Time >= tgas):
            for j in range(22,36):
                self.xmeas[j] = self.xdel[j]
                xmns = self.tesub6(self.xns[j],xmns)
                self.xmeas[j] = self.xmeas[j] + xmns
                self.xdel[j] = self.xcmp[j]
            tgas = tgas+0.1
            
        if (Time >= tprod):
            for j in range(36,41):
                self.xmeas[j] = self.xdel[j]
                xmns = self.tesub6(self.xns[j],xmns)
                self.xmeas[j] = self.xmeas[j] + xmns
                self.xdel[j] = self.xcmp[j]
            tprod = tprod+0.25
        
        #Calculation of state derivative
        for j in range(0,8):
            yp[j] = self.fcm[j,6] - self.fcm[j,7] + self.crxr[j]
            yp[j+9] = self.fcm[j,7] - self.fcm[j,8] - self.fcm[j,9] - self.fcm[j,10]
            yp[j+18] = self.fcm[j,11] - self.fcm[j,12]
            yp[j+27] = self.fcm[j,0] + self.fcm[j,1] + self.fcm[j,2] + \
                        self.fcm[j,4] + self.fcm[j,8] - self.fcm[j,5]

        yp[8] = self.hst[6]*self.ftm[6] - self.hst[7]*self.ftm[7] + rh + qur
        # Here is the "correct" version of the separator energy balance: 
        # .YP(18)=HST(8)*FTM(8)-
        # .(HST(9)*FTM(9)-cpdh)-
        # .HST(10)*FTM(10)-
        # .HST(11)*FTM(11)+s
        # .QUS
        
        #Here is the original version
        yp[17] = self.hst[7]*self.ftm[7] - self.hst[8]*self.ftm[8] - self.hst[9]*self.ftm[9] - \
                self.hst[10]*self.ftm[10] + qus
        yp[26] = self.hst[3]*self.ftm[3] + self.hst[10]*self.ftm[10] - \
                self.hst[4]*self.ftm[4] - self.hst[12]*self.ftm[12] + quc
        yp[35] = self.hst[0]*self.ftm[0] + self.hst[1]*self.ftm[1] + \
                self.hst[2]*self.ftm[2] + self.hst[4]*self.ftm[4] + \
                self.hst[8]*self.ftm[8] - self.hst[5]*self.ftm[5]
        yp[36] = (fwr*500.53*(tcwr-twr) - qur*1.0e6/1.8)/self.hwr  #
        yp[37] = (fws*500.53*(tcws-tws) - qus*1.0e6/1.8)/self.hws  #
        self.ivst[9] = idv[13]
        self.ivst[10] = idv[14]
        self.ivst[4] = idv[18]
        self.ivst[6] = idv[18]
        self.ivst[7] = idv[18]
        self.ivst[8] = idv[18]
        for j in range(0,12):
            if (Time==0.0 or (abs(self.vcv[j] - xmv[j]) > self.vst[j]*self.ivst[j])):
                self.vcv[j] = xmv[j]
            
            #Constraints of manipulated variables
            if (self.vcv[j] < 0.0):
                self.vcv[j] = 0.0
            if (self.vcv[j] > 100.0):
                self.vcv[j] = 100.0
            yp[j+38] = (self.vcv[j] - self.vpos[j])/self.vtau[j]

        if ((Time>0.0) and (isd != 0)):
            for j in range(0,nn):
                yp[j]=0.0
        
        self.xmv = xmv
        self.idv = idv
        return yp


    def teinit(self, nn, Time, yy, yp):

    #C
    #C       Initialization
    #C
    #C         Inputs:
    #C
    #C           nn   = Number of differential equations
    #C
    #C         Outputs:
    #C
    #C           Time = Current Time(hrs)
    #C           yy   = Current state values
    #C           yp   = Current derivative values

        # Component data
        self.xmw[0]=2.0
        self.xmw[1]=25.4
        self.xmw[2]=28.0
        self.xmw[3]=32.0
        self.xmw[4]=46.0
        self.xmw[5]=48.0
        self.xmw[6]=62.0
        self.xmw[7]=76.0
        
        self.avp[0]=0.0
        self.avp[1]=0.0
        self.avp[2]=0.0
        self.avp[3]=15.92
        self.avp[4]=16.35
        self.avp[5]=16.35
        self.avp[6]=16.43
        self.avp[7]=17.21
        self.bvp[0]=0.0
        self.bvp[1]=0.0
        self.bvp[2]=0.0
        self.bvp[3]=-1444.0
        self.bvp[4]=-2114.0
        self.bvp[5]=-2114.0
        self.bvp[6]=-2748.0
        self.bvp[7]=-3318.0
        self.cvp[0]=0.0
        self.cvp[1]=0.0
        self.cvp[2]=0.0
        self.cvp[3]=259.0
        self.cvp[4]=265.5
        self.cvp[5]=265.5
        self.cvp[6]=232.9
        self.cvp[7]=249.6
        
        self.ad[0]=1.0
        self.ad[1]=1.0
        self.ad[2]=1.0
        self.ad[3]=23.3
        self.ad[4]=33.9
        self.ad[5]=32.8
        self.ad[6]=49.9
        self.ad[7]=50.5
        self.bd[0]=0.0
        self.bd[1]=0.0
        self.bd[2]=0.0
        self.bd[3]=-0.0700
        self.bd[4]=-0.0957
        self.bd[5]=-0.0995
        self.bd[6]=-0.0191
        self.bd[7]=-0.0541
        self.cd[0]=0.0
        self.cd[1]=0.0
        self.cd[2]=0.0
        self.cd[3]=-0.0002
        self.cd[4]=-0.000152
        self.cd[5]=-0.000233
        self.cd[6]=-0.000425
        self.cd[7]=-0.000150
        
        self.ah[0]=1.0e-6
        self.ah[1]=1.0e-6
        self.ah[2]=1.0e-6
        self.ah[3]=0.960e-6
        self.ah[4]=0.573e-6
        self.ah[5]=0.652e-6
        self.ah[6]=0.515e-6
        self.ah[7]=0.471e-6
        self.bh[0]=0.0
        self.bh[1]=0.0
        self.bh[2]=0.0
        self.bh[3]=8.70e-9
        self.bh[4]=2.41e-9
        self.bh[5]=2.18e-9
        self.bh[6]=5.65e-10
        self.bh[7]=8.70e-10
        self.ch[0]=0.0
        self.ch[1]=0.0
        self.ch[2]=0.0
        self.ch[3]=4.81e-11
        self.ch[4]=1.82e-11
        self.ch[5]=1.94e-11
        self.ch[6]=3.82e-12
        self.ch[7]=2.62e-12
        
        self.av[0]=1.0e-6
        self.av[1]=1.0e-6
        self.av[2]=1.0e-6
        self.av[3]=86.7e-6
        self.av[4]=160.e-6
        self.av[5]=160.e-6
        self.av[6]=225.e-6
        self.av[7]=209.e-6
        self.ag[0]=3.411e-6
        self.ag[1]=0.3799e-6
        self.ag[2]=0.2491e-6
        self.ag[3]=0.3567e-6
        self.ag[4]=0.3463e-6
        self.ag[5]=0.3930e-6
        self.ag[6]=0.170e-6
        self.ag[7]=0.150e-6
        self.bg[0]=7.18e-10
        self.bg[1]=1.08e-9
        self.bg[2]=1.36e-11
        self.bg[3]=8.51e-10
        self.bg[4]=8.96e-10
        self.bg[5]=1.02e-9
        self.bg[6]=0.0
        self.bg[7]=0.0
        self.cg[0]=6.0e-13
        self.cg[1]=-3.98e-13
        self.cg[2]=-3.93e-14
        self.cg[3]=-3.12e-13
        self.cg[4]=-3.27e-13
        self.cg[5]=-3.12e-13
        self.cg[6]=0.0
        self.cg[7]=0.0
        
        # Initial state of process (Mode 1)
        yy[0]=10.40491389
        yy[1]=4.363996017
        yy[2]=7.570059737
        yy[3]=0.4230042431
        yy[4]=24.15513437
        yy[5]=2.942597645
        yy[6]=154.3770655
        yy[7]=159.1865960
        yy[8]=2.808522723
        yy[9]=63.75581199
        yy[10]=26.74026066
        yy[11]=46.38532432
        yy[12]=0.2464521543
        yy[13]=15.20484404
        yy[14]=1.852266172
        yy[15]=52.44639459
        yy[16]=41.20394008
        yy[17]=0.5699317760
        yy[18]=0.4306056376
        yy[19]=7.9906200783e-3
        yy[20]=0.9056036089
        yy[21]=1.6054258216e-2
        yy[22]=0.7509759687
        yy[23]=8.8582855955e-2
        yy[24]=48.27726193
        yy[25]=39.38459028
        yy[26]=0.3755297257
        yy[27]=107.7562698
        yy[28]=29.77250546          # false yy[25]
        yy[29]=88.32481135
        yy[30]=23.03929507
        yy[31]=62.85848794
        yy[32]=5.546318688
        yy[33]=11.92244772
        yy[34]=5.555448243
        yy[35]=0.9218489762
        # Cooling water outlet temperatures
        yy[36]=94.59927549
        yy[37]=77.29698353
        # Valve position
        yy[38]=63.05263039
        yy[39]=53.97970677
        yy[40]=24.64355755
        yy[41]=61.30192144
        yy[42]=22.21000000
        yy[43]=40.06374673
        yy[44]=38.10034370
        yy[45]=46.53415582
        yy[46]=47.44573456
        yy[47]=41.10581288
        yy[48]=18.11349055
        yy[49]=50.00000000
        
        for i in range(0,12):
            self.xmv[i] = yy[i+38]
            self.vcv[i] = self.xmv[i]
            self.vst[i] = 2.00
            self.ivst[i] = 0
        
        # Nominal flowrate through valve
        self.vrng[0]=400.00
        self.vrng[1]=400.00
        self.vrng[2]=100.00
        self.vrng[3]=1500.00
        self.vrng[6]=1500.00
        self.vrng[7]=1000.00
        self.vrng[8]=0.03
        self.vrng[9]=1000.
        self.vrng[10]=1200.0
        
        # Volume of Vessels
        self.vtr=1300.0
        self.vts=3500.0
        self.vtc=156.5
        self.vtv=5000.0
        
        self.htr[0]=0.068993810540
        self.htr[1]=0.050
        self.hwr=7060.
        self.hws=11138.
        self.sfr[0]=0.99500
        self.sfr[1]=0.99100
        self.sfr[2]=0.99000
        self.sfr[3]=0.91600
        self.sfr[4]=0.93600
        self.sfr[5]=0.93800
        self.sfr[6]=5.80000e-2
        self.sfr[7]=3.01000e-2
        
        # Stream 2
        self.xst[0,0]=0.0
        self.xst[1,0]=0.0001
        self.xst[2,0]=0.0
        self.xst[3,0]=0.9999
        self.xst[4,0]=0.0
        self.xst[5,0]=0.0
        self.xst[6,0]=0.0
        self.xst[7,0]=0.0
        self.tst[0]=45.
        
        # Stream 3
        self.xst[0,1]=0.0
        self.xst[1,1]=0.0
        self.xst[2,1]=0.0
        self.xst[3,1]=0.0
        self.xst[4,1]=0.9999
        self.xst[5,1]=0.0001
        self.xst[6,1]=0.0
        self.xst[7,1]=0.0
        self.tst[1]=45.
        
        # Stream 1
        self.xst[0,2]=0.9999
        self.xst[1,2]=0.0001
        self.xst[2,2]=0.0
        self.xst[3,2]=0.0
        self.xst[4,2]=0.0
        self.xst[5,2]=0.0
        self.xst[6,2]=0.0
        self.xst[7,2]=0.0
        self.tst[2]=45.
        
        # Stream 4
        self.xst[0,3]=0.4850
        self.xst[1,3]=0.0050
        self.xst[2,3]=0.5100
        self.xst[3,3]=0.0
        self.xst[4,3]=0.0
        self.xst[5,3]=0.0
        self.xst[6,3]=0.0
        self.xst[7,3]=0.0
        self.tst[3]=45.
        
        self.cpflmx=280275.
        self.cpprmx=1.3
        
        # Time constant of valves
        self.vtau[0]=8.
        self.vtau[1]=8.
        self.vtau[2]=6.
        self.vtau[3]=9.
        self.vtau[4]=7.
        self.vtau[5]=5.
        self.vtau[6]=5.
        self.vtau[7]=5.
        self.vtau[8]=120.
        self.vtau[9]=5.
        self.vtau[10]=5.
        self.vtau[11]=5.
        for i in range(0,12):
            self.vtau[i] = self.vtau[i]/3600.
        
        # Random number
        self.G=4651207995.0
    #C	d00_tr_new: self.G=5687912315.D0       
    #C      original: self.G=1431655765.D0
    #C        d00_tr: self.G=4243534565.D0
    #C        d01_tr: self.G=7854912354.D0
    #C        d02_tr: self.G=3456432354.D0
    #C        d03_tr: self.G=1731738903.D0
    #C        d04_tr: self.G=4346024432.D0
    #C        d05_tr: self.G=5784921734.D0
    #C        d06_tr: self.G=6678322168.D0
    #C        d07_tr: self.G=7984782901.D0
    #C        d08_tr: self.G=8934302332.D0
    #C        d09_tr: self.G=9873223412.D0
    #C        d10_tr: self.G=1089278833.D0
    #C        d11_tr: self.G=1940284333.D0
    #C        d12_tr: self.G=2589274931.D0
    #C        d13_tr: self.G=3485834345.D0
    #C        d14_tr: self.G=4593493842.D0
    #C        d15_tr: self.G=5683213434.D0
    #C        d16_tr: self.G=6788343442.D0
    #C        d17_tr: self.G=1723234455.D0
    #C        d18_tr: self.G=8943243993.D0
    #C       dd18_tr: self.G=1234567890.D0

    #C        d19_tr: self.G=9445382439.D0
    #C        d20_tr: self.G=9902234324.D0
    #C        d21_tr: self.G=2144342545.D0
    #C        d22_tr: self.G=3433249064.D0
    #C        d23_tr: self.G=4356565463.D0
    #C        d24_tr: self.G=8998485332.D0
    #C        d25_tr: self.G=7654534567.D0
    #C        d26_tr: self.G=5457789234.D0

    #C        d00_te: self.G=1254545354.D0
    #C        d01_te: self.G=2994833239.D0
    #C        d02_te: self.G=2891123453.D0
    #C        d03_te: self.G=3420494299.D0
    #C        d04_te: self.G=4598956239.D0
    #C        d05_te: self.G=5658678765.D0
    #C        d06_te: self.G=6598593453.D0
    #C        d07_te: self.G=7327843434.D0
    #C        d08_te: self.G=8943242344.D0
    #C        d09_te: self.G=9343430004.D0
    #C        d10_te: self.G=1039839281.D0
    #C        d11_te: self.G=1134345551.D0
    #C        d12_te: self.G=2232323236.D0
    #C        d13_te: self.G=3454354353.D0
    #C        d14_te: self.G=4545445883.D0
    #C        d15_te: self.G=5849489384.D0
    #C        d16_te: self.G=6284545932.D0
    #C        d17_te: self.G=4342232344.D0
    #C        d18_te: self.G=5635346588.D0
    #C        d19_te: self.G=9090909232.DO
    #C        d20_te: self.G=8322308324.D0
    #C        d21_te: self.G=2132432423.D0
    #C        d22_te: self.G=5454589923.D0
    #C        d23_te: self.G=6923255678.D0
    #C        d24_te: self.G=8493323434.D0
    #C        d25_te: self.G=9338398429.D0
    #C        d26_te: self.G=1997072199.D0
        
        # Amplitudes of measurement noise
        self.xns[0]=0.00120
        self.xns[1]=18.0000
        self.xns[2]=22.0000
        self.xns[3]=0.05000
        self.xns[4]=0.20000
        self.xns[5]=0.21000
        self.xns[6]=0.30000
        self.xns[7]=0.50000
        self.xns[8]=0.01000
        self.xns[9]=0.00170
        self.xns[10]=0.01000
        self.xns[11]=1.00000
        self.xns[12]=0.30000
        self.xns[13]=0.12500
        self.xns[14]=1.00000
        self.xns[15]=0.30000
        self.xns[16]=0.11500
        self.xns[17]=0.01000
        self.xns[18]=1.15000
        self.xns[19]=0.20000
        self.xns[20]=0.01000
        self.xns[21]=0.01000
        self.xns[22]=0.2500
        self.xns[23]=0.1000
        self.xns[24]=0.2500
        self.xns[25]=0.1000
        self.xns[26]=0.2500
        self.xns[27]=0.0250
        self.xns[28]=0.2500
        self.xns[29]=0.1000
        self.xns[30]=0.2500
        self.xns[31]=0.1000
        self.xns[32]=0.2500
        self.xns[33]=0.0250
        self.xns[34]=0.0500
        self.xns[35]=0.0500
        self.xns[36]=0.0100
        self.xns[37]=0.0100
        self.xns[38]=0.0100
        self.xns[39]=0.5000
        self.xns[40]=0.5000
        
        # inititialization of disturbanve flags
        for i in range(0,20):
            self.idv[i]=0
        
        # Data of disturbance processes
        self.hspan[0]=0.20
        self.hzero[0]=0.50
        self.sspan[0]=0.030
        self.szero[0]=0.4850
        self.spspan[0]=0.0
        
        self.hspan[1]=0.70
        self.hzero[1]=1.00
        self.sspan[1]=.0030
        self.szero[1]=.0050
        self.spspan[1]=0.0
        
        self.hspan[2]=0.250
        self.hzero[2]=0.50
        self.sspan[2]=10.0
        self.szero[2]=45.0
        self.spspan[2]=0.0
        
        self.hspan[3]=0.70
        self.hzero[3]=1.00
        self.sspan[3]=10.0
        self.szero[3]=45.0
        self.spspan[3]=0.0
        
        self.hspan[4]=0.150
        self.hzero[4]=0.250
        self.sspan[4]=10.0
        self.szero[4]=35.0
        self.spspan[4]=0.0
        
        self.hspan[5]=0.150
        self.hzero[5]=0.250
        self.sspan[5]=10.0
        self.szero[5]=40.0
        self.spspan[5]=0.0
        
        self.hspan[6]=1.0
        self.hzero[6]=2.0
        self.sspan[6]=0.250
        self.szero[6]=1.00
        self.spspan[6]=0.0
        
        self.hspan[7]=1.0
        self.hzero[7]=2.0
        self.sspan[7]=0.250
        self.szero[7]=1.00
        self.spspan[7]=0.0
        
        self.hspan[8]=0.40
        self.hzero[8]=0.50
        self.sspan[8]=0.250
        self.szero[8]=0.00
        self.spspan[8]=0.0
        
        self.hspan[9]=1.50
        self.hzero[9]=2.00
        self.sspan[9]=0.00
        self.szero[9]=0.00
        self.spspan[9]=0.0
        
        self.hspan[10]=2.00
        self.hzero[10]=3.00
        self.sspan[10]=0.00
        self.szero[10]=0.00
        self.spspan[10]=0.0
        
        self.hspan[11]=1.50
        self.hzero[11]=2.00
        self.sspan[11]=0.00
        self.szero[11]=0.00
        self.spspan[11]=0.0
        
        # Initialization of disturbance processes parameters
        for i in range(0,12):
            self.tlast[i]=0.0
            self.tnext[i]=0.10
            self.adist[i]=self.szero[i]
            self.bdist[i]=0.0
            self.cdist[i]=0.0
            self.ddist[i]=0.0
        Time=0.0
        xmv_ = self.xmv
        idv_ = self.idv
        yp = self.tefunc(nn, Time,yy, yp, xmv_,idv_)
        #self.yy = yy
        #self.yp = yp
        return yy, yp

    
    def tesub1(self,z,t,h,ity):
        if (ity==0):
            h = 0.00
            for i in range(0,8):
                hi = t*(self.ah[i] + self.bh[i]*t/2.0 + self.ch[i]*(t**2)/3.0)
                hi = 1.80*hi
                h = h + z[i]*self.xmw[i]*hi

        else:
            h = 0.00
            for i in range(0,8):
                hi = t*(self.ag[i] + self.bg[i]*t/2.0 + \
                    self.cg[i]*(t**2)/3.0)
                hi = 1.80*hi
                hi = hi + self.av[i]
                h = h + z[i]*self.xmw[i]*hi
     
        if (ity==2):
            r = 3.576960e-6       #3.576960/1.e6
            h = h - r*(t + 273.15)

        return h


    def tesub2(self,z,t,h,ity):
        tin = t
        htest = dh = float()
        for j in range(0,100):
            htest = self.tesub1(z,t,htest,ity)
            err = htest-h
            dh = self.tesub3(z,t,dh,ity)
            dt = -err/dh
            t = t + dt
            if (abs(dt)<1e-12):   
                return t          
                                
        t = tin                 
        return t              


    def tesub3(self,z,t,dh,ity):
        if (ity==0):
            dh = 0.00
            for i in range(0,8):
                d = t
                dhi = self.ah[i] + self.bh[i]*t + self.ch[i]*(d**2)
                dhi = 1.80*dhi
                dh = dh + z[i]*self.xmw[i]*dhi
                
        else:
            dh = 0.00
            d = t
            for i in range(0,8):
                dhi = self.ag[i] + self.bg[i]*t + self.cg[i]*(d**2)
                dhi = 1.80*dhi
                dh = dh + z[i]*self.xmw[i]*dhi
     
        if (ity==2):
            r = 3.576960/1.e6
            dh = dh - r

        return dh


    def tesub4(self,x,t,r):
        #v = 0.0
        v = float(0)
        for i in range(0,8):
            v = v + x[i]*self.xmw[i]/(self.ad[i] + (self.bd[i] + self.cd[i]*t)*t)

        #r = 1.0/v
        r = float(1)/v
        return r


    def tesub5(self,s,sp,adist,bdist,cdist,ddist,tlast,
            tnext,hspan,hzero,sspan,szero,spspan,idvflag):
        i = -1
        h = hspan*self.tesub7(i)+hzero
        s1 = sspan*self.tesub7(i)*idvflag + szero
        s1p = spspan*self.tesub7(i)*idvflag
        adist = s
        bdist = sp
        cdist = (3.0*(s1-s) - h*(s1p + 2.0*sp))/(h**2)
        ddist = (2.0*(s-s1) + h*(s1p+sp))/(h**3)
        tnext = tlast + h
        return adist, bdist, cdist, ddist, tnext


    def tesub6(self,std,x):
        x = 0.0      # berpengaruh pada self.xmeas[1],self.xmeas[2],self.xmeas[3]
        for i in range(1,12+1):         # akhir pada i = 12
            x = x + self.tesub7(i)
        x = (x - 6.0)*std
        return x


    def tesub7(self,i):

        self.G = self.G*9228907.0
        #self.G = np.remainder(self.G,4294967296.0)
        self.G = self.G - int(self.G/4294967296.0)*4294967296.0
        #self.G = self.G % 4294967296.0
        
        # generation of random numbers for measurement noise
        if (i>=0):
            tesub7 = self.G/4294967296.0
        # generation of random numbers for process disturbances
        if (i<0):
            tesub7 = self.G*2.0/4294967296.0 - 1.0
        return tesub7


    def tesub8(self,i,t):
        h = t-self.tlast[i-1]
        tesub8 = self.adist[i-1] + h*(self.bdist[i-1] + h*(self.cdist[i-1] + h*self.ddist[i-1]))
        return tesub8

     
    def getOutput(self):
        return self.xmeas, self.xmv, self.idv