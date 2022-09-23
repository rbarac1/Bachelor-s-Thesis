#from textwrap import dedent
import numpy as np


c = 299792458   #speed of light


#CERN Linac4
Vmax = 3e5  #max source voltage
omega = 352e6*2*np.pi   #AC frequency

#protons and ions
# Vmax = 5e7
# omega = 300e6

#T=2*pi/omega=1/3e8

#dtc = 1e-12     #time interval at c

def gamma(v):   #Lorentz factor
    #print(v)
    return 1/np.sqrt(1-np.dot(v,v)/c/c)

# def dgamma(v):  #dgamma/dv
#     return 1/np.sqrt(1-v*v/c/c)**3*v/c/c

def En(v,m):    #particle energy
    return gamma(v)*m*c*c

def Ekin(v, m):     #particle kinetic energy in MeV
    return (gamma(v)-1)*m*c*c*(6.242e12)

def v2fromE(E, m):
    return c**2*(1-(m*m*c**4)/(E*E))

def vratio(q, m, d, vold, vnew):
    return (4*np.abs(q)*Vmax/(omega*d*m)+vold*gamma(vold))/(vnew*gamma(vnew))

#particle object class
class particle:
    def __init__(self, q, m, v0, E=[0.,0.,0.], B=[0.,0.,0.], dt=1e-7, dt0=0):
        self.m = m  #mass
        self.q = q  #charge
        self.E = np.array(E)    #electric field
        self.B = np.array(B)    #magnetic field

        self.v = np.array([v0])     #velocity

        self.r = np.array([[0,0,0]])    #position

        self.t = [0]    #time
        self.dtc = 1e-9     #time interval at c
        self.dt = self.dtc    #time step
        self.t0 = 0     #time after the first collimator
        self.dtpos = [1,1,1,1,1]

        #self.segment = 0

        self.starts = 0   #number of passed starts of a linac segment
        self.stops = 0   #number of passed stops of a linac segment
        self.d = 0  #distance of electrodes in the accelerator

        self.En = 0     #total energy of the particle (in joules)

        self.E_MeV = []     #kinetic energy of a particle in MeV

        #self.x = 0  #start of the segment with an electric field


    def E_field(self): #calculates and returns the electric field
        # if t==0:
        #     t=self.t[-1]
        
        self.E = np.array([0,0,np.sign(self.q)*Vmax*2/self.d*np.sin(omega*(self.t[-1]-self.t0))])   #constat frequency source
        return self.E

    def Emax(self):
        return Vmax*2/self.d    #max electric field strength

    def __ac(self, v):      #acceleration
        if self.starts > self.stops:    #acceleration inside the segments
            return 0
        else:       #acceleration between segments
            #self.E = np.array([self.q*Vmax*2/self.d*np.cos(omega*self.t[-1]),self.q*Vmax*2/self.d*np.cos(omega*self.t[-1]),0])
            #self.E = np.array([0,0,np.sign(self.q)*Vmax*2/self.d*np.sin(omega*(self.t[-1]-self.t0))])
            #self.E_field()
            #self.E = np.array([0,0,0])
            F_clas = self.q*(self.E_field()+np.cross(v, self.B))
            #print("AAAAAAAAA")
            acc = 1/self.m/gamma(v)*(F_clas-np.dot(F_clas,v)*v/c/c)
            #print("acc=", acc)
            return acc
        #return (self.q*(self.E+np.cross(v, self.B))/self.m) / (dgamma(v)*v+gamma(v))

    def move(self, met="RK4"):      #numerically move the particle by a step
        #variable timesteps
        if self.v[-1,2]>2e8:
            self.dt = self.dtc
            if self.dtpos[4] != 1:
                self.dtpos[4] = len(self.t)
                
        elif self.v[-1,2]>2e7:
            self.dt = self.dtc*3
            if self.dtpos[3] != 1:
                self.dtpos[3] = len(self.t)
                
        elif self.v[-1,2]>2e6:
            self.dt = self.dtc*10
            if self.dtpos[2] != 1:
                self.dtpos[2] = len(self.t)
                
        elif self.v[-1,2]>2e5:
            self.dt = self.dtc*100
            if self.dtpos[1] != 1:
                self.dtpos[1] = len(self.t)

        elif self.v[-1,2]>2e4:
            self.dt = self.dtc*1000
            if self.dtpos[0] != 1:
                self.dtpos[0] = len(self.t) 
        

        if met=="energy":   #currently only works with the E field in the z-direction 
            if self.En == 0:
                self.En = En(self.v[-1], self.m)
                print("Calculating energy!")

            dEn = self.v[-1,2]*self.dt*self.q*self.E_field()[2]
            self.En += dEn

            #dvz = np.sqrt(v2fromE(self.En, self.m)-self.v[-1,0]**2-self.v[-1,1]**2)-self.v[-1,2]

            dvz = self.m*self.m*c**6/self.En**3/self.v[-1,2]*dEn
            #print(self.En, dEn, self.E[2])
            #print(vz)

            # print("dE = ", dEn)
            # print("dv = ", dvz)

            self.v = np.append(self.v, [self.v[-1]+np.array([0,0,dvz])],0)
            self.r = np.append(self.r, [self.r[-1]+self.v[-1]*self.dt],0)
            self.t.append(self.t[-1]+self.dt)


        elif met=="RK4":      #runge-kutta 4 method
            k1v = self.__ac(self.v[-1])*self.dt
            k1r = self.v[-1]*self.dt
            # if self.starts==self.stops:
            #     print("k1v=",k1v)

            k2v = self.__ac(self.v[-1]+k1v/2)*self.dt
            k2r = (self.v[-1]+k1v/2)*self.dt
            # if self.starts==self.stops:
            #     print("k2v=",k2v)

            k3v = self.__ac(self.v[-1]+k2v/2)*self.dt
            k3r = (self.v[-1]+k2v/2)*self.dt

            k4v = self.__ac(self.v[-1]+k3v)*self.dt
            k4r = (self.v[-1]+k3v)*self.dt

            self.v = np.append(self.v, [self.v[-1]+ (k1v+2*k2v+2*k3v+k4v)/6], 0)    #velocity

            self.r = np.append(self.r, [self.r[-1]+ 1/6*(k1r+2*k2r+2*k3r+k4r)], 0)  #position
                
            self.t += [self.t[-1]+self.dt]  #time
        
        else:       #Euler's method
            #print(self.v)
            a = self.__ac(self.v[-1])
            self.v = np.append(self.v, [self.v[-1]+a*self.dt],0)
            self.r = np.append(self.r, [self.r[-1]+self.v[-1]*self.dt],0)
            self.t += [self.t[-1]+self.dt]

    def set_Efield(self, E):    #set the electric field
        self.E = np.array(E)

    def reset(self):        #reset the motion
        self.v = np.array([self.v[0]])
        self.r = np.array([self.r[0]])
        self.t = [0]

    def remove_steps(self, N):      #remove the previous N steps of the motion
        self.t = self.t[:-N]
        while N>0:
            self.r = np.delete(self.r, -1, 0)
            self.v = np.delete(self.v, -1, 0)
            N -= 1

    # def segment(self, rad, len):
    #     self.rad = rad
    #     self.len = len

    #def motion(self, T, method="RK4", part="collimator"):
     #   if part=="collimator":




    # else:
    #     while(self.t[-1]<T):
    #         #print(self.r)
    #         self.move(method)

        #return self.r



class electron(particle):       #electron subclass
    def __init__(self, v0, E=[0.,0.,0.], B=[0.,0.,0.], dt=1e-7):
        self.m = 9.1093837e-31      #electron mass
        self.q = -1.60217663e-19    #electron charge
        super(electron, self).__init__(self.q, self.m, v0, E, B, dt)

class proton(particle):         #proton subclass
    def __init__(self, v0, E=[0.,0.,0.], B=[0.,0.,0.], dt=1e-7):
        self.m = 1.67262192e-27     #proton mass
        self.q = 1.60217663e-19     #proton charge
        super(proton, self).__init__(self.q, self.m, v0, E, B, dt)

class ion(particle):         #ion subclass
    def __init__(self, v0, C, A, E=[0.,0.,0.], B=[0.,0.,0.], dt=1e-7):
        self.m = A*1.66e-27     #ion mass
        self.q = C*1.60217663e-19     #ion charge
        super(ion, self).__init__(self.q, self.m, v0, E, B, dt)
###############################################################################################################
###############################################################################################################





###############################################################################################################
###############################################################################################################
#accelerator parts
class segment:      #segment class
    def __init__(self, l, R, name):
        self.l = l      #segment length
        self.R = R      #segment radius (cylindrical segments)

        if(name=="collimator"):
            self.name = "collimator"

        elif(name=="electrode"):
            self.name = "electrode"

class collimator(segment):      #collimator subclass
    def __init__(self, l, R):   #length, radius
        super().__init__(l, R, "collimator")

class electrode(segment):       #electroce subclass
    def __init__(self, l, R):
        super().__init__(l, R, "electrode")
###############################################################################################################
###############################################################################################################




###############################################################################################################
###############################################################################################################
#accelerator object class
class accelerator:
    def __init__(self):
        self.particles = []     #all created particles
        self.beam = []          #particles which still remain in the beam
        self.segments = []      #list of segments
        self.seg_positions1 = []        #segment start positions
        self.seg_positions2 = []        #segment stop positions
        self.d = []     #segment distances
        self.R = 0  #electrode radii
        self.t0 =  0    #average time to get to the end of the first collimator

        self.v1 = []
        self.E1 = []
        
        self.v2 = []
        self.E2 = []

    def reset(self, collimator="keep"):     #reset the accelerator (deletes all segments and particles)
        self.particles = []
        self.beam = []
        if collimator=="keep":
            self.segments = [self.segments[0]]
            self.seg_positions1 = [self.seg_positions1[0]]
            self.seg_positions2 = [self.seg_positions2[0]]
        else:
            self.segments = []
            self.seg_positions1 = []
            self.seg_positions2 = []
        self.d = []

    def remove_particles(self):     #remove the particles but keep the segments
        self.particles = []
        self.beam = []

    def add_particle(self, particle):      #add a new particle
        self.particles.append(particle)
        self.beam.append(particle)

    def add_segment(self, seg):         #add a new segment
        self.segments.append(seg)
        if(len(self.seg_positions1)==0):
            self.seg_positions1.append(0)
            self.seg_positions2.append(seg.l)
        else:
            self.seg_positions1.append(self.seg_positions2[-1]+self.d[-1])
            self.seg_positions2.append(self.seg_positions1[-1]+seg.l)

        if(seg.name=="collimator"):
            seg.efac = 0

        elif(seg.name=="electrode"):
            seg.efac = (-1)**(len(self.segments))*np.sign(self.particles[0].q)

    def pop_segment(self):      #delete the last segment
        self.segments.pop()
        self.seg_positions1.pop()
        self.seg_positions2.pop()
        self.d.pop()

    def change_to_collimator(self, l, R):   #changes the last segment to a collimator
        self.seg_positions2[-1] = self.seg_positions2[-1]-self.segments[-1].l+l
        self.segments[-1].l = l
        self.segments[-1].R = R


    def set_electrode_distance(self, d):        #add a new distance between the electrodes
        self.d.append(d)

    def set_electrode_radius(self, R):      #set the segment radius
        self.R = R

    # def particle_position(self, current_particle):
    #     for i, x in enumerate(self.seg_positions2):
    #         if current_particle.r[-1, 2]>x:
    #             current_particle.segment = i+1


###############################################################################################################
#evolve functions
    #find the median velocity after the first collimator
    def vtest_evolve(self): 
        steps = 0
        while self.beam[0].r[-1, 2]<1.5*self.segments[0].l:
            N = len(self.beam)
            j = 0
            while j<N:
                if steps==0:
                    self.beam[j].starts = 1
                    self.beam[j].dtc *= 10
                #moving the particle
                self.beam[j].move()

                #checking if the particle hit the boundary
                if self.beam[j].stops<len(self.segments):
                    if self.beam[j].r[-1,0]**2+self.beam[j].r[-1,1]**2>self.segments[self.beam[j].stops].R**2:
                        self.beam.pop(j)
                        N -= 1
                        continue
                #self.particle_position(self.beam[j])
                j += 1
            steps += 1
            if(len(self.beam)==0):
                break
        # v_part = np.array([])
        # for particle in self.beam:
        return np.median([self.beam[i].v[-1,2] for i in range(len(self.beam))])


    #construct the accelerator from the median particle
    # def testrun_evolve(self, Nel):
    #     steps = 0
        
    #     while 1:
    #         N = len(self.beam)
    #         j = 0
    #         while j<N:
    #             #checking where the particle is
    #             if steps==0:
    #                 self.beam[j].starts = 1
    #                 self.beam[j].d = self.d[0]
    #             else:
    #                 if self.beam[j].starts != len(self.seg_positions1):
    #                     if self.beam[j].r[-1,2]>self.seg_positions1[self.beam[j].starts] and self.beam[j].r[-2,2]<self.seg_positions1[self.beam[j].starts]:
    #                         self.beam[j].starts += 1
    #                 if self.beam[j].stops != len(self.seg_positions2):     
    #                     if self.beam[j].r[-1,2]>self.seg_positions2[self.beam[j].stops] and self.beam[j].r[-2,2]<self.seg_positions2[self.beam[j].stops]:
    #                         #self.x = self.seg_positions1[self.beam[j].stops]
    #                         self.beam[j].d = self.d[self.beam[j].stops]
    #                         self.beam[j].stops += 1
    #                         checkpoint = len(self.beam[j].t)

    #             #moving the particle
    #             self.beam[j].move()
    #         steps += 1
            
    #     return steps

    def testrun_evolve(self, Nel):
        eps1 = 0.035      #field error
        eps2 = 0.03      #separation error

        #1st collimator
        self.beam[0].starts = 1
        self.beam[0].stops = 0

        if len(self.beam)!=1:
            print("Error, the number of particles should be 1! it is currently: ", len(self.beam))

        print(self.beam[0].r[-1, 2], self.beam[0].v[-1, 2])
        print(self.seg_positions2[-1])

        while self.beam[0].r[-1, 2]<self.seg_positions2[-1]:
            #if len(self.beam[0].t)%10000==0:
                #print(self.beam[0].r[2,-1], self.beam[0].v[2,-1])
            self.beam[0].move()
            if self.beam[0].r[-1, 2]<self.beam[0].r[-2, 2]:
                print("ERROR")

        
        self.beam[0].dtc *= 0.01   #lowering the dt in case E is very strong
        dtcstart = self.beam[0].dtc     #saving it for later (variable dt)
        
        print(self.beam[0].r[-1, 2])
        self.t0 = self.beam[0].t0 = self.beam[0].t[-1]
        self.beam[0].stops = 1

        #1st optimal electrode distance guess
        v = 0.1*(9.7*self.beam[0].v[-1,2]+0.3*c)
        d = v*np.pi/omega

        #current number of electrodes
        Nc = 0


        #need an E for E_prev in the loop
        self.beam[0].d = d
        E = self.beam[0].E_field()[2]

        
        #velocities for comparing to analytical solutions
        v_antests = [self.beam[0].v[-1,2]]


        while Nc<Nel:
            #adding a test segment
            self.beam[0].d = d
            self.d.append(d)
            self.add_segment(electrode(1, self.R))
            

            steps = 0
            Emax = self.beam[0].Emax()

            # if steps%1000==0:
            #     print("z%1000 =", self.beam[0].r[-1,2])
            #print(self.seg_positions1[-1])

            print("d=", d)
            print("E/Emax = ", E/Emax)

            self.beam[0].En = 0
            
            self.beam[0].dtc = dtcstart     #reset dtc

            while 1:
                #print("z_before=",self.beam[0].r[-1,2])
                #self.beam[0].move(met="energy")
                self.beam[0].move()
                #self.beam[0].move(met="euler")
                #print("z_after=",self.beam[0].r[-1,2])
                steps += 1

                E_prev = E
                E = self.beam[0].E_field()[2]


                if self.beam[0].v[-1,2] > 0.9*c:
                    self.beam[0].dtc = 100*dtcstart     #increase the step

                #distance to the end of the electrode as a proportion of the electrode length
                d_end = np.abs((self.beam[0].r[-1,2]-self.seg_positions1[-1])/d)

                if self.beam[0].r[-1, 2] > self.seg_positions1[-1]:
                    d = 1.1*d
                    print("z=",self.beam[0].r[-1,2]-self.seg_positions2[-2])
                    self.pop_segment()
                    print("v=", self.beam[0].v[-1,2])
                    print("dt=", self.beam[0].t[-1]-self.beam[0].t[-2])
                    print("steps=", steps)
                    print("electrode: ", len(self.segments))
                    self.beam[0].remove_steps(steps)

                    #for E_prev in the next loop
                    E = self.beam[0].E_field()[2]

                    print("TOO SHORT!")
                    break

                if np.abs(E/Emax) < eps1 and d_end < eps2:
                    #setting the correct length of the electrode with E=0 inside

                    while self.beam[0].r[-1, 2] < self.seg_positions1[-1]:
                        self.beam[0].move()

                    #OLD SEGMENTS
                    # print("OLD SEGMENTS:")
                    # print("starts: ", self.seg_positions1)
                    # print("stops: ", self.seg_positions2)

                    self.pop_segment()

                    #AFTER POP
                    #print("AFTER POP:")
                    #print("starts: ", self.seg_positions1)
                    #print("stops: ", self.seg_positions2)
                    #d = (1-d_end)*d

                    self.d.append(d)
                    self.add_segment(electrode(self.beam[0].v[-1, 2]*np.pi/omega, self.R))

                    # NEW SEGMENTS
                    # print("\nNEW SEGMENTS:")
                    # print("starts: ", self.seg_positions1)
                    # print("stops: ", self.seg_positions2)


                    Nc += 1
                    # print("dt=", self.beam[0].t[-1]-self.beam[0].t[-2])
                    # print("v=", self.beam[0].v[-1,2])
                    print("d_gaps", d+self.segments[-1].l)
                    print("beta=", self.beam[0].v[-1,2]/c)

                    #preparing for the next electrode
                    v = 0.1*(9.7*self.beam[0].v[-1,2]+0.3*c)
                    d = v*np.pi/omega

                    self.beam[0].starts += 1

                    v_antests.append(self.beam[0].v[-1,2])

                    while self.beam[0].r[-1, 2] < self.seg_positions2[-1]:
                        self.beam[0].move()

                    self.beam[0].stops += 1
                    
                    #print("beta_new_el", self.beam[0].v[-1,2]/c)
                    #for E_prev in the next loop
                    E = self.beam[0].E_field()[2]

                    print("steps=", steps)
                    print("SUCCESS!\n\n\n")

                    self.beam[0].E_MeV.append(Ekin(self.beam[0].v[-1,2], self.beam[0].m))

                    break

                if E*self.beam[0].q < 0 and E_prev*self.beam[0].q > 0:
                    d = 0.9*d
                    print("z=",self.beam[0].r[-1,2]-self.seg_positions2[-2])
                    self.pop_segment()
                    print("v=", self.beam[0].v[-1,2])
                    #print("E=",E)
                    #print("Emax=", self.beam[0].Emax())
                    #print("q_sign=", np.sign(self.beam[0].q))
                    # print("t0=", self.beam[0].t0)
                    # print("t=", self.beam[0].t[-1])
                    print("dt=", self.beam[0].t[-1]-self.beam[0].t[-2])
                    print("steps=", steps)
                    print("electrode: ", len(self.segments))
                    print("TOO LONG!")
                    self.beam[0].remove_steps(steps)

                    #for E_prev in the next loop
                    E = self.beam[0].E_field()[2]

                    break

                

        print("velocities", v_antests)
        print("d", self.d)
        print("seg.l", [self.segments[i].l for i in range(len(self.segments))])
        for i in range(len(self.d)):
            print("vtest ratio: ", vratio(self.beam[0].q,self.beam[0].m,self.d[i], v_antests[i], v_antests[i+1]))



    def evolve(self):
        steps = 0
        mod_dt = 0      #dt modification
        N_gaps = len(self.d)
        #print("EVOLVE1")
        #print(self.beam[0].r[-1, 2], self.seg_positions1[-1], self.seg_positions2[-1])
        while self.beam[0].r[-1, 2]<self.seg_positions2[-1]+(self.seg_positions2[-1]-self.seg_positions1[-1])*4:      #steps loop
            N = len(self.beam)

            if mod_dt==0 and self.beam[0].t[-1]>0.95*self.t0:
                mod_dt = 1
                for j in range(N):
                    self.beam[j].dtc *= 0.01
                
            elif mod_dt==1 and self.beam[0].v[-1,2]>0.92*c:
                mod_dt = 2
                for j in range(N):
                    self.beam[j].dtc *= 100

            j = 0

            #print("EVOLVE2")
            while j<N:      #particles loop
                #checking where the particle is
                if steps==0:
                    self.beam[j].starts = 1
                    self.beam[j].d = self.d[0]
                    self.beam[j].t0 = self.t0
                else:
                    if self.beam[j].starts != len(self.seg_positions1):
                        if self.beam[j].r[-1,2]>self.seg_positions1[self.beam[j].starts] and self.beam[j].r[-2,2]<self.seg_positions1[self.beam[j].starts]:
                            self.beam[j].starts += 1

                            # if self.beam[j].starts==len(self.segments):
                            #     self.v2.append(np.sqrt(np.dot(self.beam[j].v[-1], self.beam[j].v[-1])))
                            #     self.E2.append(Ekin(self.v2[-1], self.beam[j].m))
                            #     self.v2[-1] *= 1/c
                            #print("A")

                    if self.beam[j].stops != len(self.seg_positions2):     
                        if self.beam[j].r[-1,2]>self.seg_positions2[self.beam[j].stops] and self.beam[j].r[-2,2]<self.seg_positions2[self.beam[j].stops]:
                            #self.x = self.seg_positions1[self.beam[j].stops]
                            if self.beam[j].stops<N_gaps:
                                self.beam[j].d = self.d[self.beam[j].stops]
                            self.beam[j].stops += 1

                            if self.beam[j].stops==1:
                                self.v1.append(np.sqrt(np.dot(self.beam[j].v[-1], self.beam[j].v[-1])))
                                self.E1.append(Ekin(self.v1[-1], self.beam[j].m))
                                self.v1[-1] *= 1/c

                            if self.beam[j].stops==len(self.segments):
                                self.v2.append(np.sqrt(np.dot(self.beam[j].v[-1], self.beam[j].v[-1])))
                                self.E2.append(Ekin(self.v2[-1], self.beam[j].m))
                                self.v2[-1] *= 1/c
                            #print("B")

            
                #moving the particle
                self.beam[j].move()
                #print("EVOLVE3")

                #checking if the particle hit the boundary
                if self.beam[j].stops<len(self.segments):
                    if self.beam[j].r[-1,0]**2+self.beam[j].r[-1,1]**2>self.segments[self.beam[j].stops].R**2:

                        #to make the plots more realistic
                        self.beam[j].r[-1] = (self.beam[j].r[-1]+3*self.beam[j].r[-2])/4
                        self.beam[j].r[-1] *= (self.beam[j].r[-1, 2]**2+self.segments[self.beam[j].stops].R**2)/np.dot(self.beam[j].r[-1], self.beam[j].r[-1])
                        # self.beam[j].r = np.delete(self.beam[j].r, -1, 0)
                        #print(self.beam[j].r[-1])
                        self.beam.pop(j)
                        N -= 1
                        continue
                #self.particle_position(self.beam[j])
                j += 1
            steps += 1
            if(len(self.beam)==0):
                break
        return steps
