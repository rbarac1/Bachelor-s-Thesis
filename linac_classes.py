import numpy as np

c = 299792458

Vmax = 1e11  #max voltage
omega = 100   #AC frequency

def gamma(v):
    return 1/np.sqrt(1-v*v/c/c)

def dgamma(v):  #dgamma/dv
    return 1/np.sqrt(1-v*v/c/c)**3*v/c/c


#particle object class
class particle:
    def __init__(self, q, m, v0, E=[0.,0.,0.], B=[0.,0.,0.], dt=1e-7):
        self.m = m
        self.q = q
        self.E = np.array(E)
        self.B = np.array(B)

        self.v = np.array([v0])

        self.r = np.array([[0,0,0]])

        self.t = [0]
        self.dt = dt

        self.segment = 0

        self.starts = 0   #number of passed starts of a linac segment
        self.stops = 0   #number of passed stops of a linac segment
        self.d = 0  #distance of electrodes in the accelerator

        #self.x = 0  #start of the segment with an electric field

    def __ac(self, v):
        if self.starts > self.stops:
            return 0
        else:
            #self.E = np.array([self.q*Vmax*2/self.d*np.cos(omega*self.t[-1]),self.q*Vmax*2/self.d*np.cos(omega*self.t[-1]),0])
            self.E = np.array([0,0,self.q*Vmax*2/self.d*np.cos(omega*self.t[-1])])
            F_clas = self.q*(self.E+np.cross(v, self.B))
            #print("AAAAAAAAA")
            return 1/self.m/gamma(v)*(F_clas-np.dot(F_clas,v)*v/c/c)
        #return (self.q*(self.E+np.cross(v, self.B))/self.m) / (dgamma(v)*v+gamma(v))

    def move(self, met="RK4"):
        if met=="RK4":
            k1v = self.__ac(self.v[-1])*self.dt
            k1r = self.v[-1]*self.dt

            k2v = self.__ac(self.v[-1]+k1v/2)*self.dt
            k2r = (self.v[-1]+k1v/2)*self.dt

            k3v = self.__ac(self.v[-1]+k2v/2)*self.dt
            k3r = (self.v[-1]+k2v/2)*self.dt

            k4v = self.__ac(self.v[-1]+k3v)*self.dt
            k4r = (self.v[-1]+k3v)*self.dt

            self.v = np.append(self.v, [self.v[-1]+ (k1v+2*k2v+2*k3v+k4v)/6], 0)

            self.r = np.append(self.r, [self.r[-1]+ 1/6*(k1r+2*k2r+2*k3r+k4r)], 0)
                
            self.t += [self.t[-1]+self.dt]
        
        else:
            #print(self.v)
            a = self.__ac(self.v[-1])
            self.v = np.append(self.v, [self.v[-1]+a*self.dt],0)
            self.r = np.append(self.r, [self.r[-1]+self.v[-1]*self.dt],0)
            self.t += [self.t[-1]+self.dt]

    def set_Efield(self, E):
        self.E = np.array(E)

    def reset(self):
        self.v = np.array([self.v[0]])
        self.r = np.array([self.r[0]])
        self.t = [0]

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



class electron(particle):
    def __init__(self, v0, E=[0.,0.,0.], B=[0.,0.,0.], dt=1e-7):
        self.m = 9.1093837e-31
        self.q = -1.60217663e-19
        super(electron, self).__init__(self.q, self.m, v0, E, B, dt)

class proton(particle):
    def __init__(self, v0, E=[0.,0.,0.], B=[0.,0.,0.], dt=1e-7):
        self.m = 1.67262192e-27
        self.q = 1.60217663e-19
        super(electron, self).__init__(self.q, self.m, v0, E, B, dt)



#accelerator parts
class segment:
    def __init__(self, l, R, name):
        self.l = l
        self.R = R

        if(name=="collimator"):
            self.name = "collimator"

        elif(name=="electrode"):
            self.name = "electrode"

class collimator(segment):
    def __init__(self, l, w):
        super().__init__(l, w, "collimator")

class electrode(segment):
    def __init__(self, l, w):
        super().__init__(l, w, "electrode")



#accelerator object class
class accelerator:
    def __init__(self):
        self.particles = []
        self.beam = []
        self.segments = []
        self.seg_positions1 = []
        self.seg_positions2 = []
    
    def add_particle(self, particle):
        self.particles.append(particle)
        self.beam.append(particle)

    def add_segment(self, seg):
        self.segments.append(seg)
        if(len(self.seg_positions1)==0):
            self.seg_positions1.append(0)
            self.seg_positions2.append(seg.l)
        else:
            self.seg_positions1.append(self.seg_positions2[-1]+self.d)
            self.seg_positions2.append(self.seg_positions1[-1]+seg.l)

        if(seg.name=="collimator"):
            seg.efac = 0

        elif(seg.name=="electrode"):
            seg.efac = (-1)**(len(self.segments))*np.sign(self.particles[0].q)



    def set_electrode_distance(self, d):
        self.d = d

    def particle_position(self, current_particle):
        for i, x in enumerate(self.seg_positions2):
            if current_particle.r[-1, 2]>x:
                current_particle.segment = i+1

    def evolve(self):
        steps = 0
        while self.beam[0].r[-1, 2]<self.seg_positions2[-1]+(self.seg_positions2[-1]-self.seg_positions1[-1])*0.5:
            N = len(self.beam)
            j = 0
            while j<N:
                #checking where the particle is
                if steps==0:
                    self.beam[j].starts = 1
                    self.beam[j].d = self.d
                else:
                    if self.beam[j].starts != len(self.seg_positions1):
                        if self.beam[j].r[-1,2]>self.seg_positions1[self.beam[j].starts] and self.beam[j].r[-2,2]<self.seg_positions1[self.beam[j].starts]:
                            self.beam[j].starts += 1
                            print("A")

                    if self.beam[j].stops != len(self.seg_positions2):     
                        if self.beam[j].r[-1,2]>self.seg_positions2[self.beam[j].stops] and self.beam[j].r[-2,2]<self.seg_positions2[self.beam[j].stops]:
                            #self.x = self.seg_positions1[self.beam[j].stops]
                            self.beam[j].stops += 1
                            print("B")

                #moving the particle
                self.beam[j].move()

                #checking if the particle hit the boundary
                if self.beam[j].segment<len(self.segments):
                    if self.beam[j].r[-1,0]**2+self.beam[j].r[-1,1]**2>self.segments[self.beam[j].segment].R**2:
                        self.beam.pop(j)
                        N -= 1
                        continue
                self.particle_position(self.beam[j])
                j += 1
            steps += 1
            if(len(self.beam)==0):
                break
        return steps
