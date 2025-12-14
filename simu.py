import numpy as np
from scipy.optimize import *
import matplotlib.pyplot as plt
import math

m = 1.75
g = 9.81
Ixx = 0.0007309
Iyy = 0.0006644
Izz = 0.0012558
Ir = 3.357*10**(-5)
dx = 0.125 
dy = 0.125#機体中心からプロペラ中心までの距離
DTR = 1/57.3; RTD = 57.3

# Simulation time and model parameters
tstep = 0.03            # Sampling time (sec)
simulation_time = 12.78# Length of time to run simulation (sec)
t = np.arange(0,simulation_time,tstep)   # time array
max_angle_x = math.pi*35/180
max_angle_y = math.pi*5/180
max_angle_z = math.pi*60/180


# Model size
n_states = 12  # Number of states
n_inputs = 4   # Number of inputs

R = 0.12
A = np.pi * R **2

kt = 1.12e-7# 0.1*1.225*(2*0.12)**4/3600
bt = 1.35e-8#0.05*1.225*(2*0.12)**5/3600
cq = bt

# Initialize State Conditions
x = np.zeros((n_states,np.size(t)))  # time history of state vectors
tau = np.zeros((3,np.size(t)))  # time history of state vectors
tu = np.zeros((4,np.size(t)))
th = np.zeros((4,np.size(t)))
pos = np.zeros((3,np.size(t)))
speeds = np.zeros((3,np.size(t)))
# Initial height
#x[6,0] = -max_angle_x
#x[7,0] = max_angle_y
#x[8,0] = 0

x[9,0] = 1.0
x[10,0] =1.0
x[11,0] = 1.5

pos[0,0] = 1.0
pos[1,0] = 1.3
pos[2,0] = 1.5




# Initialize inputs
u = np.zeros((n_inputs,np.size(t)))  # time history of input vectors
# Initial control inputs
u[:,0] = np.zeros(4)

def ThrustEqn(vi, *prop_params):
    #unpack parameters
    R,A,rho,a,b,c,eta, theta0, theta1,U,V,W, Omega = prop_params
    
    Vprime = np.sqrt(U**2 + V**2 + (W-vi)**2)
    
    Thrust = 1/4 * rho * a * b * c * R * \
        ( (W - vi) * Omega * R + 2/3 * (Omega * R)**2 * (theta0 + 3/4 * theta1) + \
          (U**2 + V**2) * (theta0 + 1/2 * theta1) )

    residual = eta * 2 * vi * rho * A * Vprime -Thrust
    
    return residual


def Fthrust(x, u, dx, dy):
    #Propeller configuration parameters
    R = 0.12
    A = np.pi * R **2
    rho = 1.225
    a = 5.7
    b = 2
    c = 0.0243
    eta = 1
    
    #Manufacturer propeller length x pitch
    p_diameter = 9.4488 #inches
    p_pitch = 5 #inches
    
    theta0 = 2 * np.arctan2(p_pitch, (2 * np.pi * 3/4 * p_diameter/2))
    theta1 = -4/3 * np.arctan2(p_pitch, 2 * np.pi * 3/4 * p_diameter/2)
    
    ub,vb,wb = x[0],x[1],x[2]
    p, q, r = x[3], x[4], x[5]
    
    #Transform velocity to local propeller location
    
    U = ub -r * dy
    V = vb + r * dx
    W = wb-q * dx + p * dy

    #Convert commaneded RPM to rad/s
    Omega = 2 * np.pi / 60 * u
    
    #Collecy propeller config, state, and input parameters
    prop_params = (R, A, rho, a, b, c, eta, theta0, theta1, U, V, W, Omega)
    
    #Numrically solve for propeller induced velocity, vi
    #using nonlinear root finder, fsolve and prop_params
    vi0 = 0.01 # initial guess for vi
    vi = fsolve(ThrustEqn, vi0, args=prop_params,maxfev=1000)
    #vi = broyden1(ThrustEqn,vi0,prop_params)
    #vi = ThrustEqn(vi0, prop_params)
    Vprime = np.sqrt(U**2 + V**2 + (W-vi)**2)
    Thrust = eta * 2 * vi * rho * A * Vprime
    
    return Thrust


# Torque function
def T(F,dx,dy):
    # Returns torque about cg given thrust force and dx,dy distance from cg
    
    #### PLACEHOLDER ####
    return F*cq

def stateDerivative(x, u):
    
    #Store variables is  a readble format
    ub = x[0]
    vb = x[1]
    wb = x[2]
    p = x[3]
    q = x[4]
    r = x[5]
    phi = x[6]
    theta = x[7]
    psi = x[8]
    xE = x[9]
    yE = x[10]
    hE = x[11]
    
    #calculate forces propeller inputs
    #F1 = Fthrust(x, u[0], dx, dy)
    #F2 = Fthrust(x, u[1],dx,-dy)
    #F3 = Fthrust(x, u[2],-dx,-dy)
    #F4 = Fthrust(x, u[3], -dx, dy)
    F1 = kt*u[0]**2
    F2 = kt*u[1]**2
    F3 = kt*u[2]**2
    F4 = kt*u[3]**2
    #print([F1,F2,F3,F4])
    Fz = F1 + F2 + F3 + F4
    L = (F2 + F3) * dy - (F1 + F4) *dy#tau phi
    M = (F1 + F2) * dx - (F3 + F4) *dx#tau theta 
    #L = (F1) * dy - (F4) *dy#tau phi
    #M = (F3) * dx - (F1) * dx#tau theta 
    #Tってなんの関数?推進力->プロペラの推力とその半径によって回転方向にトルクを与えるものを関数Tとして->ヨーモーメントを表見してるらしい
    #N = -T(F1, dx, dy) + T(F2, dx, dy) - T(F3, dx, dy) + T(F4, dx, dy) #tau psi
    omega = -u[0] + u[1]-u[2]+u[3]
    N = cq*(-F1 + F2 - F3 + F4)/kt #tau psi
    #print("F1 "+ str(F1)+ " F2 "+str(F2)+" F3 "+str(F3)+" F4 " +str(F4))
    #print("L "+ str(L)+ " M "+str(M)+" N "+str(N))
    global tau
    tau[0] = L
    tau[1] = M
    tau[2] = N

    #L *=100
    #M *=100
    #N *=100
    cphi = np.cos(phi)
    sphi = np.sin(phi)
    cthe = np.cos(theta)
    sthe = np.sin(theta)
    cpsi = np.cos(psi)
    spsi = np.sin(psi)
    
    x_dot = np.zeros(12)
    
    #x_dot[0] = (Fz/m) *(cphi * sthe * cpsi + sphi * spsi)-0.25*ub/m
    #x_dot[1] = (Fz/m) *(cphi * sthe * spsi - sphi * cpsi)- 0.25*vb/m
    #x_dot[2] =   -g+ (Fz/m) * (cphi * cthe)-0.25*wb/m 
    x_dot[0] = -g * sthe + r * vb - q * wb -0.25*ub/m#u_dot
    x_dot[1] = g * sphi * cthe - r * ub + p * wb - 0.2*vb/m# v_dot
    x_dot[2] = 1/m * (-Fz) + g * cphi *cthe + q *ub - p * vb-0.25*wb/m # w_dot
    
    x_dot[3] = 1/Ixx *(L + (Iyy - Izz)  * q * r)#- q*Ir/Ixx #p_dot
    x_dot[4] = 1/Iyy *(M + (Izz - Ixx) * p * r)# + p*Ir/Iyy #q_dot
    x_dot[5] = 1/Izz *(N + (Ixx- Iyy) * p * q) #r_dot
    
    x_dot[6] = p + (q * sphi + r*cphi) * sthe / cthe #phi_dot
    x_dot[7] = q * cphi - r *sphi #theta_dot
    x_dot[8] = (q * sphi + r * cphi) / cthe #psi_dot
    
    x_dot[9] = cthe*cpsi * ub + (-cphi*spsi + sphi * sthe*cpsi)*vb + \
                    (sphi*spsi + cphi*sthe*cpsi) * wb#xE_dot
    x_dot[10] = cthe * spsi * ub +(cphi *cpsi + sphi * sthe *spsi) * vb + \
                    (-sphi*cpsi + cphi*sthe*spsi) *wb #yE_dot
    x_dot[11] =  -1 * (-sthe * ub + sphi*cthe * vb + cphi*cthe *wb) #hE_dot

    #x_dot[9] = ub
    #x_dot[10] = vb
    #x_dot[11] = wb
    return x_dot

# 4th Order Runge Kutta Calculation
def RK4(x,u,dt):
    # Inputs: x[k], u[k], dt (time step, seconds)
    # Returns: x[k+1]
    
    # Calculate slope estimates
    K1 = stateDerivative(x, u)
    K2 = stateDerivative(x + K1 * dt / 2, u)
    K3 = stateDerivative(x + K2 * dt / 2, u)
    K4 = stateDerivative(x + K3 * dt, u)
    
    # Calculate x[k+1] estimate using combination of slope estimates
    x_next = x + 1/6 * (K1 + 2*K2 + 2*K3 + K4) * dt
    
    return x_next

def Huristic(a,b,t):
    t = t*0.02
    f = 0
    if 0<=t and t <= b:
        f = a*np.sin(np.pi*t/b)
    elif b <=t and t <=3*b:
        f = -a*np.sin(np.pi*t/b - np.pi)
    elif 3*b <=t and t <=4*b:
        f = a*np.sin(np.pi*t/b - 3*np.pi)
    #f = t*1.2
    return f

class PID:
    def __init__(self, kp, ki, kd, dt):
        self.kp = kp
        self.ki = ki
        self.kd = kd
        self.dt = dt
        self.integral = 0
        self.p_error = 0
        self.p_output = 0
        self.p_deriv = 0
        self.low_pass_deriv =0
        self.p_y = 0

    def update(self, ref, y, dt):
        error = ref - y
        self.integral +=error*dt
        derrivative =(y - self.p_y)/dt
        self.low_pass_deriv += (derrivative - self.low_pass_deriv)/8
        output = self.kp * error + self.ki * self.integral - self.kd *self.low_pass_deriv
        self.p_error = error
        self.p_y = y
        self.p_output = output
        """
        prop = error - self.p_error
        deriv = prop - self.p_prop
        du =  self.kp * prop + self.ki * error * dt - self.kd * deriv
        self.p_error = error
        self.p_prop = prop
        self.p_output += du*dt
        """
        return self.p_output
    
    def pd(self, error, dt):
        deriv = (error - self.p_error)/dt
        prop = (deriv - self.p_deriv)/dt
        du =  self.kp * error +self.ki*prop + self.kd * deriv 
        self.p_error = error
        self.p_deriv = deriv
        self.p_output = du
        return self.p_output
class Controller:
    def __init__(self):
        Kp_pos = [0.7, 1., 0.015] # proportional [x,y,z]
        Ki_pos = [0.0, 0.00, 0.0013]  # integral [x,y,z]
        Kd_pos = [3.2, 4., 0.2] # derivative [x,y,z]

        # Gains for 
        # angle controller
        Kp_ang= [5, 4, 6] # proportional [x,y,z]
        Ki_ang = [0.00, 0.000, 0.01]  # integral [x,y,z]
        Kd_ang = [3, 3.5, 5] # derivative [x,y,z]
        self.limit = 0.05
        self.flag = 0
        self.ux_t =0
        self.uy_t =0
        self.uz_t =0
        self.roll = 0
        self.pitch = 0
        self.yaw = 0
        self.position = np.array([0, 0, 0.45])
        self.attitude = np.array([0.0, 0.0, 0])
        self.outer_pid_x = PID(Kp_pos[0], Ki_pos[0], Kd_pos[0], 0.02)
        self.outer_pid_y = PID(Kp_pos[1], Ki_pos[1], Kd_pos[1], 0.02)
        self.inner_pid_z = PID(Kp_pos[2], Ki_pos[2], Kd_pos[2], 0.02)
        self.inner_pid_phi = PID(Kp_ang[0], Ki_ang[0], Kd_ang[0], 0.02)
        self.inner_pid_psi = PID(Kp_ang[1], Ki_ang[1], Kd_ang[1], 0.02)
        self.inner_pid_theta = PID(Kp_ang[2], Ki_ang[2], Kd_ang[2],0.02)
    def controller(self,u,x,k,dt):
            global speeds
            ub = x[0,k]
            vb = x[1,k]
            wb = x[2,k]
            phi = x[6,k]
            theta = x[7,k]
            psi = x[8,k]
            cphi = np.cos(phi)
            sphi = np.sin(phi)
            cthe = np.cos(theta)
            sthe = np.sin(theta)
            cpsi = np.cos(psi)
            spsi = np.sin(psi)
            speeds[0,k] = cthe*cpsi * ub + (-cphi*spsi + sphi * sthe*cpsi)*vb + \
                            (sphi*spsi + cphi*sthe*cpsi) * wb#xE_dot
            speeds[1,k] = cthe * spsi * ub +(cphi *cpsi + sphi * sthe *spsi) * vb + \
                            (-sphi*cpsi + cphi*sthe*spsi) *wb #yE_dot
            speeds[2,k] =  -1 * (-sthe * ub + sphi*cthe * vb + cphi*cthe *wb) #hE_dot
            #self.position[0] = np.clip(-x[9,0]*(k+1)/300 +x[9,0],0,np.inf)
            #self.position[1] = np.clip(-x[10,0]*(k+1)/300 +x[10,0],0,np.inf)
            #print(self.position)
            #print("k ="+ str(k))
            #pos[0,k] = self.position[0]
            #pos[1,k] = self.position[1]
            if x[11, k] < 0.45 and np.abs(x[9,k]) <self.limit and np.abs(x[10,k]) <self.limit and self.flag == 0:
                #self.position[1] = x[9,k]
                #self.position[1] = x[10,k]
                #t_ux = x[0,k]
                #t_uy = x[1,k]
                #self.ux_t = self.position[0] - x[9,k]
                #self.ux_t = self.position[1] - x[10,k]
                self.flag = 1
                self.position[2] = 0
                #self.uz_t = (self.position[2] - x[11,k])/10

            error_x = self.position[0] - x[9,k]
            error_y = self.position[1] - x[10,k]
            error_z = self.position[2] - x[11,k]
            #ud = error_x*np.cos(x[8,k]) + error_y*np.sin(x[8,k])
            #vd = error_y*np.cos(x[8,k]) - error_x*np.sin(x[8,k])
            #print("speeds")
            #print(x[:3,k])
            #print([ud,vd,error_z])
                
            #error_ud = ud - x[0,k]
            #error_vd = vd - x[1,k]
            if self.flag == 2:
                ux = self.ux_t
                uy = self.uy_t
                uz =self.uz_t
            else:
                ux = self.outer_pid_x.update(self.position[0],x[9,k], dt)
                uy = self.outer_pid_y.update(self.position[1],x[10,k],dt)
                uz = self.inner_pid_z.update(self.position[2],x[11,k],dt)
            if self.flag ==2:
                self.ux_t = ux
                self.uy_t = uy
                #self.uz_t = uz
                self.flag +=1
            if self.flag == 1:
                uz -=0.2
            #dpsi = np.arccos(np.sqrt((self.position[0]-x[9,k])**2+(self.position[1]-x[10,k])**2)/np.sqrt(ux**2+uy**2))
            dpsi = self.attitude[2]
            if x[11, k] >= 0.45 and np.abs(x[9,k]) >= self.limit and np.abs(x[10,k]) >= self.limit:
                dpsi = np.arccos((np.abs(self.position[1] - x[10,k]))/np.sqrt(((self.position[0]-x[9,k])**2+ (self.position[1]-x[10,k])**2 + (self.position[2]-x[11,k])**2)))
            #dpsi = Huristic(0.08,0.01,k)
            #dpsi = np.arccos(np.sqrt((self.position[0]-x[9,k])**2+(self.position[1]-x[10,k])**2)/np.sqrt(ux**2+uy**2))
            #dpsi = np.sin(k)+np.cos(k)
            #dpsi = np.arccos(ux/np.sqrt(ux**2+uy**2))
            #dphi = np.arcsin((ux*np.sin(x[8,k])-uy*np.cos(x[8,k]))/(ux**2+uy**2+(uz+g)**2))
            #dtheta = np.arctan((ux*np.cos(x[8,k])+uy*np.sin(x[8,k]))/(uz+g))
            dphi = np.arcsin((ux*np.sin(dpsi)-uy*np.cos(dpsi))/(ux**2+uy**2+(uz+g)**2))
            dtheta = np.arctan((ux*np.cos(dpsi)+uy*np.sin(dpsi))/(uz+g))
            
            #dpsi = np.arccos(np.sqrt((self.position[0]-x[9,k])**2+(self.position[1]-x[10,k])**2)/np.sqrt(ux**2+uy**2))
            #self.attitude[0] = ud*np.sin(x[8,k])-vd*np.cos(x[8,k])
            #self.attitude[1] = -1*(ud*np.cos(x[8,k])+vd*np.sin(x[8,k]))
            #self.attitude[0] = ux
            #self.attitude[1] = -uy
            self.attitude[0] = -dphi
            self.attitude[1] = -dtheta
            self.attitude[2] = dpsi

            
            
            mag_angle_des = np.linalg.norm(self.attitude)
            if mag_angle_des > max_angle_x:
                self.attitude[0] = (self.attitude[0] / mag_angle_des) * max_angle_x
            if mag_angle_des > max_angle_y:
                self.attitude[1] = (self.attitude[1] / mag_angle_des) * max_angle_y
            
            if mag_angle_des > max_angle_z:
                self.attitude[2] = (self.attitude[2] / mag_angle_des) * max_angle_z
            
            #else:
            #    self.attitude[0] = 0
            #    self.attitude[1] = 0
            #print(self.attitude)
            error_phi = self.attitude[0] - x[6,k]
            error_theta = self.attitude[1] - x[7,k]
            error_psi = self.attitude[2] - x[8,k]
            #T = m*(ux*(np.sin(x[7,k])*np.cos(x[8,k])*np.cos(x[6,k])+np.sin(x[6,k])*np.sin(x[8,k]))+uy*(np.sin(x[7,k])*np.sin(x[8,k])*np.cos(x[6,k])-np.cos(x[8,k])*np.sin(x[6,k]))+(uz+g)*np.cos(x[7,k])*np.cos(x[6,k]))
            #thrust = np.clip(T,0.0,13)
            thrust = np.clip((g+uz)*m/(np.cos(self.attitude[0])*np.cos(self.attitude[1])),0.0,17.5)
            if self.flag == 2:
                torque_x = self.roll
                torque_y = self.pitch
                torque_z = self.yaw
            else:
                torque_x = self.inner_pid_phi.update(self.attitude[0],x[6,k],dt)*Ixx
                torque_y = self.inner_pid_psi.update(self.attitude[1],x[7,k],dt)*Iyy
                torque_z = self.inner_pid_theta.update(self.attitude[2],x[8,k],dt)*Izz
            
            if self.flag == 1:
                self.roll = torque_x
                self.pitch = torque_y
                self.yaw = torque_z
            
            #torque_z = 0
            #x[6,k] +=0.0001*k
            global tu
            tu[:,k] = [thrust,torque_x, torque_y, torque_z]
            
            
            l = dx
            
            motor_torque_1 = np.clip(0.25*thrust/kt - 0.25*torque_x /(l*kt) + 0.25*torque_y/
                                    (l*kt) - 0.25 *torque_z/bt, 0, np.inf)
            motor_torque_2 = np.clip(0.25*thrust/kt + 0.25*torque_x /(l*kt) + 0.25*torque_y/
                                    (l*kt) + 0.25 *torque_z/bt,0, np.inf)
            motor_torque_3 = np.clip(0.25*thrust/kt + 0.25*torque_x /(l*kt) - 0.25*torque_y/
                                    (l*kt) - 0.25 *torque_z/bt,0, np.inf)
            motor_torque_4 = np.clip(0.25*thrust/kt - 0.25*torque_x /(l*kt) - 0.25*torque_y/
                                    (l*kt) + 0.25 *torque_z/bt,0, np.inf)
            motor_speeds = [motor_torque_1, motor_torque_2, motor_torque_3, motor_torque_4]

            for i in range(4):
                 motor_speeds[i] = np.clip(motor_speeds[i]*kt, 0, 4.5)
                 motor_speeds[i] /=kt
            #motor_speeds = np.array(motor_speeds) * (2*np.pi/60)**2
            
            motor_speed_1 = np.clip(np.power(motor_speeds[0],1/2), 0, np.inf)
            motor_speed_2 = np.clip(np.power(motor_speeds[1],1/2), 0, np.inf)
            motor_speed_3 = np.clip(np.power(motor_speeds[2],1/2), 0, np.inf)
            motor_speed_4 = np.clip(np.power(motor_speeds[3],1/2), 0, np.inf)

            u[0,k] = motor_speed_1
            u[1,k] = motor_speed_2
            u[2,k] = motor_speed_3
            u[3,k] = motor_speed_4
            #u[:,k] = np.array(u[:,k]) * (2*np.pi/60)
            #print(u[:,k])
            #print("=========")
            #print(x[6:9,k])
            global th
            th[:,k] =kt*u[:,k]**2
            #print(th[:,k])
    
            return u

# March through time array and numerically solve for vehicle states
#vertvel = np.array([0,0,1] + 9*[0])
cont = Controller()
max_t  =[1.0,1.3,1.5]
for k in range(0, np.size(t) -1): 
        
    # Determine control inputs based on current state
    #wu[:,k] = controlInputs(x[:,k], t[k])
    
    # Predict state after one time step
    print(x[9:,k])
    #if x[11,k] <=0.45:
    #    print(k*0.03)
    u = cont.controller(u,x,k,tstep)
    #if k == 0:
    #    u[:,k] = 117.3
    #print(u[:,k])
    x[:,k+1] = RK4(x[:,k], u[:,k], tstep)
    #x[7,k+1] = np.clip(x[7,k+1]*RTD,-3,3)/RTD
    
    for i in range(3):
        if max_t[i] > x[9+i,k+1]:
            max_t[i] = x[9+i,k+1]
    
    #pos[:2,k+1] = pos[:2,k]+x[:2,k+1]*0.02
    #pos[2,k+1] = pos[2,k]-x[2,k+1]*0.02

    #print(x[9:,k])
    
    if  x[11,k+1] <= 0 :
        x[11,k+1] =0
        break
    
    #print(tau[:,k])
    #if x[11,k+1] < 0.0:
    #    break
#print(max_t)
    
    

plt.figure(1, figsize=(8,8))
plt.subplot(311)
plt.plot(t,x[9,:],'r',label='X')
plt.plot(t,x[10,:],'b',label='Y')
plt.plot(t,x[11,:],'g',label='Z')

#plt.ylim(-0.05, 0.1)
#plt.xlim(0, 3)
plt.legend(loc='best')
plt.ylabel('Position (m)')
plt.xlabel('Time (s)')
#plt.legend(loc='best')
#plt.text(0.5, -0.38, '(a)', ha='center', va='center', transform=plt.gca().transAxes)  # 下にタイトル

plt.subplot(312)
#plt.plot(t,tu[1,:],'r',label='torque_x')
#plt.plot(t,tu[2,:],'b',label='torque_y')
plt.plot(t,x[8,:]*RTD,'g',label='Psi')
#plt.plot(t,x[9,:],'r',label='x')
#plt.xlim(0, 1)
plt.legend(loc='best')
#plt.ylabel('tau (deg)')
plt.ylabel('Euler Angle (degrees)')
plt.xlabel('Time (s)')

plt.subplot(313)
plt.plot(t,x[6,:]*RTD,'r',label='Phi')
plt.plot(t,x[7,:]*RTD,'b',label='Theta')

#plt.xlim(0, 3)
plt.legend(loc='best')
plt.ylabel('Euler Angles (degrees)')
plt.xlabel('Time (s)')
plt.subplots_adjust(hspace=0.5)
"""
plt.figure(2, figsize=(8,8))
ax = plt.subplot(1,1,1)
plt.plot(x[9,0:-1:20],x[11,0:-1:20],'bo-',label='x')
plt.text(x[9,0] + 0.1, x[11,0],'START')
plt.text(x[9,-1], x[11,-1],'END')
plt.ylabel('h [m]'); plt.xlabel('x [m]')
ax.axis('equal')
#plt.legend(loc='best')
plt.title('Vertical Profile')
"""
plt.figure(3, figsize=(8,4))
plt.plot(t[0:-1],u[0,0:-1],'b',label='T1')
plt.plot(t[0:-1],u[1,0:-1],'g',label='T2')
plt.plot(t[0:-1],u[2,0:-1],'r',label='T3')
plt.plot(t[0:-1],u[3,0:-1],'y',label='T4')
#plt.xlim(0, 1)
plt.xlabel('Time (sec)')
plt.ylabel('Propeller RPM')
plt.legend(loc='best')
plt.title('Time History of Control Inputs')


plt.figure(4, figsize=(8,8))
plt.plot(t,speeds[0,:],'r',label='u')
plt.plot(t,speeds[1,:],'b',label='v')
plt.plot(t,speeds[2,:],'g',label='w')
#plt.xlim(0, 10)
#plt.ylim(-10, 2)
plt.xlabel('Time (sec)')
plt.ylabel('[m/s]')
plt.legend(loc='best')
plt.title('speeds',y=-0.25)

plt.figure(5, figsize=(8,8))
plt.plot(t,tu[0,:],'r',label='thrust')
#plt.xlim(0, 1)
plt.xlabel('Time (sec)')
plt.ylabel('T[N]')
plt.legend(loc='best')
plt.title('Time History of Control Inputs')


fig = plt.figure(6, figsize=(4,5))
axes = fig.add_subplot(1, 1, 1, projection='3d')
axes.plot(x[9,:], x[10,:], x[11,:])
axes.set_title('Flight Path')
axes.set_xlabel('x (m)')
axes.set_ylabel('y (m)')
axes.set_zlabel('z (m)')
#axes.set_xlim(0,1.8)
#axes.set_ylim(-0.1,1.35)
#axes.set_zlim(0,1.55)

plt.figure(7, figsize=(8,4))
plt.plot(t[0:-1],th[0,0:-1],'b',label='T1')
plt.plot(t[0:-1],th[1,0:-1],'g',label='T2')
plt.plot(t[0:-1],th[2,0:-1],'r',label='T3')
plt.plot(t[0:-1],th[3,0:-1],'y',label='T4')
#plt.xlim(0, 1)
plt.xlabel('Time (sec)')
plt.ylabel('Propeller Thrust')
plt.legend(loc='best')
plt.title('Time History of Control Inputs')
"""
plt.figure(8, figsize=(8,8))
plt.plot(t,pos[0,:],'r',label='x')
plt.plot(t,pos[1,:],'b',label='y')
plt.plot(t,pos[2,:],'g',label='z')
#plt.xlim(0, 10)
#plt.ylim(-10, 2)
plt.xlabel('Time (sec)')
plt.ylabel('m')
plt.legend(loc='best')
plt.title('position',y=-0.25)
"""

plt.show()