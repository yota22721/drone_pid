import numpy as np
from scipy.optimize import *
import matplotlib.pyplot as plt
import math

m = 1.28
g = 9.81
Ixx = 0.0007309
Iyy = 0.0006644
Izz = 0.0012558
dx = 0.125 
dy = 0.125#機体中心からプロペラ中心までの距離
DTR = 1/57.3; RTD = 57.3

# Simulation time and model parameters
tstep = 0.02            # Sampling time (sec)
simulation_time = 30   # Length of time to run simulation (sec)
t = np.arange(0,simulation_time,tstep)   # time array
max_angle = math.pi/4


# Model size
n_states = 12  # Number of states
n_inputs = 4   # Number of inputs

kt = 1e-7

# Initialize State Conditions
x = np.zeros((n_states,np.size(t)))  # time history of state vectors
tau = np.zeros((3,np.size(t)))  # time history of state vectors
# Initial height
#x[6,0] = -0.784
#x[7,0] = 0.784
#x[8,0] = 0

x[9,0] = 1.0
x[10,0] =1.3
x[11,0] = 1.5




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
    return F*dx

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
    F1 = Fthrust(x, u[0], dx, dy)
    F2 = Fthrust(x, u[1],-dx,-dy)
    F3 = Fthrust(x, u[2],dx,-dy)
    F4 = Fthrust(x, u[3], -dx, dy)
    #F1 = kt*u[0]**2
    #F2 = kt*u[1]**2
    #F3 = kt*u[2]**2
    #F4 = kt*u[3]**2
    #print([F1,F2,F3,F4])
    Fz = F1 + F2 + F3 + F4
    L = (F2 + F3) * dy - (F1 + F4) *dy#tau phi
    M = (F1 + F3) * dx - (F2 + F4) *dx#tau theta 
    #L = (F1) * dy - (F4) *dy#tau phi
    #M = (F3) * dx - (F1) * dx#tau theta 
    #Tってなんの関数?推進力->プロペラの推力とその半径によって回転方向にトルクを与えるものを関数Tとして->ヨーモーメントを表見してるらしい
    #N = -T(F1, dx, dy) + T(F2, dx, dy) - T(F3, dx, dy) + T(F4, dx, dy) #tau psi
    N = -T(F1, dx, dy) - T(F2, dx, dy) + T(F3, dx, dy) + T(F4, dx, dy) #tau psi
    #print("F1 "+ str(F1)+ " F2 "+str(F2)+" F3 "+str(F3)+" F4" +str(F4))
    #print("L "+ str(L)+ " M "+str(M)+" N "+str(N))
    global tau
    tau[0] = L
    tau[1] = M
    tau[2] = N
    
    cphi = np.cos(phi)
    sphi = np.sin(phi)
    cthe = np.cos(theta)
    sthe = np.sin(theta)
    cpsi = np.cos(psi)
    spsi = np.sin(psi)
    
    x_dot = np.zeros(12)
    
    x_dot[0] = -g * sthe + r * vb - q * wb #u_dot
    x_dot[1] = g * sphi * cthe - r * ub + p * wb# v_dot
    x_dot[2] = 1/m * (-Fz) + g * cphi *cthe + q *ub - p * vb # w_dot
    #x_dot[0] = -g + (Fz/m) *(cphi * sthe * cpsi + sphi * spsi)
    #x_dot[1] = -g + (Fz/m) *(cphi * sthe * spsi - sphi * cpsi)
    #x_dot[2] =  - g+ (Fz/m) * (cphi * cthe)
    
    x_dot[3] = 1/Ixx *(L + (Iyy - Izz)  * q * r) #p_dot
    x_dot[4] = 1/Iyy *(M + (Izz - Ixx) * p * r) #q_dot
    x_dot[5] = 1/Izz *(N + (Ixx- Iyy) * p * q) #r_dot
    
    x_dot[6] = p + (q * sphi + r*cphi) * sthe / cthe #phi_dot
    x_dot[7] = q * cphi - r *sphi #theta_dot
    x_dot[8] = (q * sphi + r * cphi) / cthe #psi_dot
    x_dot[8] = 0
    
    x_dot[9] = cthe*cpsi * ub + (-cphi*spsi + sphi * sthe*cpsi)*vb + \
                    (sphi*spsi + cphi*sthe*cpsi) * wb#xE_dot
    x_dot[10] = cthe * spsi * ub +(cphi *cpsi + sphi * sthe *spsi) * vb + \
                    (-sphi*cpsi + cphi*sthe*spsi) *wb #yE_dot
    x_dot[11] = -1 * (-sthe * ub + sphi*cthe * vb + cphi*cthe *wb) #hE_dot
    
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


class PID:
    def __init__(self, kp, ki, kd, dt):
        self.kp = kp
        self.ki = ki
        self.kd = kd
        self.dt = dt
        self.integral = 0
        self.p_error = 0
        self.p_output = 0

    def update(self, error, dt):
        self.integral +=error*dt
        derrivative =(self.p_output - self.p_error)/dt
        output = self.kp * error + self.ki * self.integral - self.kd * derrivative
        self.p_error = error
        self.p_output = output
        return output
class Controller:
    def __init__(self):
        Kp_pos = [0.31, 0.32, 0.33] # proportional [x,y,z]
        Ki_pos = [0.12, 0.3, 0.11]  # integral [x,y,z]
        Kd_pos = [8, 8, 8] # derivative [x,y,z]

        # Gains for angle controller
        Kp_ang= [0.35, 0.3574, 0.35] # proportional [x,y,z]
        Ki_ang = [0.155, 0.15, 0.15]  # integral [x,y,z]
        Kd_ang = [8, 8, 8.00] # derivative [x,y,z]
        self.position = np.array([0.0, 0.0, 0.0])
        self.attitude = np.array([0.0, 0.0, 0.0])
        self.outer_pid_x = PID(Kp_pos[0], Ki_pos[0], Kd_pos[0], 0.125)
        self.outer_pid_y = PID(Kp_pos[1], Ki_pos[1], Kd_pos[1], 0.125)
        self.inner_pid_z = PID(Kp_pos[2], Ki_pos[2], Kd_pos[2], 0.125)
        self.inner_pid_phi = PID(Kp_ang[0], Ki_ang[0], Kd_ang[0], 0.125)
        self.inner_pid_psi = PID(Kp_ang[1], Ki_ang[1], Kd_ang[1], 0.125)
        self.inner_pid_theta = PID(Kp_ang[2], Ki_ang[2], Kd_ang[2],0.125)
    def controller(self,u,x,k,dt):
        
            error_x = self.position[0] - x[9,k]
            error_y = self.position[1] - x[10,k]
            #if x[11, k] > 0.45:
            t_phi = self.outer_pid_x.update(error_x, dt)
            t_psi = self.outer_pid_y.update(error_y,dt)
            self.attitude[0] = t_phi
            self.attitude[1] = t_psi
            #else:
            #    self.attitude[0] = 0
            #    self.attitude[1] = 0
            
            
            #print(self.attitude)
            error_z = self.position[2] - x[11,k]
            error_phi = self.attitude[0] - x[6,k]
            error_psi = self.attitude[1] - x[8,k]
            error_theta = self.attitude[2] - x[7,k]
            thrust = self.inner_pid_z.update(error_z,dt)
            torque_x = self.inner_pid_phi.update(error_phi,dt)
            torque_y = self.inner_pid_psi.update(error_psi,dt)
            torque_z = self.inner_pid_theta.update(error_theta,dt)
            #torque_z = 0
            #x[6,k] +=0.0001*k
            print([torque_x, torque_y, torque_z])
            #print([x[6,k], x[8,k], x[7,k]])
            
            
            trim = 1702
            l = dx
            
            motor_torque_1 = np.clip(0.25*thrust + 0.5*torque_x /(l) - 0.5*torque_y/
                                    (l) + 0.5 *torque_z, 0, np.inf)
            motor_torque_2 = np.clip(0.25*thrust - 0.5*torque_x /(l) - 0.5*torque_y/
                                    (l) - 0.5 *torque_z,0, np.inf)
            motor_torque_3 = np.clip(0.25*thrust - 0.5*torque_x /(l) + 0.5*torque_y/
                                    (l) + 0.5 *torque_z,0, np.inf)
            motor_torque_4 = np.clip(0.25*thrust + 0.5*torque_x /(l) + 0.5*torque_y/
                                    (l) - 0.5 *torque_z,0, np.inf)
            motor_speeds = [motor_torque_1, motor_torque_2, motor_torque_3, motor_torque_4]
            """
            maxT = 17.5 #  max thrust from any single motor, N
            minT = .5 # min thrust from any single motor, N 
            # Ensure that desired thrust is within overall min and max of all motors
            thrust_all = np.array(motor_speeds) * (kt)
            #print(thrust_all)
            over_max = np.argwhere(thrust_all > maxT)
            under_min = np.argwhere(thrust_all < minT)

            if over_max.size != 0:
                for i in range(over_max.size):
                    motor_speeds[over_max[i][0]] = maxT / (kt)
            if under_min.size != 0:
                for i in range(under_min.size):
                    motor_speeds[under_min[i][0]] = minT / (kt)
            """        
            motor_speed_1 = np.clip(np.power(motor_speeds[0],1/2), 0, 5000)
            motor_speed_2 = np.clip(np.power(motor_speeds[1],1/2), 0, 5000)
            motor_speed_3 = np.clip(np.power(motor_speeds[2],1/2), 0, 5000)
            motor_speed_4 = np.clip(np.power(motor_speeds[3],1/2), 0, 5000)
            
            
            #motor_speed_1 = np.clip(np.sqrt(motor_torque_1), 0, 5000)
            #motor_speed_2 = np.clip(np.sqrt(motor_torque_2), 0, 5000)
            #motor_speed_3 = np.clip(np.sqrt(motor_torque_3), 0, 5000)
            #motor_speed_4 = np.clip(np.sqrt(motor_torque_4), 0, 5000)
            
            """
            e1 = torque_x * Ixx
            e2 = torque_y * Iyy
            e3 = torque_z * Izz

            #less typing
            n = 4

            # Thrust desired converted into motor speeds
            weight_speed = thrust / (n*kt)
            #print([e1, e2, e3, weight_speed])

            # Thrust differene in each motor to achieve needed torque on body
            motor_speeds = []
            motor_speeds.append(weight_speed + (e2/((n/2)*kt*dx)) - (e3/(n*kt)))
            motor_speeds.append(weight_speed - (e1/((n/2)*kt*dx)) - (e3/(n*kt)))
            motor_speeds.append(weight_speed - (e2/((n/2)*kt*dx)) + (e3/(n*kt)))
            motor_speeds.append(weight_speed + (e1/((n/2)*kt*dx)) + (e3/(n*kt)))
            #print(motor_speeds)
            maxT = 12.5 #  max thrust from any single motor, N
            minT = .5 # min thrust from any single motor, N 
            # Ensure that desired thrust is within overall min and max of all motors
            thrust_all = np.array(motor_speeds) * (kt)
            #print(thrust_all)
            over_max = np.argwhere(thrust_all > maxT)
            under_min = np.argwhere(thrust_all < minT)

            if over_max.size != 0:
                for i in range(over_max.size):
                    motor_speeds[over_max[i][0]] = maxT / (kt)
            if under_min.size != 0:
                for i in range(under_min.size):
                    motor_speeds[under_min[i][0]] = minT / (kt)
            #print(motor_speeds)
            
         
            motor_speed_1 = np.clip(np.power(motor_speeds[0],1/2), 0, 5000)
            motor_speed_2 = np.clip(np.power(motor_speeds[1],1/2), 0, 5000)
            motor_speed_3 = np.clip(np.power(motor_speeds[2],1/2), 0, 5000)
            motor_speed_4 = np.clip(np.power(motor_speeds[3],1/2), 0, 5000)
            """
            u[0,k] = motor_speed_1 +trim
            u[1,k] = motor_speed_2 +trim
            u[2,k] = motor_speed_3 +trim
            u[3,k] = motor_speed_4 +trim
            print(u[:,k])
            print("=========")
            #print(x[6:9,k])
        
    
            return u

# March through time array and numerically solve for vehicle states
#vertvel = np.array([0,0,1] + 9*[0])
cont = Controller()
for k in range(0, np.size(t) - 1): 
        
    # Determine control inputs based on current state
    #wu[:,k] = controlInputs(x[:,k], t[k])
    
    # Predict state after one time step
    print(x[9:,k])
    mag_angle_des = np.linalg.norm(x[6:9,k])
            
    if mag_angle_des > max_angle:
        x[6:9,k] = (x[6:9,k] / mag_angle_des) * max_angle
    print("-------------")
    #print(u[:,k])
    #print("L "+ str(tau[0,k])+ " M "+str(tau[1,k])+" N "+str(tau[2,k]))
    u = cont.controller(u,x,k,tstep)
    #print(u[:,k])
    x[:,k+1] = RK4(x[:,k], u[:,k], tstep)
    #print(x[9:,k])
    #if  x[11,k+1] <= 0 :
    #print(tau[:,k])
    if x[11,k+1] < 0:
        break
    

plt.figure(1, figsize=(8,8))
plt.subplot(311)
plt.plot(t,x[9,:],'r',label='x')
plt.plot(t,x[10,:],'b',label='y')
plt.plot(t,x[11,:],'g',label='z')

#plt.ylim(-0.05, 0.1)
plt.xlim(0, 1.5)
plt.legend(loc='best')
plt.ylabel('z (m)')
#plt.xlabel('Time (sec)')
#plt.legend(loc='best')
plt.title('Time History of Height, X Position, and Pitch')

fig = plt.figure(1, figsize=(8,8))
axes = fig.add_subplot(2, 4, 1, projection='3d')
axes.plot(x[9,:], x[10,:], x[11,:])
axes.set_title('Flight Path')
axes.set_xlabel('x (m)')
axes.set_ylabel('y (m)')
axes.set_zlabel('z (m)')

plt.subplot(312)
plt.plot(t,tau[0,:],'r',label='roll')
plt.plot(t,tau[1,:],'b',label='pitch')
plt.plot(t,tau[2,:],'g',label='yaw')
#plt.plot(t,x[9,:],'r',label='x')
plt.xlim(0, 1)
plt.legend(loc='best')
#plt.ylabel('tau (deg)')
plt.ylabel('tau (m)')
#plt.xlabel('Time (sec)')

plt.subplot(313)
plt.plot(t,x[6,:]*RTD,'r',label='phi',marker=".")
plt.plot(t,x[7,:]*RTD,'b',label='theta',marker=".")
plt.plot(t,x[8,:]*RTD,'g',label='psi',marker=".")
plt.xlim(0, 1)
plt.legend(loc='best')
plt.ylabel('Theta (deg)')
plt.xlabel('Time (sec)')

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
plt.xlim(0, 1)
plt.xlabel('Time (sec)')
plt.ylabel('Propeller RPM')
plt.legend(loc='best')
plt.title('Time History of Control Inputs')
"""
plt.figure(4, figsize=(8,8))
plt.plot(t[0:-1],tau[0,0:-1],'r',label='roll')
plt.plot(t[0:-1],tau[1,0:-1],'b',label='pitch')
plt.plot(t[0:-1],tau[2,0:-1],'g',label='yaw')
plt.xlim(0, 1)
plt.xlabel('Time (sec)')
plt.ylabel('tau')
plt.legend(loc='best')
plt.title('Time History of Control Inputs')
"""
"""
i = 0
while True:
    print(tau[:,i])
    if i> 20:
        break
    else:
        i +=1
"""
plt.show()

    