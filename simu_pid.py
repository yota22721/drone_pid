import numpy as np
import matplotlib.pyplot as plt
import plotly.express as px
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

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

# Model size
n_states = 12  # Number of states
n_inputs = 4   # Number of inputs


# Initialize State Conditions
x = np.zeros((n_states,np.size(t)))  # time history of state vectors

# Initial height
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
    vi = fsolve(ThrustEqn, vi0, args=prop_params)
    Vprime = np.sqrt(U**2 + V**2 + (W-vi)**2)
    Thrust = eta * 2 * vi * rho * A * Vprime
    
    return Thrust
# Torque function
def T(F,dx,dy):
    # Returns torque about cg given thrust force and dx,dy distance from cg
    
    #### PLACEHOLDER ####
    return F*dx

class PIDController:
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
    
class UAV:
    def __init__(self, i_moment, i_position, i_attitude,dt):
        self.i_moment = i_moment
        self.i_position = i_position
        self.i_attitude = i_attitude
        self.dt = dt
        self.position = i_position   #u, v, w
        self.attitude = i_attitude #phi theta psi
        self.euler_v  = i_moment          #p, q, r
        self.Ix = 0.0007309
        self.Iy = 0.0006644
        self.Iz = 0.0012558
        self.m = 1.28
        self.g = 9.81

    def thrust_torques(self, t_position, t_attitude):
        thrust, torques = controller.c_control(t_position, t_attitude, self.position,self.attitude,
                            self.dt)
        return thrust, torques
    
    def rk4_step(self, dt, thrust, torques):
        k1 = self.c_derivatives(self.position, self.attitude, self.euler_v, thrust, torques)
        k2 = self.c_derivatives(self.position + dt/2 * k1[:3], self.attitude + dt/2 * k1[3:6],
                                self.euler_v + dt/2 * k1[6:],thrust, torques)
        k3 = self.c_derivatives(self.position + dt/2 * k2[:3],self.attitude + dt/2 * k2[3:6],
                                self.euler_v + dt/2 * k2[6:],thrust, torques)
        k4 = self.c_derivatives(self.position + dt * k3[:3], self.attitude + dt * k3[3:6],
                                self.euler_v + dt * k3[6:],thrust, torques)
        self.position +=dt/6 *(k1[:3]+2*k2[:3]+2*k3[:3]+k4[:3])
        self.attitude +=dt/6 *(k1[3:6]+2*k2[3:6]+2*k3[3:6]+k4[3:6])
        self.euler_v +=dt/6 *(k1[6:]+2*k2[6:]+2*k3[6:]+k4[6:])

    def c_derivatives(self, x, u):
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
        Fz = F1 + F2 + F3 + F4
        L = (F2 + F3) * dy - (F1 + F4) *dy#tau phi
        M = (F1 + F3) * dx - (F2 + F4) * dx#tau theta 
        #Tってなんの関数?推進力->プロペラの推力とその半径によって回転方向にトルクを与えるものを関数Tとして->ヨーモーメントを表見してるらしい
        N = -T(F1, dx, dy) - T(F2, dx, dy) + T(F3, dx, dy) + T(F4, dx, dy) #tau psi
        
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
        
        x_dot[3] = 1/Ixx *(L + (Iyy - Izz)  * q * r) #p_dot
        x_dot[4] = 1/Iyy *(M + (Izz - Ixx) * p * r) #q_dot
        x_dot[5] = 1/Izz *(N + (Ixx- Iyy) * p * q) #r_dot
        
        x_dot[6] = p + (q * sphi + r*cphi) * sthe / cthe #phi_dot
        x_dot[7] = q * cphi - r *sphi #theta_dot
        x_dot[8] = (q * sphi + r * cphi) / cthe #psi_dot
        
        x_dot[9] = cthe*cpsi * ub + (-cphi*spsi + sphi * sthe*cpsi)*vb + \
                        (sphi*spsi + cphi*sthe*cpsi) * wb #xE_dot
        x_dot[10] = cthe * spsi * ub +(cphi *cpsi + sphi * sthe *spsi) * vb + \
                        (-sphi*cpsi + cphi*sthe*spsi) *wb #yE_dot
        x_dot[11] = -1 * (-sthe * ub + sphi*cthe * vb + cphi*cthe *wb) #hE_dot
        
        return x_dot

class Controller:
    def __init__(self,dt):
        self.dt = dt
        self.outer_pid_x = PIDController(0.35, 0.15, 8, 0.125)
        self.outer_pid_y = PIDController(0.35, 0.15, 8, 0.125)
        self.inner_pid_z = PIDController(0.35, 0.15, 8, 0.125)
        self.inner_pid_phi = PIDController(0.35, 0.15, 8, 0.125)
        self.inner_pid_psi = PIDController(0.35, 0.15, 8, 0.125)
        self.inner_pid_theta = PIDController(0.35, 0.15, 8, 0.125)

    def c_control(self, t_position, t_attitude, c_position, c_attitude,dt):
        error_x = t_position[0] - c_position[0]
        error_y = t_position[1] - c_position[1]
        t_phi = self.outer_pid_x.update(error_x, dt)
        t_theta = self.outer_pid_y.update(error_y,dt)
        t_attitude[0] = t_phi
        t_attitude[1] = t_theta
        error_z = t_position[2] - c_position[2]
        error_phi = t_attitude[0] - c_attitude[0]
        error_psi = t_attitude[1] - c_attitude[1]
        error_theta = t_attitude[2] - c_attitude[2]
        thrust = self.inner_pid_z.update(error_z,dt)
        torque_phi = self.inner_pid_phi.update(error_phi,dt)
        torque_psi = self.inner_pid_psi.update(error_psi,dt)
        torque_theta = self.inner_pid_theta.update(error_theta,dt)
        torques = [torque_phi, torque_psi, torque_theta]
        return thrust, torques

class Rotation:
    def __init__(self):
        self.l = 1

    def motor_speeds(self, thrust, torques):
        torque_x = torques[0]
        torque_y = torques[1]
        torque_z = torques[2]
        motor_torque_1 = np.clip(0.25*thrust + 0.5*torque_x /self.l - 0.5*torque_y/
                                 self.l + 0.5 *torque_z,0, np.inf)
        motor_torque_2 = np.clip(0.25*thrust - 0.5*torque_x /self.l - 0.5*torque_y/
                                 self.l - 0.5 *torque_z,0, np.inf)
        motor_torque_3 = np.clip(0.25*thrust - 0.5*torque_x /self.l + 0.5*torque_y/
                                 self.l + 0.5 *torque_z,0, np.inf)
        motor_torque_4 = np.clip(0.25*thrust + 0.5*torque_x /self.l + 0.5*torque_y/
                                 self.l - 0.5 *torque_z,0, np.inf)
        
        motor_torque_1 = max(motor_torque_1,0)
        motor_torque_2 = max(motor_torque_2,0)
        motor_torque_3 = max(motor_torque_3,0)
        motor_torque_4 = max(motor_torque_4,0)
        motor_speed_1 = np.clip(np.sqrt(motor_torque_1), 0, np.inf)
        motor_speed_2 = np.clip(np.sqrt(motor_torque_2), 0, np.inf)
        motor_speed_3 = np.clip(np.sqrt(motor_torque_3), 0, np.inf)
        motor_speed_4 = np.clip(np.sqrt(motor_torque_4), 0, np.inf)
        return motor_speed_1, motor_speed_2, motor_speed_3, motor_speed_4
    
i_position = np.array([1.0, 1.3, 1.5])
i_attitude = np.array([0.0, 0.0, 0.0])
t_position = np.array([0.0, 0.0, 0.0])
t_attitude = np.array([0.0, 0.0, 0.0])
i_moment =np.array([0.0, 0.0, 0.0]) 
time = 30
dt = 0.1
uav = UAV(i_moment, i_position, i_attitude, dt)
controller = Controller(dt)
#rotation = Rotation()
for t in np.arange(0, time, dt):
    print(i_position)
    if uav.position[2] <= 0:
        break
    thrust, torques = controller.c_control(t_position, t_attitude, uav.position,
                                            uav.attitude,dt)
    print("torques, thrust")
    print(torques)
    print(thrust)
    motor_speeds = rotation.motor_speeds(thrust, torques)
    uav.rk4_step(dt, thrust, torques)

