import numpy as np
import matplotlib.pyplot as plt
#import plotly.express as px

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

    def c_derivatives(self, position, attitude, euler_v, thrust, torques):
        phi, theta, psi = euler_v
        p, q, r = euler_v
        dphi = p + (q * np.sin(phi) + r * np.cos(phi)) * np.tan(theta)
        dtheta = q *np.cos(phi) - r*np.sin(phi)
        dpsi = (q * np.sin(phi) + r * np.cos(phi))/np.cos(theta)
        dp = (torques[0] + (self.Iz - self.Iy) * q * r)/self.Ix
        dq = (torques[1] + (self.Ix - self.Iz) * p * r)/self.Iy
        dr = (torques[2] + (self.Iy - self.Ix) * p * q)/self.Iz
        du = -self.g + (thrust/self.m) *(np.cos(phi) * np.sin(theta) * np.cos(psi)
                                         + np.sin(phi) * np.sin(psi))
        dv = -self.g + (thrust/self.m) *(np.cos(phi) * np.sin(theta) * np.sin(psi)
                                         - np.sin(phi) * np.cos(psi))
        dw = (thrust/self.m) * (np.cos(phi) * np.cos(theta)) - self.g
        return np.array([dphi, dtheta, dpsi, dp, dq, dr,du, dv, dw])

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
        #thrust_2 = self.inner_pid_z.update(error_y,dt)
        #thrust_3 = self.inner_pid_z.update(error_z,dt)
        torque_phi = self.inner_pid_phi.update(error_phi,dt)
        torque_psi = self.inner_pid_psi.update(error_psi,dt)
        torque_theta = self.inner_pid_theta.update(error_theta,dt)
        #thrust = [thrust_1, thurst_2, thurst_3, thrust_4]
        torques = [torque_phi, torque_psi, torque_theta]
        return thrust, torques

class Rotation:
    def __init__(self):
        self.l = 0.125
    def motor_speeds(self, thrust, torques):
        torque_x = torques[0]
        torque_y = torques[1]
        torque_z = torques[2]
        motor_torque_1 = np.clip(0.25*thrust + 0.5*torque_x /self.l - 0.5*torque_y/
                                 self.l + 0.5 *torque_z, 0, np.inf)
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
dt = 0.02
uav = UAV(i_moment, i_position, i_attitude, dt)
controller = Controller(dt)
rotation = Rotation()
#t_position = i_position
for t in np.arange(0, time, dt):
    #print(i_position)
    #print(t_position)
    #print(t_attitude)
    print(uav.position)
    if uav.position[2] <= 0:
        break
    thrust, torques = controller.c_control(t_position, t_attitude, uav.position,
                                            uav.attitude,dt)
    print("torques, thrust, position")
    #print(torques)
    #print(thrust)
    motor_speeds = rotation.motor_speeds(thrust, torques)
    print(motor_speeds)
    torques = []
    uav.rk4_step(dt, thrust, torques)

