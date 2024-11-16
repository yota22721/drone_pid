import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

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


class QuadCopter:
    def __init__(self):
        kt = 1e-9
        self.m = 1.28
        self.g = 9.81
        self.Ixx = 0.0007309
        self.Iyy = 0.0006644
        self.Izz = 0.0012558
        self.dx = 0.125 
        self.dy = 0.125#機体中心からプロペラ中心までの距離
        self.kt = kt # proportionality constant to convert motor rotational speed into thrust (T=kt*omega^2), N/(rpm)^2
        self.b_prop = 1e-9# proportionality constant to convert motor speed to torque (torque = b*omega^2), (N*m)/(rpm)^2

    def T(T,dx,dy):
        return 0

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
        #F2 = Fthrust(x, u[1],-dx,-dy)
        #F3 = Fthrust(x, u[2],dx,-dy)
        #F4 = Fthrust(x, u[3], -dx, dy)
        F1 = kt*u[0]**2
        F2 = kt*u[1]**2
        F3 = kt*u[2]**2
        F4 = kt*u[3]**2
        Fz = F1 + F2 + F3 + F4
        L = (F2 + F3) * self.dy - (F1 + F4) *self.dy#tau phi
        M = (F1 + F3) * self.dx - (F2 + F4) * self.dx#tau theta 
        #Tってなんの関数?推進力->プロペラの推力とその半径によって回転方向にトルクを与えるものを関数Tとして->ヨーモーメントを表見してるらしい
        N = -T(F1, self.dx, self.dy) - T(F2, self.dx, self.dy) + T(F3, self.dx, self.dy) + T(F4, self.dx, self.dy) #tau psi
        #print("F1 "+ str(F1)+ " F2 "+str(F2)+" F3 "+str(F3)+" F4" +str(F4))
        #print("L "+ str(L)+ " M "+str(M)+" N "+str(N))
        
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
        
        x_dot[9] = cthe*cpsi * ub + (-cphi*spsi + sphi * sthe*cpsi)*vb + \
                        (sphi*spsi + cphi*sthe*cpsi) * wb #xE_dot
        x_dot[10] = cthe * spsi * ub +(cphi *cpsi + sphi * sthe *spsi) * vb + \
                        (-sphi*cpsi + cphi*sthe*spsi) *wb #yE_dot
        x_dot[11] = -1 * (-sthe * ub + sphi*cthe * vb + cphi*cthe *wb) #hE_dot
        
        return x_dot

    def des2speeds(self,thrust_des, tau_des):
        ''' finds speeds of motors to achieve a desired thrust and torque '''

        # Needed torque on body
        e1 = tau_des[0] * self.Ixx
        e2 = tau_des[1] * self.Iyy
        e3 = tau_des[2] * self.Izz

        #less typing
        n = self.num_motors

        # Thrust desired converted into motor speeds
        weight_speed = thrust_des / (n*self.kt)

        # Thrust differene in each motor to achieve needed torque on body
        motor_speeds = []
        motor_speeds.append(weight_speed - (e2/((n/2)*self.kt*self.L)) - (e3/(n*self.b_prop)))
        motor_speeds.append(weight_speed - (e1/((n/2)*self.kt*self.L)) + (e3/(n*self.b_prop)))
        motor_speeds.append(weight_speed + (e2/((n/2)*self.kt*self.L)) - (e3/(n*self.b_prop)))
        motor_speeds.append(weight_speed + (e1/((n/2)*self.kt*self.L)) + (e3/(n*self.b_prop)))

        # Ensure that desired thrust is within overall min and max of all motors
        thrust_all = np.array(motor_speeds) * (self.kt)
        over_max = np.argwhere(thrust_all > self.maxT)
        under_min = np.argwhere(thrust_all < self.minT)

        if over_max.size != 0:
            for i in range(over_max.size):
                motor_speeds[over_max[i][0]] = self.maxT / (self.kt)
        if under_min.size != 0:
            for i in range(under_min.size):
                motor_speeds[under_min[i][0]] = self.minT / (self.kt)
        
        self.speeds = motor_speeds