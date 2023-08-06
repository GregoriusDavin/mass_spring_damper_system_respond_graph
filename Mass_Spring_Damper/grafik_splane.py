import numpy as np
import matplotlib.pyplot as plt

def mass_spring_damper(m, c, k, Kp, Ki, Kd, x0, v0, t):

    dt = t[1] - t[0]  

    x = np.zeros_like(t)
    v = np.zeros_like(t)
    e = np.zeros_like(t)
    integral = 0.0

    for i in range(1, len(t)):
        e[i] = x0 - x[i-1]
        p = Kp * e[i]
        integral += e[i] * dt
        d = Kd * (e[i] - e[i-1]) / dt

        f = p + Ki * integral + d  
        a = (f - c * v[i-1] - k * x[i-1]) / m  

        v[i] = v[i-1] + a * dt  
        x[i] = x[i-1] + v[i] * dt  

    return x

m = 1.0  
c = 0.5  
k = 2.0  

Kp = 1.0  
Ki = 0.5 
Kd = 0.1  

dt = 0.01  
t = np.arange(0, 10, dt)  

x0 = 0.0  
v0 = 0.0  

x_p = mass_spring_damper(m, c, k, Kp, 0.0, 0.0, x0, v0, t)

x_pi = mass_spring_damper(m, c, k, Kp, Ki, 0.0, x0, v0, t)

x_pid = mass_spring_damper(m, c, k, Kp, Ki, Kd, x0, v0, t)

plt.plot(t, x_p, label='P')
plt.plot(t, x_pi, label='PI')
plt.plot(t, x_pid, label='PID')
plt.xlabel('Time (s)')
plt.ylabel('Displacement (m)')
plt.title('Transient Response of Mass-Spring-Damper System')
plt.legend()
plt.grid(True)
plt.show()