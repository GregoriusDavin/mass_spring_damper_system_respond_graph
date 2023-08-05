import numpy as np
import matplotlib.pyplot as plt

def mass_spring_damper(m, c, k, x0, v0, t):
    # Fungsi untuk menghitung perpindahan sistem pegas-penahan massa pada waktu t

    wn = np.sqrt(k / m)  # Frekuensi eigen
    zeta = c / (2 * np.sqrt(m * k))  # Rasio peredaman kritis

    # Mencari koefisien solusi homogen
    if zeta < 1:  # Sistem underdamped
        wd = wn * np.sqrt(1 - zeta**2)  # Frekuensi getaran redaman
        c1 = x0
        c2 = (v0 + zeta * wn * x0) / wd

        x = np.exp(-zeta * wn * t) * (c1 * np.cos(wd * t) + c2 * np.sin(wd * t))

    elif zeta == 1:  # Sistem kritis
        c1 = x0
        c2 = v0 + wn * x0

        x = (c1 + c2 * t) * np.exp(-wn * t)

    else:  # Sistem overdamped
        r1 = -wn * (zeta + np.sqrt(zeta**2 - 1))
        r2 = -wn * (zeta - np.sqrt(zeta**2 - 1))
        A = (v0 - r2 * x0) / (r1 - r2)
        B = (r1 * x0 - v0) / (r1 - r2)

        x = A * np.exp(r1 * t) + B * np.exp(r2 * t)

    return x

# Parameter sistem pegas-penahan massa
m = 1.0  # Massa (kg)
c = 0.5  # Koefisien penahan (Ns/m)
k = 2.0  # Konstanta pegas (N/m)

# Kondisi awal
x0 = 1.0  # Perpindahan awal (m)
v0 = 0.0  # Kecepatan awal (m/s)

# Waktu
dt = 0.01  # Interval waktu
t = np.arange(0.0, 10.0, dt)

# Hitung perpindahan sistem pegas-penahan massa pada waktu t
x = mass_spring_damper(m, c, k, x0, v0, t)

# Plot grafik perpindahan sistem pegas-penahan massa terhadap waktu
plt.plot(t, x)
plt.xlabel('Waktu (s)')
plt.ylabel('Perpindahan (m)')
plt.title('Grafik Respon Transien Sistem Pegas-Penahan Massa')
plt.grid(True)
plt.show()

# Mencari puncak respons
peak_value = np.max(np.abs(x))
print(f'Puncak respons transien: {peak_value:.2f} m')

# Mencari nilai setimbang
steady_state_value = np.mean(x[-int(len(x)/10):])  # Mengambil rata-rata nilai terakhir 10% data
print(f'Nilai setimbang: {steady_state_value:.2f} m')
