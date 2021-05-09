from math import *
import matplotlib.pyplot as plt
import numpy as np

### CONSTANTS ###
c = 2.997925 * pow(10, 10)                  # см/с
k_B = 1.380649 * pow(10, -16)               # эрг * К^-1 
m_aem = 15.994915                           # аем
m = m_aem * 1.66053906660 * pow(10, -24)    # г   
T = 1000                                    # К
lambda_S0 = 6.62 * pow(10, -14)             # см
lambda_0 = 792.23 * pow(10, -8)             # см
T_e = 10000                                 # K
N_e = pow(10, 16)                           # см^-3
#################

def Calculate_Freq_from_lambda(wave_lenght):
    return 1 / wave_lenght       # было 100000000

HWHM_Doppler = (lambda_0 / c) * sqrt(2 * k_B * log(2) * T / m)
freq_HWHM_Doppler = Calculate_Freq_from_lambda(HWHM_Doppler)


HWHM_Stark = (lambda_S0 * pow(T_e / pow(10, 4), (1/3))) * (N_e / pow(10, 16))
freq_HWHM_Stark = Calculate_Freq_from_lambda(HWHM_Stark)


freq_Voigt = 0.5 * freq_HWHM_Stark + sqrt(1/4 * freq_HWHM_Stark**2 + freq_HWHM_Doppler**2)


C2 = freq_HWHM_Stark / freq_Voigt
C3 = 2 * pow(10, -4) * freq_Voigt * (1.065 + 0.447 * C2 + 0.058 * C2**2)
C1 = (1 - C2) / C3
C4 = C2 / C3

freq_0 = Calculate_Freq_from_lambda(lambda_0)


def Calculate_b_v(freq):
    freq = freq - freq_0
    D = abs(Calculate_Freq_from_lambda(lambda_0) - freq) / (2 * freq_Voigt)
    D2 = pow(D, 2.25)

    first = (C1 * exp(-2.882 * (D**2)))
    second = (C4 / (1 + 4 * (D**2)))
    third = ((0.016 * C4 * (1 - freq_HWHM_Stark / freq_Voigt)) * (exp(-0.4 * D2) - (1 / (1 + 0.1 * D2))))

    b_v = 0.0001 * ( first + second + third)
    return b_v

x_axes = np.linspace(0, 3000000000, 100) #1500000000000000
y_axes = np.array([Calculate_b_v(x) for x in x_axes])

fig, ax = plt.subplots()                                    # будет 1 график, на нем:
ax.plot(x_axes, y_axes, color="blue", label="b_v(freq)")    # функция y1(x), синий, надпись y(x)
ax.set_xlabel("freq")                                       # подпись у горизонтальной оси х
ax.set_ylabel("I")                                          # подпись у вертикальной оси y
# ax.set_ylim(pow(10, -16) * 3.465, pow(10, -16) * 3.466)
ax.legend()                                                 # показывать условные обозначения

plt.show()   



    