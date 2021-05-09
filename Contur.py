from math import *
import matplotlib.pyplot as plt
import numpy as np

### CONSTANTS ###
c = 299792458                           # m/s
k_B = 8.617333262 * pow(10, -5)         # эВ * К^-1 
m = 15.994915                           # аем   
T = 1000                                # К
lambda_S0_CONST = 6.62 * pow(10, -6)    # A
lambda_0_CONST = 792.23                 # A
T_e = 10000                             # K
N_e = pow(10, 16)                       # см^-3
#################

class Contur:

    def Calculate_Freq_from_lambda(self, wave_lenght):
        return pow(10, 8) / wave_lenght

    # конструктор
    def __init__(self, input_lambda_0, input_lambda_S0):
        if(input_lambda_0 != None):
            self.lambda_0 = input_lambda_0  
        else:
            self.lambda_0 = lambda_0_CONST

        if(input_lambda_S0 != None):
            self.lambda_S0 = input_lambda_S0
        else:
            self.lambda_S0 = lambda_S0_CONST


        self.HWHM_Doppler = (self.lambda_0 / c) * sqrt(2 * k_B * log(2) * T / m)
        self.freq_HWHM_Doppler = self.Calculate_Freq_from_lambda(self.HWHM_Doppler)


        self.HWHM_Stark = (self.lambda_S0 * pow(T_e / pow(10, 4), (1/3))) * (N_e / pow(10, 16))
        self.freq_HWHM_Stark = self.Calculate_Freq_from_lambda(self.HWHM_Stark)


        self.freq_Voigt = 0.5 * self.freq_HWHM_Stark + sqrt(1/4 * self.freq_HWHM_Stark**2 + self.freq_HWHM_Doppler**2)


        self.C2 = self.freq_HWHM_Stark / self.freq_Voigt
        self.C3 = 2 * pow(10, -4) * self.freq_Voigt * (1.065 + 0.447 * self.C2 + 0.058 * self.C2**2)
        self.C1 = (1 - self.C2) / self.C3
        self.C4 = self.C2 / self.C3

        self.freq_0 = self.Calculate_Freq_from_lambda(self.lambda_0)


    def Calculate_b_v(self, freq):
        freq = freq - self.freq_0
        D = abs(self.Calculate_Freq_from_lambda(self.lambda_0) - freq) / (2 * self.freq_Voigt)
        D2 = pow(D, 2.25)

        first = (self.C1 * exp(-2.882 * (D**2)))
        second = (self.C4 / (1 + 4 * (D**2)))
        third = ((0.016 * self.C4 * (1 - self.freq_HWHM_Stark / self.freq_Voigt)) * (exp(-0.4 * D2) - (1 / (1 + 0.1 * D2))))

        b_v = 0.0001 * ( first + second + third)
        return b_v


    #функция "Нарисовать график" - считает и рисует
    def Draw_Plot(self, From, To, Number_of_points = 100):
        x_axes = np.linspace(From, To, Number_of_points) #1500000000000000
        y_axes = np.array([self.Calculate_b_v(x) for x in x_axes])

        fig, ax = plt.subplots()                                    # будет 1 график, на нем:
        ax.plot(x_axes, y_axes, color="blue", label="b_v(freq)")    # функция y1(x), синий, надпись y(x)
        ax.set_xlabel("freq")                                       # подпись у горизонтальной оси х
        ax.set_ylabel("I")                                          # подпись у вертикальной оси y
        # ax.set_ylim(pow(10, -16) * 3.465, pow(10, -16) * 3.466)
        ax.legend()                                                 # показывать условные обозначения

        plt.show()   
        