import tkinter
from tkinter import *
from tkinter import ttk
from tkinter import filedialog
from tkinter import messagebox
from tkinter import PhotoImage
from PIL import Image, ImageTk
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
from matplotlib.figure import Figure
import scipy.optimize as opt
import scipy.integrate as inte
import scipy.signal as sing
import struct as st

'''Contador de Picos para el HR'''

'''Interfaz'''

'''Definimos los Frames y sus configuraciones para organizar la GUI'''

# Ventana
root = Tk()
root.title("ECG Waveform Generator - Programación científica")
root.geometry("816x616")
root.config(cursor="heart")
root.config(bg="blue")
root.config(bd=25)
root.config(relief="ridge")

# Principal
frame = Frame(root, width=480, height=320)
frame.pack(fill='both', expand=1)
frame.config(cursor="arrow")
frame.config(bg="lightblue")
frame.config(bd=25)
frame.config(relief="sunken")

# ECG
frame1 = Frame(master=root)
frame1.place(x=100, y=65)
frame1.config(bg="#F4D03F", width=250, height=50, relief=FLAT, bd=0)

# Funciones
framefun = Frame(master=root)
framefun.place(x=440, y=300)
framefun.config(bg="#A9CCE3", width=300, height=250, relief=RIDGE, bd=8)

# Parametros de entrada
frameent = Frame(master=root)
frameent.place(x=440, y=40)
frameent.config(bg="#A9CCE3", width=300, height=250, relief=RIDGE, bd=8)

# PQRST
framePQRST = Frame(master=root)
framePQRST.place(x=40, y=420)
framePQRST.config(bg="#A9CCE3", width=385, height=130, relief=RIDGE, bd=8)

# Grafica
frameGrafica = Frame(master=root)
frameGrafica.place(x=40, y=135)
frameGrafica.config(bg="#A9CCE3", width=385, height=230, relief=RIDGE, bd=8)

'''Se definen todas las funciones requeridas'''


# Función para cerrar la aplicación
def CerrarAplicacion():
    MsgBox = messagebox.askquestion('Cerrar Aplicación', '¿Está seguro que desea cerrar la aplicación?', icon='warning')
    if MsgBox == 'yes':
        root.destroy()
    else:
        messagebox.showinfo('Retornar', 'Será retornado a la aplicación')


# Función para cargar datos
def CargarDatos():
    arch = filedialog.askopenfilename(initialdir="/", title="Seleccione un archivo",
                                      filetypes=(("Archivos de texto", "*.txt"), ("Todos los archivos", "*.*")))
    fobj = open(arch, "rb")
    freadig = fobj.readline()
    fread = fobj.read()
    fobj.close()
    data = np.array(st.unpack('d' * int(len(fread) / 8), fread))
    fobj.close()
    a = str(freadig).replace(", ", " = ")
    encabezado = a.split(" = ")
    ValorFrecuenciaCardiaca.set(float(encabezado[1]))
    ValorLatidos.set(float(encabezado[3]))
    ValorFrecuenciaMuestreo.set(float(encabezado[5]))
    ValorFactorRuido.set(float(encabezado[7]))
    ValoraiP.set(float(encabezado[9].lstrip("[")))
    ValoraiQ.set(float(encabezado[10]))
    ValoraiR.set(float(encabezado[11]))
    ValoraiS.set(float(encabezado[12]))
    ValoraiT.set(float(encabezado[13].rstrip("]")))
    ValorbiP.set(float(encabezado[15].lstrip("[")))
    ValorbiQ.set(float(encabezado[16]))
    ValorbiR.set(float(encabezado[17]))
    ValorbiS.set(float(encabezado[18]))
    ValorbiT.set(float(encabezado[19].rstrip("]")))
    activos = 0
    if np.size(encabezado) > 26:
        eulerfor.set(int(encabezado[26]))
        eulerback.set(int(encabezado[27]))
        eulermod.set(int(encabezado[28]))
        rk2.set(int(encabezado[29]))
        rk4.set(int(encabezado[30]))
        if int(encabezado[26]) == 1:
            activos += 1
        if int(encabezado[27]) == 2:
            activos += 1
        if int(encabezado[28]) == 3:
            activos += 1
        if int(encabezado[29]) == 4:
            activos += 1
        if int(encabezado[30]) == 5:
            activos += 1

    y1 = np.array([])
    x1 = np.array([])
    y2 = np.array([])
    x2 = np.array([])
    y3 = np.array([])
    x3 = np.array([])
    y4 = np.array([])
    x4 = np.array([])
    y5 = np.array([])
    x5 = np.array([])
    lon = 0

    if activos != 0:
        lon = np.size(data) // (2 * activos)
    else:
        lon = np.size(data) // 2
    for i in range(0, lon):
        x1 = np.append(x1, data[i])
        y1 = np.append(y1, data[lon + i])
        if activos > 1:
            x2 = np.append(x2, data[(2 * lon) + i])
            y2 = np.append(y2, data[(3 * lon) + i])
        if activos > 2:
            x3 = np.append(x3, data[(4 * lon) + i])
            y3 = np.append(y3, data[(5 * lon) + i])
        if activos > 3:
            x4 = np.append(x4, data[(6 * lon) + i])
            y4 = np.append(y4, data[(7 * lon) + i])
        if activos > 4:
            x5 = np.append(x5, data[(8 * lon) + i])
            y5 = np.append(y5, data[(9 * lon) + i])

    gra.plot(x1, y1)
    gra.plot(x2, y2)
    gra.plot(x3, y3)
    gra.plot(x4, y4)
    gra.plot(x5, y5)
    HRnum1 = 0
    HRnum2 = 0
    HRnum3 = 0
    HRnum4 = 0
    HRnum5 = 0
    h = ActualizarH(ActualizarFs())
    T = ActualizarT(ActualizarTf(), h)
    diferencias1 = np.array([])
    PosPicos1 = sing.find_peaks(y1, height=0.02)
    for i in range(1, np.size(PosPicos1[0])):
        dif = T[PosPicos1[0][i]] - T[PosPicos1[0][i - 1]]
        diferencias1 = np.append(diferencias1, dif)
    HRnum1 = np.mean(diferencias1) * 60
    diferencias2 = np.array([])
    if np.size(y2) != 0:
        PosPicos2 = sing.find_peaks(y2, height=0.02)
        for i in range(1, np.size(PosPicos2[0])):
            dif = T[PosPicos2[0][i]] - T[PosPicos2[0][i - 1]]
            diferencias2 = np.append(diferencias2, dif)
        HRnum2 = np.mean(diferencias2) * 60
    diferencias3 = np.array([])
    if np.size(y3) != 0:
        PosPicos3 = sing.find_peaks(y3, height=0.02)
        for i in range(1, np.size(PosPicos3[0])):
            dif = T[PosPicos3[0][i]] - T[PosPicos3[0][i - 1]]
            diferencias3 = np.append(diferencias3, dif)
        HRnum3 = np.mean(diferencias3) * 60
    diferencias4 = np.array([])
    if np.size(y4) != 0:
        PosPicos4 = sing.find_peaks(y4, height=0.02)
        for i in range(1, np.size(PosPicos4[0])):
            dif = T[PosPicos4[0][i]] - T[PosPicos4[0][i - 1]]
            diferencias4 = np.append(diferencias4, dif)
        HRnum4 = np.mean(diferencias4) * 60
    diferencias5 = np.array([])
    if np.size(y5) != 0:
        PosPicos5 = sing.find_peaks(y5, height=0.02)
        for i in range(1, np.size(PosPicos5[0])):
            dif = T[PosPicos5[0][i]] - T[PosPicos5[0][i - 1]]
            diferencias5 = np.append(diferencias5, dif)
        HRnum5 = np.mean(diferencias5) * 60
    if activos != 0:
        prom = (HRnum5 + HRnum4 + HRnum3 + HRnum2 + HRnum1) / activos
        HRtxt.set(prom)
        messagebox.showinfo('Advertencia', 'Se ha calculado el HR promedio para evitar interferencia de datos')
    if activos == 0 and np.size(y1) != 0:
        HRtxt.set(HRnum1)
        messagebox.showinfo('Advertencia', 'Se ha calculado el HR promedio para evitar interferencia de datos')
    if np.size(y1) == 0 and np.size(y2) == 0 and np.size(y3) == 0 and np.size(y4) == 0 and np.size(y5) == 0:
        messagebox.showinfo('Error al graficar',
                            'El archivo no contiene valores del EGG')


# Función para exportar datos
def ExportarDatos():
    enunciado = "Frecuencia cardiaca media = " + str(
        ValorFrecuenciaCardiaca.get()) + ", Número de latidos = " + str(
        ValorLatidos.get()) + ", Frecuencia de muestreo = " + str(
        ValorFrecuenciaMuestreo.get()) + ",  Factor de ruido = " + str(ValorFactorRuido.get()) + ", ai = " + str(
        ai) + ", bi = " + str(bi) + ", theta_i = " + str(theta_i) + ", " + str(eulerfor.get()) + ", " + str(
        eulerback.get()) + ", " + str(eulermod.get()) + ", " + str(rk2.get()) + ", " + str(rk4.get()) + ", " + "\n"
    guard = filedialog.asksaveasfilename(initialdir="/", title="Exportar datos", defaultextension=".txt",
                                         filetypes=(("Archivos de texto", ".txt"), ("Todos los arcivos", ".")))
    archivo = open(guard, "w")
    archivo.write(enunciado)
    archivo.close()
    y1 = np.array([])
    x1 = np.array([])
    y2 = np.array([])
    x2 = np.array([])
    y3 = np.array([])
    x3 = np.array([])
    y4 = np.array([])
    x4 = np.array([])
    y5 = np.array([])
    x5 = np.array([])
    si = FALSE
    if eulerfor.get() == 1:
        x1, y1 = EulerFor()
        xPack = st.pack("d" * np.size(x1), *x1)
        yPack = st.pack("d" * np.size(y1), *y1)
        archivo = open(guard, "ab")
        archivo.write(xPack)
        archivo.write(yPack)
        archivo.close()
        archivo.close()
        si = TRUE
    if eulerback.get() == 2:
        x2, y2 = EulerBack()
        xPack = st.pack("d" * np.size(x2), *x2)
        yPack = st.pack("d" * np.size(y2), *y2)
        archivo = open(guard, "ab")
        archivo.write(xPack)
        archivo.write(yPack)
        archivo.close()
        archivo.close()
        si = TRUE
    if eulermod.get() == 3:
        x3, y3 = EulerModificado()
        xPack = st.pack("d" * np.size(x3), *x3)
        yPack = st.pack("d" * np.size(y3), *y3)
        archivo = open(guard, "ab")
        archivo.write(xPack)
        archivo.write(yPack)
        archivo.close()
        archivo.close()
        si = TRUE
    if rk2.get() == 4:
        x4, y4 = RK2()
        xPack = st.pack("d" * np.size(x4), *x4)
        yPack = st.pack("d" * np.size(y4), *y4)
        archivo = open(guard, "ab")
        archivo.write(xPack)
        archivo.write(yPack)
        archivo.close()
        archivo.close()
        si = TRUE
    if rk4.get() == 5:
        x5, y5 = RK4()
        xPack = st.pack("d" * np.size(x5), *x5)
        yPack = st.pack("d" * np.size(y5), *y5)
        archivo = open(guard, "ab")
        archivo.write(xPack)
        archivo.write(yPack)
        archivo.close()
        archivo.close()
        si = TRUE
    if si == TRUE:
        messagebox.showinfo('Exportar datos', 'El archivo fue creado con exito')
    else:
        messagebox.showinfo('Error', 'No hay datos registrados')


# Función para hallar el HR
def HRcalculo():
    y1 = np.array([])
    y2 = np.array([])
    y3 = np.array([])
    y4 = np.array([])
    y5 = np.array([])
    x = np.array([])
    HRnum1 = 0
    HRnum2 = 0
    HRnum3 = 0
    HRnum4 = 0
    HRnum5 = 0
    activos = 0
    if eulerfor.get() == 1:
        x, y1 = EulerFor()
        PosPicos = sing.find_peaks(y1, height=0.02)
        diferencias = np.array([])
        for i in range(1, np.size(PosPicos[0])):
            dif = T[PosPicos[0][i]] - T[PosPicos[0][i - 1]]
            diferencias = np.append(diferencias, dif)
        HRnum1 = np.mean(diferencias) * 60
        activos += 1
    if eulerback.get() == 2:
        x, y2 = EulerBack()
        PosPicos = sing.find_peaks(y2, height=0.02)
        diferencias = np.array([])
        for i in range(1, np.size(PosPicos[0])):
            dif = T[PosPicos[0][i]] - T[PosPicos[0][i - 1]]
            diferencias = np.append(diferencias, dif)
        HRnum2 = np.mean(diferencias) * 60
        activos += 1
    if eulermod.get() == 3:
        x, y3 = EulerModificado()
        PosPicos = sing.find_peaks(y3, height=0.02)
        diferencias = np.array([])
        for i in range(1, np.size(PosPicos[0])):
            dif = T[PosPicos[0][i]] - T[PosPicos[0][i - 1]]
            diferencias = np.append(diferencias, dif)
        HRnum3 = np.mean(diferencias) * 60
        activos += 1
    if rk2.get() == 4:
        x, y4 = RK2()
        PosPicos = sing.find_peaks(y4, height=0.02)
        diferencias = np.array([])
        for i in range(1, np.size(PosPicos[0])):
            dif = T[PosPicos[0][i]] - T[PosPicos[0][i - 1]]
            diferencias = np.append(diferencias, dif)
        HRnum4 = np.mean(diferencias) * 60
        activos += 1
    if rk4.get() == 5:
        x, y5 = RK4()
        PosPicos = sing.find_peaks(y5, height=0.02)
        diferencias = np.array([])
        for i in range(1, np.size(PosPicos[0])):
            dif = T[PosPicos[0][i]] - T[PosPicos[0][i - 1]]
            diferencias = np.append(diferencias, dif)
        HRnum5 = np.mean(diferencias) * 60
        activos += 1
    if activos == 0:
        messagebox.showinfo('Error', 'No hay datos registrados')
    if (HRnum1 != 0 and HRnum2 + HRnum3 + HRnum4 + HRnum5 > 0) or (
            HRnum2 != 0 and HRnum1 + HRnum3 + HRnum4 + HRnum5 > 0) or (
            HRnum3 != 0 and HRnum2 + HRnum1 + HRnum4 + HRnum5 > 0) or (
            HRnum4 != 0 and HRnum2 + HRnum3 + HRnum1 + HRnum5 > 0) or (
            HRnum5 != 0 and HRnum2 + HRnum3 + HRnum4 + HRnum1 > 0):
        MsgBox = messagebox.askquestion('Advertencia',
                                        'Ha seleccionado más de un método.\n' + "¿Desea hallar el HR promedio de los métodos seleccionados?",
                                        icon='warning')
        if MsgBox == 'yes':
            suma = HRnum1 + HRnum2 + HRnum3 + HRnum4 + HRnum5
            prom = suma / activos
            HRtxt.set(prom)
        else:
            messagebox.showinfo('Advertencia', 'Seleccione un solo un método')
    else:
        suma = HRnum1 + HRnum2 + HRnum3 + HRnum4 + HRnum5
        HRtxt.set(suma)


'''Sistema de Ecuaciones ECG'''

# Condiciones Iniciales:
ValoraiP = DoubleVar(value=1.2)
ValoraiQ = DoubleVar(value=-5.0)
ValoraiR = DoubleVar(value=30.0)
ValoraiS = DoubleVar(value=-7.5)
ValoraiT = DoubleVar(value=0.75)
ValorbiP = DoubleVar(value=0.25)
ValorbiQ = DoubleVar(value=0.1)
ValorbiR = DoubleVar(value=0.1)
ValorbiS = DoubleVar(value=0.1)
ValorbiT = DoubleVar(value=0.4)
ValorFrecuenciaCardiaca = DoubleVar(value="80.0")
ValorLatidos = DoubleVar(value=30.0)
ValorFrecuenciaMuestreo = DoubleVar(value="200.0")  # Linea149
ValorFactorRuido = DoubleVar(value=0.0025)

X0 = 1
Y0 = 0
Z0 = 0.03
theta_i = [(-1 / 3) * np.pi, (-1 / 12) * np.pi, 0, (1 / 12) * np.pi, (1 / 2) * np.pi]
To = 0


def Actualizarai():
    ai = [ValoraiP.get(), ValoraiQ.get(), ValoraiR.get(), ValoraiS.get(), ValoraiT.get()]
    return ai


def Actualizarbi():
    bi = [ValorbiP.get(), ValorbiQ.get(), ValorbiR.get(), ValorbiS.get(), ValorbiT.get()]
    return bi


def ActualizarHR():
    HR = ValorFrecuenciaCardiaca.get()
    return HR


def ActualizarPromedioHR(HR):
    promedioHR = 60 / HR
    return promedioHR


def ActualizarFs():
    fs = ValorFrecuenciaMuestreo.get()
    return fs


def ActualizarH(fs):
    h = 1 / fs
    return h


def ActualizarTf():
    Tf = ValorLatidos.get()
    return Tf


def ActualizarT(Tf, h):
    T = np.arange(To, Tf + h, h)
    return T


def ActualizarRR(promedioHR):
    rr = np.random.normal(promedioHR, promedioHR * 0.05)
    return rr


def ActualizarRuido(T):
    ruido = np.random.normal(0, ValorFactorRuido.get(), np.size(T))
    return ruido


T = ActualizarT(ActualizarTf(), ActualizarH(ActualizarFs()))
ai = Actualizarai()
bi = Actualizarbi()


# Ecuaciones para (x, y, z)

def F1(y1, y2, rr):
    alpha = 1.0 - np.sqrt(y1 ** 2 + y2 ** 2)
    omega = 2.0 * np.pi / rr
    return (alpha * y1) - (omega * y2)


def F2(y1, y2, rr):
    alpha = 1.0 - np.sqrt(y1 ** 2 + y2 ** 2)
    omega = 2.0 * np.pi / rr
    return (alpha * y2) + (omega * y1)


def F3(y1, y2, y3, T):
    theta = np.arctan2(y2, y1)
    sumaZ = 0
    z0 = (0.15 * np.sin(2.0 * np.pi * 0.25 * T))
    ai = Actualizarai()
    bi = Actualizarbi()
    for x in range(np.size(ai)):
        deltatt = np.fmod(theta - theta_i[x], 2.0 * np.pi)
        sumaZ = ai[x] * deltatt * np.exp(-(deltatt ** 2) / (2 * bi[x] ** 2)) + sumaZ

    return -sumaZ - (y3 - (0.0015 * np.sin(2 * np.pi * T * 0.25)))


def EulerFor():
    h = ActualizarH(ActualizarFs())
    T = ActualizarT(ActualizarTf(), h)
    rr = ActualizarRR(ActualizarPromedioHR(ActualizarHR()))
    ruido = ActualizarRuido(T)
    Y1EulerFor = np.zeros(len(T))
    Y2EulerFor = np.zeros(len(T))
    Y3EulerFor = np.zeros(len(T))
    Y1EulerFor[0] = X0
    Y2EulerFor[0] = Y0
    Y3EulerFor[0] = Z0
    for i in range(1, len(T)):
        # Euler hacia adelante
        Y1EulerFor[i] = Y1EulerFor[i - 1] + h * F1(Y1EulerFor[i - 1],
                                                   Y2EulerFor[i - 1], rr)
        Y2EulerFor[i] = Y2EulerFor[i - 1] + h * F2(Y1EulerFor[i - 1],
                                                   Y2EulerFor[i - 1], rr)
        Y3EulerFor[i] = Y3EulerFor[i - 1] + h * F3(Y1EulerFor[i - 1],
                                                   Y2EulerFor[i - 1], Y3EulerFor[i - 1], T[i])
    return T, Y3EulerFor + ruido


def RK2():
    h = ActualizarH(ActualizarFs())
    T = ActualizarT(ActualizarTf(), h)
    rr = ActualizarRR(ActualizarPromedioHR(ActualizarHR()))
    ruido = ActualizarRuido(T)
    Y1RK2 = np.zeros(len(T))
    Y2RK2 = np.zeros(len(T))
    Y3RK2 = np.zeros(len(T))
    Y1RK2[0] = X0
    Y2RK2[0] = Y0
    Y3RK2[0] = Z0
    for i in range(1, len(T)):
        k11 = F1(Y1RK2[i - 1], Y2RK2[i - 1], rr)
        k21 = F2(Y1RK2[i - 1], Y2RK2[i - 1], rr)
        k31 = F3(Y1RK2[i - 1], Y2RK2[i - 1], Y3RK2[i - 1], T[i])

        k12 = F1(Y1RK2[i - 1] + k11 * h, Y2RK2[i - 1] + k21 * h, rr)
        k22 = F2(Y1RK2[i - 1] + k11 * h, Y2RK2[i - 1] + k21 * h, rr)
        k32 = F3(Y1RK2[i - 1] + k11 * h, Y2RK2[i - 1] + k21 * h, Y3RK2[i - 1] + k31 * h, T[i])

        Y1RK2[i] = Y1RK2[i - 1] + (h / 2) * (k11 + k12)
        Y2RK2[i] = Y2RK2[i - 1] + (h / 2) * (k21 + k22)
        Y3RK2[i] = Y3RK2[i - 1] + (h / 2) * (k31 + k32)
    return T, Y3RK2 + ruido


def RK4():
    h = ActualizarH(ActualizarFs())
    T = ActualizarT(ActualizarTf(), h)
    rr = ActualizarRR(ActualizarPromedioHR(ActualizarHR()))
    ruido = ActualizarRuido(T)
    Y1RK4 = np.zeros(len(T))
    Y2RK4 = np.zeros(len(T))
    Y3RK4 = np.zeros(len(T))
    Y1RK4[0] = X0
    Y2RK4[0] = Y0
    Y3RK4[0] = Z0
    for i in range(1, len(T)):
        k11 = F1(Y1RK4[i - 1], Y2RK4[i - 1], rr)
        k21 = F2(Y1RK4[i - 1], Y2RK4[i - 1], rr)
        k31 = F3(Y1RK4[i - 1], Y2RK4[i - 1], Y3RK4[i - 1], T[i])

        k12 = F1(Y1RK4[i - 1] + 0.5 * k11 * h, Y2RK4[i - 1] + 0.5 * k21 * h, rr)
        k22 = F2(Y1RK4[i - 1] + 0.5 * k11 * h, Y2RK4[i - 1] + 0.5 * k21 * h, rr)
        k32 = F3(Y1RK4[i - 1] + 0.5 * k11 * h, Y2RK4[i - 1] + 0.5 * k21 * h, Y3RK4[i - 1] + 0.5 * k31 * h, T[i])

        k13 = F1(Y1RK4[i - 1] + 0.5 * k12 * h, Y2RK4[i - 1] + 0.5 * k22 * h, rr)
        k23 = F2(Y1RK4[i - 1] + 0.5 * k12 * h, Y2RK4[i - 1] + 0.5 * k22 * h, rr)
        k33 = F3(Y1RK4[i - 1] + 0.5 * k12 * h, Y2RK4[i - 1] + 0.5 * k22 * h, Y3RK4[i - 1] + 0.5 * k32 * h, T[i])

        k14 = F1(Y1RK4[i - 1] + k13 * h, Y2RK4[i - 1] + k23 * h, rr)
        k24 = F2(Y1RK4[i - 1] + k13 * h, Y2RK4[i - 1] + k23 * h, rr)
        k34 = F3(Y1RK4[i - 1] + k13 * h, Y2RK4[i - 1] + k23 * h, Y3RK4[i - 1] + k33 * h, T[i])

        Y1RK4[i] = Y1RK4[i - 1] + (h / 6) * (k11 + k12 + k13 + k14)
        Y2RK4[i] = Y2RK4[i - 1] + (h / 6) * (k21 + k22 + k23 + k24)
        Y3RK4[i] = Y3RK4[i - 1] + (h / 6) * (k31 + k32 + k33 + k34)

    return T, Y3RK4 + ruido


# Euler Modificado

def EulerModraiz(g, y1t1, y2t1, y3t1, h, tiempo1, rr):
    return [y1t1 + (h / 2.0) * (F1(y1t1, y2t1, rr) + F1(g[0], g[1], rr)) - g[0],
            y2t1 + (h / 2.0) * (F2(y1t1, y2t1, rr) + F2(g[0], g[1], rr)) - g[1],
            y3t1 + (h / 2.0) * (F3(y1t1, y2t1, y3t1, tiempo1) + F3(g[0], g[1], g[2], tiempo1)) - g[2]]


def EulerModificado():
    h = ActualizarH(ActualizarFs())
    T = ActualizarT(ActualizarTf(), h)
    rr = ActualizarRR(ActualizarPromedioHR(ActualizarHR()))
    ruido = ActualizarRuido(T)
    Y1EulerModraiz = np.zeros(len(T))
    Y2EulerModraiz = np.zeros(len(T))
    Y3EulerModraiz = np.zeros(len(T))
    Y1EulerModraiz[0] = X0
    Y2EulerModraiz[0] = Y0
    Y3EulerModraiz[0] = Z0
    for i in range(1, len(T)):
        solucion = opt.fsolve(EulerModraiz,
                              np.array([Y1EulerModraiz[i - 1],
                                        Y2EulerModraiz[i - 1],
                                        Y3EulerModraiz[i - 1]]),
                              (Y1EulerModraiz[i - 1], Y2EulerModraiz[i - 1], Y3EulerModraiz[i - 1], h, T[i - 1], rr),
                              xtol=10 ** -5)
        Y1EulerModraiz[i] = solucion[0]
        Y2EulerModraiz[i] = solucion[1]
        Y3EulerModraiz[i] = solucion[2]

    return T, Y3EulerModraiz + ruido


variabledetestable = EulerModificado()


# Euler Back

def EulerBackraiz(g, y1t1, y2t1, y3t1, h, tiempo1, rr):
    return [y1t1 + h * F1(g[0], g[1], rr) - g[0],
            y2t1 + h * F2(g[0], g[1], rr) - g[1],
            y3t1 + h * F3(g[0], g[1], g[2], tiempo1) - g[2]]


def EulerBack():
    h = ActualizarH(ActualizarFs())
    T = ActualizarT(ActualizarTf(), h)
    rr = ActualizarRR(ActualizarPromedioHR(ActualizarHR()))
    ruido = ActualizarRuido(T)
    Y1EulerBackraiz = np.zeros(len(T))
    Y2EulerBackraiz = np.zeros(len(T))
    Y3EulerBackraiz = np.zeros(len(T))
    Y1EulerBackraiz[0] = X0
    Y2EulerBackraiz[0] = Y0
    Y3EulerBackraiz[0] = Z0
    for i in range(1, len(T)):
        solucion = opt.fsolve(EulerBackraiz,
                              np.array([Y1EulerBackraiz[i - 1],
                                        Y2EulerBackraiz[i - 1],
                                        Y3EulerBackraiz[i - 1]]),
                              (Y1EulerBackraiz[i - 1], Y2EulerBackraiz[i - 1], Y3EulerBackraiz[i - 1], h, T[i], rr),
                              xtol=10 ** -2)
        Y1EulerBackraiz[i] = solucion[0]
        Y2EulerBackraiz[i] = solucion[1]
        Y3EulerBackraiz[i] = solucion[2]

    return T, Y3EulerBackraiz + ruido


F = EulerBack()
variablequeodio = F[1]

# Formato de los botones a utilizar
Style = ttk.Style()
Style.configure('1.TButton', font=('Ink Free', 5, 'bold italic'), foreground='red', background='red', padding=2)
Style.configure('2.TButton', font=('Ink Free', 10, 'bold italic'), foreground='black', background='blue',
                padding=2)
Style.configure('3.TButton', font=('Ink Free', 12, 'bold italic'), foreground='black', background='#0000FF',
                padding=2)

# Para el botón de cerrar
Style.map("1.TButton",
          foreground=[('pressed', 'red'), ('active', 'red')],
          background=[('pressed', '!disabled', 'red'), ('active', 'red')])

# Para los botones de Exportar datos y Cargar datos
Style.map("2.TButton",
          foreground=[('pressed', 'yellow'), ('active', '#34495E')],
          background=[('pressed', '!disabled', 'black'), ('active', 'white')])

# Para el botón hallar HR
Style.map("3.TButton",
          foreground=[('pressed', 'blue'), ('active', 'blue')],
          background=[('pressed', '!disabled', 'black'), ('active', 'white')])

# Creación de botones
BotonX = ttk.Button(master=root, text="X", style="1.TButton", command=CerrarAplicacion).place(x=15, y=15)
BotonExportarDatos = ttk.Button(master=root, text="Exportar datos", style="2.TButton", command=ExportarDatos).place(
    x=110, y=29)
BotonCargarDatos = ttk.Button(master=root, text="Cargar datos", style="2.TButton", command=CargarDatos).place(
    x=240, y=29)
BotonHallarHR = ttk.Button(master=root, text="Hallar HR", style="3.TButton", command=HRcalculo).place(x=110, y=380)

'''ECG'''
lbl_titulo = Label(master=frame1, bg="#F4D03F", font=('Ink Free', 15, 'bold italic'), text=f"Señal de ECG").place(
    x=55, y=10)

'''HR'''
HRtxt = StringVar(value="0")
lbl_x = Label(master=root, textvariable=HRtxt, width=15, bg='#B2BABB', fg='blue').place(x=250, y=387)

'''Graficador de funciones'''

fig = plt.figure(figsize=(3.5, 1.9), dpi=100)
gra = fig.add_subplot(111)  # subplot(filas, columnas, item)

plt.close()

canvas = FigureCanvasTkAgg(fig, master=frameGrafica)

NavigationToolbar2Tk(canvas, frameGrafica)
canvas._tkcanvas.pack()

canvas.draw()
canvas.get_tk_widget().pack()


def fun():
    if rk4.get() == 10:
        gra.clear()
    if rk2.get() == 9:
        gra.clear()
    if eulermod.get() == 8:
        gra.clear()
    if eulerback.get() == 7:
        gra.clear()
    if eulerfor.get() == 6:
        gra.clear()
    if eulerfor.get() == 1:
        x, y = EulerFor()
        gra.plot(x, y, "#FF8C00")
    if eulerback.get() == 2:
        x, y = EulerBack()
        gra.plot(x, y, "#CD1076")
    if eulermod.get() == 3:
        x, y = EulerModificado()
        gra.plot(x, y, "#EE2C2C")
    if rk2.get() == 4:
        x, y = RK2()
        gra.plot(x, y, "#00C957")
    if rk4.get() == 5:
        x, y = RK4()
        gra.plot(x, y, "#00688B")
    if eulerfor.get() == 6 and eulerback.get() == 7 and eulermod.get() == 8 and rk2.get() == 9 and rk4.get() == 10:
        messagebox.showinfo('Error', 'No ha seleccionado ningún método')


'''Métodos de solución'''
lbl_gra = Label(master=framefun, bg='#A9CCE3', font=('Ink Free', 15, 'bold italic'),
                text=f"Método de Solución ED").place(x=30, y=15)

# Se definen las variables referentes a los métodos de solución
opcion = IntVar()
eulerfor = IntVar(value=6)
eulerback = IntVar(value=7)
eulermod = IntVar(value=8)
rk2 = IntVar(value=9)
rk4 = IntVar(value=10)

# Se crean los botones de caracter checkbutton para que esto permita elegir más de una opción
eulerFor = Checkbutton(master=framefun, text="Euler adelante", variable=eulerfor, bg='#A9CCE3',
                       font=('Ink Free', 15, 'bold italic'), onvalue=1, offvalue=6, command=fun).place(x=30, y=60)
eulerBack = Checkbutton(master=framefun, text='Euler atrás', variable=eulerback, bg='#A9CCE3',
                        font=('Ink Free', 15, 'bold italic'), onvalue=2, offvalue=7, command=fun).place(x=30, y=90)
EulerMod = Checkbutton(master=framefun, text='Euler modificado', variable=eulermod, bg='#A9CCE3',
                       font=('Ink Free', 15, 'bold italic'), onvalue=3, offvalue=8, command=fun).place(x=30, y=120)
RungeKutta2 = Checkbutton(master=framefun, text='Runge-Kutta 2', variable=rk2, bg='#A9CCE3',
                          font=('Ink Free', 15, 'bold italic'), onvalue=4, offvalue=9, command=fun).place(x=30, y=150)
RungeKutta4 = Checkbutton(master=framefun, text='Runge-Kutta 4', variable=rk4, bg='#A9CCE3',
                          font=('Ink Free', 15, 'bold italic'), onvalue=5, offvalue=10, command=fun).place(x=30, y=180)

# Botón actualizar
" Actualiza los valores de ai y bi, mostrando con esto, una nueva gráfica con los nuevos valores"

Style.map("5.TButton",
          foreground=[('pressed', 'blue'), ('active', 'blue')],
          background=[('pressed', '!disabled', 'black'), ('active', 'white')])

img = Image.open('ECG-GRUPO_30/actualizaaaa.png')
img = img.resize((20, 20), Image.ANTIALIAS)  # Redimension (Alto, Ancho)
img = ImageTk.PhotoImage(img)
BotonActualizar = ttk.Button(master=framePQRST, style="5.TButton", command=fun, image=img).place(x=5, y=5)

'''Parametros de entrada'''
lbl_ent = Label(master=frameent, bg='#A9CCE3', font=('Ink Free', 15, 'bold italic'), text=f"Paramétros").place(
    x=90, y=10)
lbl_frecuenciacardiaca = Label(master=frameent, bg='#A9CCE3', font=('Ink Free', 12, 'bold italic'),
                               text=f"Frecuencia Cardiaca").place(x=15, y=53)
lbl_latidos = Label(master=frameent, bg='#A9CCE3', font=('Ink Free', 12, 'bold italic'),
                    text=f"# de latidos").place(x=45, y=93)
lbl_frecuenciamuestreo = Label(master=frameent, bg='#A9CCE3', font=('Ink Free', 12, 'bold italic'),
                               text=f"Frecuencia Muestreo").place(x=15, y=133)
lbl_factorruido = Label(master=frameent, bg='#A9CCE3', font=('Ink Free', 12, 'bold italic'),
                        text=f"Factor de Ruido").place(x=35, y=173)

ent_frecuenciacardiaca = Entry(master=frameent, textvariable=ValorFrecuenciaCardiaca, width=12).place(x=190, y=60)
ent_latidos = Entry(master=frameent, textvariable=ValorLatidos, width=12).place(x=190, y=100)
ent_frecuenciamuestreo = Entry(master=frameent, textvariable=ValorFrecuenciaMuestreo, width=12).place(x=190, y=140)
ent_factorRuido = Entry(master=frameent, textvariable=ValorFactorRuido, width=12).place(x=190, y=180)

'''PQRST'''
lbl_PQRST = Label(master=framePQRST, bg='#A9CCE3', font=('Ink Free', 15, 'bold italic'),
                  text=f"P        Q        R        S        T").place(x=50, y=7)
lbl_ai = Label(master=framePQRST, bg='#A9CCE3', font=('Ink Free', 12, 'bold italic'), text=f"ai").place(x=10, y=40)
lbl_bi = Label(master=framePQRST, bg='#A9CCE3', font=('Ink Free', 12, 'bold italic'), text=f"bi").place(x=10, y=76)

ent_aiP = Entry(master=framePQRST, textvariable=ValoraiP, width=6).place(x=40, y=43)
ent_aiQ = Entry(master=framePQRST, textvariable=ValoraiQ, width=6).place(x=110, y=43)
ent_aiR = Entry(master=framePQRST, textvariable=ValoraiR, width=6).place(x=180, y=43)
ent_aiS = Entry(master=framePQRST, textvariable=ValoraiS, width=6).place(x=250, y=43)
ent_aiT = Entry(master=framePQRST, textvariable=ValoraiT, width=6).place(x=320, y=43)
ent_biP = Entry(master=framePQRST, textvariable=ValorbiP, width=6).place(x=40, y=79)
ent_biQ = Entry(master=framePQRST, textvariable=ValorbiQ, width=6).place(x=110, y=79)
ent_biR = Entry(master=framePQRST, textvariable=ValorbiR, width=6).place(x=180, y=79)
ent_biS = Entry(master=framePQRST, textvariable=ValorbiS, width=6).place(x=250, y=79)
ent_biT = Entry(master=framePQRST, textvariable=ValorbiT, width=6).place(x=320, y=79)

root.mainloop()
