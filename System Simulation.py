import numpy as np
from bokeh.plotting import figure, show,output_file
from bokeh.models import Panel, Tabs

class System:
    def __init__(self, n, m, a, b, inp_type, t):
        self.n = n + 1
        self.m = m + 1
        self.end_t = t
        self.k = 1000
        self.h = (self.end_t - 0) / self.k
        self.a = np.array(a) / a[0]
        self.b = np.array(b) / a[0]
        self.inp_type = inp_type
        self.time = np.linspace(0, self.end_t, self.k + 1)
        self.u = np.ones(self.k + 1)
        self.Beta = np.zeros(self.n)
        self.A = np.zeros((self.n - 1, self.n - 1))
        self.B = np.zeros(self.n - 1)
        self.C = np.zeros(self.n - 1)
        self.D = 0
        self.Xs = np.zeros((self.k + 1, self.n - 1))
        self.X_temp = np.zeros(self.n - 1)
        self.y = np.zeros((self.k + 1))

    def computeStateSpace(self):

        # Putting the System in state space representation
        self.a = self.a[::-1]
        self.b = self.b[::-1]
        self.C[0] = 1
        self.D = self.a[self.n - 1] * self.b[self.n - 1]
        for i in range(self.n):
            self.Beta[i] = self.a[self.n - 1] * self.b[self.n - i - 1]
            for j in range(1, i + 1):
                self.Beta[i] -= self.a[self.n - j - 1] * self.Beta[i - j]
        self.B = self.Beta[1:]
        for i in range(self.n - 1):
            for j in range(self.n - 1):
                if (j == i + 1) and (i < self.n - 2):
                    self.A[i][j] = 1
                elif i == self.n - 2:
                    self.A[i][j] = -self.a[j]
        return self.A, self.B, self.C, self.D

    def computeStateEqautions(self):
        self.Xs[0] = np.zeros(self.n - 1)
        self.Block = np.zeros((6, self.n - 1))
        self.i = 1
        self.j = 0

        ## Calculationg Derivatives

        while (True):
            self.t = self.time[self.j]
            self.Block[0] = self.A @ self.X_temp + self.B
            self.X_mod = self.X_temp + self.h * self.Block[0]
            self.Block[1] = self.A @ self.X_mod + self.B
            self.X_mod = self.X_temp + (self.h / 2) * (self.Block[0] + self.Block[1])
            self.Block[2] = self.A @ self.X_mod + self.B
            self.X_mod = self.X_temp + (2 * self.h) * (self.Block[2])
            self.Block[3] = self.A @ self.X_mod + self.B
            self.X_mod = self.X_temp + (self.h / 12) * (5 * self.Block[0] + 8 * self.Block[2] - self.Block[3])
            self.Block[4] = self.A @ self.X_mod + self.B
            self.X_mod = self.X_temp + (self.h / 3) * (self.Block[0] + self.Block[3] + 4 * self.Block[4])
            self.Block[5] = self.A @ self.X_mod + self.B
            if self.i == self.k + 1:
                break
            self.Xs[self.i] = self.X_temp + (self.h / 12) * (5 * self.Block[0] + 8 * self.Block[2] - self.Block[3])
            self.i += 1
            self.Xs[self.i] = self.X_temp + (self.h / 3) * (self.Block[0] + 4 * self.Block[4] + self.Block[5])
            self.X_temp = self.Xs[self.i]
            self.i += 1
            self.j += 2
        if (self.inp_type == 'impulse'):
            self.u = np.zeros(self.k + 1)
            self.u[0] = 1
        return self.time, self.Xs, self.u

    def computeOutput(self):
        if self.inp_type == 'impulse':
            for j in range(self.n - 1):
                for i in range(self.k):
                    self.Xs[i, j] = (self.Xs[i + 1, j] - self.Xs[i, j]) / (self.h)
            self.Xs[-1, :] = self.Xs[-2, :]
        self.y = self.C @ self.Xs.T + self.u * self.D
        self.y[0] = self.y[1]
        return self.y

    def plotOutput(self):
        output_file("ouput.html")
        if self.inp_type == 'impulse':
            self.time = np.linspace(-self.end_t, self.end_t, self.k + 1)
            self.u = (100 * (22 / 7)) * np.exp(-(self.time * 100) ** 2)
        else:
            self.u1 = self.u.copy()
            self.u1 = np.insert(self.u1, 0, 0.5, axis=0)
            self.u1 = np.insert(self.u1, 0, 0, axis=0)
            self.u1 = np.insert(self.u1, 0, 0, axis=0)
            self.time1 = self.time.copy()
            self.time1 = np.insert(self.time1, 0, 0, axis=0)
            self.time1 = np.insert(self.time1, 0, 0, axis=0)
            self.time1 = np.insert(self.time1, 0, -100000, axis=0)

        ## Plot Input

        fig1 = figure(title='Input',
                      plot_height=400, plot_width=600,
                      x_axis_label='Time', y_axis_label='Input',
                      x_minor_ticks=2, x_range=(-5, 5), y_range=(-1, 2),
                      toolbar_location=None)
        if self.inp_type == 'impulse':
            fig1.line(x=self.time.tolist(), y=self.u.tolist(),
                      color='red', line_width=1,
                      legend_label='Input')

        else:
            fig1.line(x=self.time1.tolist(), y=self.u1.tolist(),
                      color='red', line_width=1,
                      legend_label='Input')
        tab1 = Panel(child=fig1, title="Input")

        self.time = np.linspace(0, self.end_t, self.k + 1)
        ## Plot Output

        fig2 = figure(title='System Response',
                      plot_height=400, plot_width=600,
                      x_axis_label='Time', y_axis_label='Output',
                      x_minor_ticks=2, x_range=(0, self.t),
                      y_range=(np.min(self.y) - np.mean(self.y), np.max(self.y) + np.mean(self.y)),
                      toolbar_location=None)
        fig2.line(x=self.time.tolist(), y=self.y.tolist(),
                  color='blue', line_width=1,
                  legend_label='Output')
        tab2 = Panel(child=fig2, title="Output")
        tabs = Tabs(tabs=[tab1, tab2])
        show(tabs)

    def plotStates(self):
        output_file("States.html")
        tabs = []
        for i in range(1, self.n):
            fig = figure(title='States',
                         plot_height=400, plot_width=600,
                         x_axis_label='Time', y_axis_label='State' + str(i),
                         x_minor_ticks=2, x_range=(0, self.t), y_range=(
                np.min(self.Xs[:, i - 1]) - np.abs(np.mean(self.Xs[:, i - 1])),
                np.max(self.Xs[:, i - 1]) + np.abs(np.mean(self.Xs[:, i - 1]))),
                         toolbar_location=None)
            fig.line(x=self.time.tolist(), y=self.Xs[:, i - 1].tolist(),
                     color='blue', line_width=1,
                     legend_label='State' + str(i))

            tabs.append(Panel(child=fig, title="State" + str(i)))
        tabs = Tabs(tabs=tabs)
        show(tabs)


from tkinter import *

root = Tk()
root.title("System simulator")
root.geometry("820x200")

inp_var1 = StringVar()
inp_var2 = StringVar()


def click():
    n = int(n_e.get())
    m = int(m_e.get())
    t = float(t_e.get())
    a = str(a_e.get()).split(',')
    b = str(b_e.get()).split(',')
    inp_type1 = inp_var1.get()
    inp_type2 = inp_var2.get()
    for i in range(len(a)):
        a[i] = float(a[i])
    for i in range(len(b)):
        b[i] = float(b[i])

    if n != m:
        temp = []
        for i in range(n - m):
            temp.append(0)
        b = temp + b

    def hide_frame():
        for widget in frame_3.winfo_children():
            widget.destroy()

    if inp_type1 == "unit step" and inp_type2 == "impulse":
        n2 = n;
        m2 = m;
        t2 = t;
        a2 = a;
        b2 = b
        sys1 = System(n, m, a, b, inp_type1, t)
        system2 = System(n2, m2, a2, b2, inp_type2, t2)
        A, B, C, D = sys1.computeStateSpace()
        t, Xs, inp = sys1.computeStateEqautions()
        y = sys1.computeOutput()
        A2, B2, C2, D2 = system2.computeStateSpace()
        t2, Xs2, inp2 = system2.computeStateEqautions()
        y2 = system2.computeOutput()
        hide_frame()
        button_1 = Button(frame_3, text="System Response for Unit Step   ", command=sys1.plotOutput, padx=15)
        button_1.grid(row=0, column=0)
        button_2 = Button(frame_3, text="Plot States for Unit Step       ", command=sys1.plotStates, padx=15)
        button_2.grid(row=0, column=1)
        button_3 = Button(frame_3, text="System Response for Unit Impulse", command=system2.plotOutput, padx=15)
        button_3.grid(row=0, column=2)
        button_4 = Button(frame_3, text="Plot States for Unit Impulse    ", command=system2.plotStates, padx=15)
        button_4.grid(row=0, column=3)
    elif inp_type1 == "unit step" and inp_type2 == "0":
        sys1 = System(n, m, a, b, inp_type1, t)
        A, B, C, D = sys1.computeStateSpace()
        t, Xs, inp = sys1.computeStateEqautions()
        y = sys1.computeOutput()
        hide_frame()
        button_1 = Button(frame_3, text="System Response", command=sys1.plotOutput, padx=100)
        button_1.grid(row=0, column=0)
        button_2 = Button(frame_3, text="Plot States    ", command=sys1.plotStates, padx=100)
        button_2.grid(row=0, column=1)

    elif inp_type1 == "0" and inp_type2 == "impulse":
        sys1 = System(n, m, a, b, inp_type2, t)
        A, B, C, D = sys1.computeStateSpace()
        t, Xs, inp = sys1.computeStateEqautions()
        y = sys1.computeOutput()
        hide_frame()
        button_1 = Button(frame_3, text="System Response", command=sys1.plotOutput, padx=100)
        button_1.grid(row=0, column=0)
        button_2 = Button(frame_3, text="Plot States    ", command=sys1.plotStates, padx=100)
        button_2.grid(row=0, column=1)


# root

frame_1 = LabelFrame(root, text="Input parameters", padx=5, pady=5)
frame_1.grid(row=0, column=0)

frame_2 = LabelFrame(root, text="Response type", padx=60, pady=20)
frame_2.grid(row=0, column=1)

frame_3 = LabelFrame(root, text="Visualizing", padx=5, pady=5)
frame_3.grid(row=1, column=0, columnspan=2)

# frame 1

lbl_1 = Label(frame_1, text="Output's order :           ")
lbl_1.grid(row=0, column=0)
n_e = Entry(frame_1, width=50)
n_e.grid(row=0, column=1)

lbl_2 = Label(frame_1, text="Inputs's order :             ")
lbl_2.grid(row=1, column=0)
m_e = Entry(frame_1, width=50)
m_e.grid(row=1, column=1)

lbl_3 = Label(frame_1, text="Output's coefficients : ")
lbl_3.grid(row=2, column=0)
a_e = Entry(frame_1, width=50)
a_e.grid(row=2, column=1)

lbl_4 = Label(frame_1, text="Inputs's coefficients : ")
lbl_4.grid(row=3, column=0)
b_e = Entry(frame_1, width=50)
b_e.grid(row=3, column=1)

lbl_5 = Label(frame_1, text="Maximum time :        ")
lbl_5.grid(row=4, column=0)
t_e = Entry(frame_1, width=50)
t_e.grid(row=4, column=1)

# frame 2

check_1 = Checkbutton(frame_2, text="Unit step      ", variable=inp_var1, onvalue="unit step", offvalue="0")
check_1.grid(row=0, column=0)
check_1.deselect()

check_2 = Checkbutton(frame_2, text="Unit impulse", variable=inp_var2, onvalue="impulse", offvalue="0")
check_2.grid(row=1, column=0)
check_2.deselect()

button_3 = Button(frame_2, text="Enter Parameters", command=click)
button_3.grid(row=2, column=1, rowspan=2)

mainloop()
