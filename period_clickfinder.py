import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.timeseries import LombScargle
from PyAstronomy.pyTiming import pyPDM
import scipy.optimize

class PeriodClickfinder:

    max_freq = 1
    min_freq = 0
    number_of_freq = 10000
    frq = []
    frequency_list = []
    PS = []
    found_freq = []
    result = []
    errors = []
    resids = []
    nharm = 0
    removed = []
    width = 35

    def __init__(self, root):
        self.root = root
        self.root.title("Period Clickfinder")
        self.root.geometry()

        # Create left frame for buttons
        left_frame = ttk.Frame(root, padding="10")
        left_frame.grid(row=0, column=0, sticky="nsew")

        left_miniframe = ttk.Frame(left_frame, padding = 10)
        left_miniframe.grid(row = 6, column = 0, sticky = "nsew")

        # Create buttons
        load_data_button = tk.Button(left_frame, text="Load data", command=self.load_data)
        periodogram_button = tk.Button(left_frame, text="Calculate periodogram",
                                       command=self.periodogram)
        show_frequency_list_button = tk.Button(left_frame, text="Show frequency list",
                                               command=self.show_frequency_list)
        fit_button = tk.Button(left_frame, text="Fit", command=self.fit)
        restore_original_button = tk.Button(left_frame, text="Restore original",
                                            command=self.restore_original)
        save_button = tk.Button(left_frame, text='Save', command=self.save_output)

        wide_button = tk.Button(left_miniframe, text = 'Wide', command = lambda *args:
                      self.set_width_of_click(35))
        middle_button = tk.Button(left_miniframe, text = 'Middle', command = lambda *args:
                      self.set_width_of_click(10))
        narrow_button = tk.Button(left_miniframe, text = 'Narrow', command = lambda *args:
                      self.set_width_of_click(1))

        # Arrange buttons vertically
        load_data_button.grid(row=0, column=0, sticky="ew", pady=5)
        periodogram_button.grid(row=1, column=0, sticky="ew", pady=5)
        show_frequency_list_button.grid(row=2, column=0, sticky="ew", pady=5)
        fit_button.grid(row=3, column=0, sticky="ew", pady=5)
        restore_original_button.grid(row=4, column=0, sticky="ew", pady=5)
        save_button.grid(row=5, column=0, sticky="ew", pady=5)
        wide_button.grid(row = 0, column = 0)
        middle_button.grid(row = 0, column = 1)
        narrow_button.grid(row = 0, column = 2)

        # Create right frame for Matplotlib plots
        center_frame = ttk.Frame(root, padding="10")
        center_frame.grid(row=0, column=1, sticky="nsew")

        right_frame = ttk.Frame(root, padding="10")
        right_frame.grid(row=0, column=2, sticky="nsew")

        # Create Matplotlib figures and canvas
        self.fig1, self.ax1 = plt.subplots(figsize=(6, 4), tight_layout=True)
        self.canvas1 = FigureCanvasTkAgg(self.fig1, master=center_frame)
        self.canvas1.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1, pady=10)
        self.toolbar1 = NavigationToolbar2Tk(self.canvas1, center_frame)
        self.toolbar1.pack(side=tk.TOP, fill=tk.X)

        self.fig2, self.ax2 = plt.subplots(figsize=(6, 4), tight_layout=True)
        self.canvas2 = FigureCanvasTkAgg(self.fig2, master=center_frame)
        self.canvas2.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand=1, pady=10)
        self.toolbar2 = NavigationToolbar2Tk(self.canvas2, center_frame)
        self.toolbar2.pack(side=tk.BOTTOM, fill=tk.X)

        self.fig3, self.ax3 = plt.subplots(figsize=(4, 4), tight_layout=True)
        self.canvas3 = FigureCanvasTkAgg(self.fig3, master=right_frame)
        self.canvas3.get_tk_widget().pack(side=tk.TOP, fill=tk.X, expand=1)
        self.toolbar3 = NavigationToolbar2Tk(self.canvas3, right_frame)
        self.toolbar3.pack(side=tk.TOP, fill=tk.X)

        self.output_text = tk.Text(right_frame, height=22, width=40)
        self.output_text.pack(side=tk.BOTTOM, fill=tk.BOTH, expand=1, pady=10)

        self.canvas1.mpl_connect('button_press_event', self.on_lsc_click)
        self.canvas2.mpl_connect('button_press_event', self.on_pdm_click)

        self.nharm = 0

    def save_output(self):
        with open(self.file_name+'_fit', 'w') as f:
            f.write("# A0 \n")
            f.write(f"{self.result[0]} {self.errors[0]}\n")
            f.write("# Frequency \t Freq_error \t Amplitude \t Amp_error \t T0 \t T0_error \n")
            for i in range(1, self.nharm * 3, 3):
                f.write(f"{self.result[i]}\t{self.errors[i]}\t{self.result[i + 1]}\t"
                        f"{self.errors[i + 1]}\t{self.result[i + 2]}\t{self.errors[i + 2]}\n")

        if self.t1.any():
            PeriodClickfinder.write_fs(self.file_name + '_lsc_first', self.frq, self.PS_first)
            PeriodClickfinder.write_fs(self.file_name + '_pdm_first', self.freq, self.t1_first)

        if self.t1.any():
            PeriodClickfinder.write_fs(self.file_name + '_lsc_last', self.frq, self.PS)
            PeriodClickfinder.write_fs(self.file_name + '_pdm_last', self.freq, self.t1)

    def set_width_of_click(self, value):
        self.width = value

    @staticmethod
    def write_fs(filename, data1, data2):
        with open(filename, 'w') as f:
            f.write("# Frequency \t Amplitude \n")
            for i in range(len(data1)):
                f.write("{}\t{}\n".format(data1[i], data2[i]))

    def periodogram(self):
        self.calculate_lsc()     
        self.calculate_pdm()

    def load_data(self):
        # DO ZROBIENIA: teraz jak laduje dane to nie usuwa obrazkow, trzeba wyzerowac
        self.removed = []
        file_path = filedialog.askopenfilename(title="Select a file",
                                               filetypes=[("Data files", "*_data"),
                                                          ("All files", "*.*")])
        self.file_name = os.path.basename(file_path)
        if file_path:
            self.data = np.loadtxt(file_path)
            self.data_orig = self.data.copy()
            self.file_name = os.path.basename(file_path)
            # self.original_data = np.loadtxt(file_path)
            self.root.title("Period Clickfinder: "+self.file_name)
            self.write_output(f"Data file {self.file_name} loaded")
            self.plot_data()
            self.set_freq_limits()

    def set_freq_limits(self):
        self.max_freq = self.nyquist()
        self.min_freq = 1./(0.5*(np.max(self.data[:, 0])-np.min(self.data[:, 0])))
        self.frq = np.linspace(self.min_freq, self.max_freq, self.number_of_freq)

    def restore_original(self):
        self.data = self.data_orig.copy()
        self.nharm = 0
        self.frequency_list = []

    def nyquist(self):
        """
        fny = nyquist(t)
        Estimate effective Nyquist frequency for uneven sampling
        """
        d = np.sort(self.data[:, 0])
        return 0.5/np.median(d[1:]-d[:-1])

    def write_output(self, message):
        self.output_text.insert(tk.END, message + "\n")
        self.output_text.see(tk.END)
    
    def plot_data(self):
        self.ax3.clear()
        self.ax3.invert_yaxis()
        self.ax3.scatter(self.data_orig[:, 0], self.data_orig[:, 1])
        self.ax3.set_title("Data")
        self.canvas3.draw()

    def plot_phased(self):
        epoch = 0
        phases = ((np.array(self.data[:, 0]) - epoch) / (1. / self.found_freq[0])) % 1.0
        self.ax3.clear()
        self.ax3.set_title("Phased")
        self.ax3.invert_yaxis()
        self.ax3.scatter(phases, self.data[:, 1], color='blue')
        self.ax3.scatter(phases + 1, self.data[:, 1], color='blue')
        self.canvas3.draw()

    def sum_sines(self, time, *params):
        # print('NEW RUN')
        if isinstance(params[0], np.ndarray):
            params = params[0]
        a0 = params[0]
        y = a0
        for i in range(self.nharm):
            y += PeriodClickfinder.single_sine(time, params[1+i*3:4+i*3])
        return y

    @staticmethod
    def single_sine(time, params):
        f = params[0]
        a = params[1]
        phi = params[2]

        y = a * np.cos(2. * np.pi * f * (time - phi))
        return y
    
    def calculate_lsc(self):
        self.PS = LombScargle(self.data[:, 0], self.data[:, 1], self.data[:, 2], fit_mean=True)\
            .power(self.frq)
        if (self.data == self.data_orig).all():
            self.PS_first = self.PS
        self.write_output(f"Found periodicity (LS): {1. / self.frq[np.argmax(self.PS)]}")
        self.plot_lsc()

    def plot_lsc(self):
        self.ax1.clear()

        # Calculate moving average
        window_size = 101
        moving_average = np.convolve(self.PS, np.ones(window_size) / window_size, mode='valid')

        # Expand moving average array to match the length of self.PS
        moving_average = np.concatenate((np.full(window_size - 1, np.nan), moving_average))

        self.ax1.plot(self.frq, self.PS)
        self.ax1.plot(np.roll(self.frq, int((window_size - 1) / 2)), 4. * moving_average,
                      label='Moving Average', linestyle='--', color='red')
        self.ax1.set_title("Lomb-Scargle Periodogram")
        self.ax1.set_xlim(self.ax1.get_xlim()[0], self.ax1.get_xlim()[1])
        self.ax1.set_ylim(self.ax1.get_ylim()[0], self.ax1.get_ylim()[1])
        self.canvas1.draw()

    def on_lsc_click(self, event):
        if not self.toolbar1.mode and not self.toolbar2.mode:
            delta = self.width / (np.max(self.data[:, 0])-np.min(self.data[:, 0]))
            self.write_output(f'Clicked at x={event.xdata}, y={event.ydata}')

            freq_near_click = (self.frq >= (event.xdata - delta)) & \
                              (self.frq <= (event.xdata + delta))
            valid_indices = np.where(freq_near_click)[0]
            if len(valid_indices) > 0:
                highest_point_flattened_index = np.argmax(self.PS[freq_near_click])
                highest_point_index = valid_indices[highest_point_flattened_index]

                self.write_output(f"Frequency: {self.frq[highest_point_index]}")
                self.write_output(f"Power Spectrum Value: {self.PS[highest_point_index]}")
            else:
                self.write_output("No local maximum found in the specified frequency range.")

            if hasattr(self, 'scatter_point'):
                self.scatter_point.remove()
                self.marked_point.remove()

            self.scatter_point = self.ax1.scatter(self.frq[highest_point_index],
                                                  self.PS[highest_point_index], c='magenta',
                                                  zorder=-1)
            self.canvas1.draw()
            self.marked_point = self.ax2.axvline(self.frq[highest_point_index], c='lime',
                                                 linestyle='dotted')
            self.canvas2.draw()
            self.found_freq = [self.frq[highest_point_index], self.PS[highest_point_index], 0]
            self.plot_phased()
            self.show_prompt()

    def show_prompt(self):
        result = messagebox.askyesno("Prompt", "Do you want to include this frequency proceed?")
        if result:
            self.write_output("Frequency added to the list")
            self.frequency_list.append(self.found_freq)
            self.nharm = self.nharm + 1
            self.write_output(f"Number of signals for fitting: {self.nharm}")
        else:
            self.write_output("Frequency discarded")

    def show_frequency_list(self):
        self.frequency_window = tk.Toplevel(root)
        self.frequency_window.title("Frequency list")

        # Create a button in the new window
        button_in_frequency_window = tk.Button(self.frequency_window, text="Clear frequency list",
                                               command=self.clear_frequency_list)
        button_in_frequency_window.pack()

        if self.frequency_list:
            for which_freq, entry in enumerate(self.frequency_list):
                entry_frame = tk.Frame(self.frequency_window)
                entry_frame.pack(padx=20, pady=10, anchor='w')

                delete_button = tk.Button(entry_frame, text='Delete',
                                          command=lambda number=which_freq:
                                          self.delete_frequency(number))
                delete_button.pack(side='left')
                text_label = tk.Label(entry_frame, text=str(entry))
                text_label.pack(side='right')

        else:
            text_label = tk.Label(self.frequency_window, text="Frequency list is empty.")
            text_label.pack(padx=20, pady=20)

    def delete_frequency(self, number):
        self.write_output(f'Removing frequency: {self.frequency_list[number]}')
        self.frequency_list.pop(number)
        self.nharm = self.nharm-1
        self.frequency_window.destroy()
        self.show_frequency_list()
        self.fit()

    def clear_frequency_list(self):
        self.frequency_list = []
        self.nharm = 0
        self.frequency_window.destroy()
        self.show_frequency_list()
        self.restore_original()

    def fit(self):
        guess_amp = 0.5 * (np.max(self.data[:, 1]) - np.min(self.data[:, 1]))
        guess_offset = np.mean(self.data[:, 1])
        guess_freq = self.frequency_list[-1][0]
        guess_phase = self.data[-1, 0] - (0.75 / self.frequency_list[0][0])

        if self.nharm == 1:
            guess = np.array([guess_offset, guess_freq, guess_amp, guess_phase])
        elif self.nharm > 1:
            guess = np.concatenate((self.result, [guess_freq, guess_amp, guess_phase]))

        lower_bounds = [-np.inf] + [0, 0, self.data[-1, 0] -
                                    (1.5 / self.frequency_list[0][0])] * self.nharm
        upper_bounds = [np.inf] + [np.inf, np.inf, self.data[-1, 0]] * self.nharm

        print(guess)
        print(lower_bounds)
        print(upper_bounds)
        popt, pcov = scipy.optimize.curve_fit(self.sum_sines, self.data_orig[:, 0],
                                              self.data_orig[:, 1], p0=guess,
                                              bounds=(tuple(lower_bounds), tuple(upper_bounds)))

        self.result = popt
        self.errors = np.sqrt(np.diag(pcov))
        print("RESULT: ", self.result)
        print("ERRORS: ", self.errors)

        print('Getting residuals')
        self.resids = self.data_orig[:, 1] - self.sum_sines(self.data_orig[:, 0], self.result)
        self.data[:, 1] = self.resids
        # self.data[:, 1] = np.sqrt(self.resids ** 2 / self.data_orig[:, 2] ** 2)


        self.update_frequency_list()
        self.calculate_lsc()
        self.calculate_pdm()
        self.plot_fit()

    def plot_fit(self):
        epoch = 0
        fit = self.sum_sines(self.data_orig[:, 0], self.result)
        phases = ((np.array(self.data_orig[:, 0]) - epoch) / (1./self.result[1])) % 1.0
        sorted_idx = np.argsort(phases)
        self.ax3.clear()
        self.ax3.invert_yaxis()
        self.ax3.set_title("Phased")
        self.ax3.scatter(phases, self.data_orig[:, 1], color='blue')
        self.ax3.scatter(phases+1, self.data_orig[:, 1], color='blue')
        self.ax3.plot(phases[sorted_idx], fit[sorted_idx], c='red')
        self.ax3.plot(phases[sorted_idx] + 1, fit[sorted_idx], c='red')
        self.canvas3.draw()

    def on_pdm_click(self, event):
        if not self.toolbar1.mode and not self.toolbar2.mode:
            delta = self.width / (np.max(self.data[:, 0]) - np.min(self.data[:, 0]))
            self.write_output(f'Clicked at x={event.xdata}, y={event.ydata}')

            freq_near_click = (self.freq >=
                               (event.xdata - delta)) & \
                              (self.freq <= (event.xdata + delta))
            valid_indices = np.where(freq_near_click)[0]
            if len(valid_indices) > 0:
                lowest_point_flattened_index = np.argmin(self.t1[freq_near_click])
                lowest_point_index = valid_indices[lowest_point_flattened_index]

                self.write_output(f"Frequency: {self.freq[lowest_point_index]}")
                self.write_output(f"Power Spectrum Value: {self.t1[lowest_point_index]}")
            else:
                self.write_output("No local maximum found in the specified frequency range.")

            if hasattr(self, 'scatter_point'):
                self.scatter_point.remove()
                self.marked_point.remove()

            self.scatter_point = self.ax2.scatter(self.freq[lowest_point_index],
                                                  self.t1[lowest_point_index], c='red')
            self.canvas2.draw()

            self.marked_point = self.ax1.axvline(self.freq[lowest_point_index], c='lime',
                                                 linestyle='dotted')
            self.canvas1.draw()

            self.found_freq = [self.freq[lowest_point_index], self.t1[lowest_point_index], 0]
            self.plot_phased()
            self.show_prompt()

    def update_frequency_list(self):
        for i in range(self.nharm):
            self.frequency_list[i][0] = self.result[1 + 3 * i]
            self.frequency_list[i][1] = self.result[2 + 3 * i]
            self.frequency_list[i][2] = self.result[3 + 3 * i]

    def calculate_pdm(self):
        S = pyPDM.Scanner(minVal=self.min_freq, maxVal=self.max_freq, dVal=0.0001,
                          mode="frequency")
        P = pyPDM.PyPDM(self.data[:, 0], self.data[:, 1])
        self.freq, self.t1 = P.pdmEquiBinCover(10, 3, S)

        if (self.data == self.data_orig).all():
            self.t1_first = self.t1

        self.plot_pdm()

    def plot_pdm(self):
        self.ax2.clear()
        self.ax2.plot(self.freq, self.t1)
        self.ax2.set_title("Period Dispersion Minimization")
        self.ax2.set_xlim(self.ax2.get_xlim()[0], self.ax2.get_xlim()[1])
        self.ax2.set_ylim(self.ax2.get_ylim()[0], self.ax2.get_ylim()[1])
        self.canvas2.draw()


if __name__ == "__main__":
    root = tk.Tk()
    app = PeriodClickfinder(root)
    root.geometry("1300x950")
    root.mainloop()
