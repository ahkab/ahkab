import matplotlib.pyplot as plt
import numpy as np

# size
axis_fontsize = 17
legend_size = 15
size_ticks = 13
size_labels = 18 

output_voltage_style = dict(c = 'darkolivegreen', lw=2, label=r'$V_{O}$')
input_voltage_style = dict(c = 'limegreen', lw=2, label=r'$V_{in}$')
condct_time_style = dict(c = 'navy', lw=2, label=r'$g(t)$')
current_res_style = dict(c = 'darkolivegreen', lw=2, label=r'$I_r$')
current_written_style = dict(c = 'cadetblue', lw=2, label=r'$I_w$')

def plot_voltages(ax, result, plot_type):

    if plot_type=="output_voltage":

        ax.plot(result['tran']['T'], result['tran']['VN2'], **output_voltage_style)

        # print(result['tran']['T'], result['tran']['VN2'])

    if plot_type=="input_voltage":

        ax.plot(result['tran']['T'], result['tran']['VN1'], **input_voltage_style)


    ax.grid(ls=':')
    # ax.set_xlabel(r't[ms]', fontsize = axis_fontsize)

    ax.tick_params('y', labelsize=size_ticks)
    ax.tick_params('x', labelsize=size_ticks)
    ax.legend()


def plot_conductance(ax, result, plot_type):
    

    if plot_type=="analytical_sol":

        data = np.loadtxt(f"../data/solution_euler_g.txt", unpack=True)
        time = data[0]  
        conductance = data[1]

        ax.plot(time, conductance, **condct_time_style)

    if plot_type=="time_dep":

        conductance = result['tran']['I(V1)']/result['tran']['VN1']
        print(1/conductance)
        # print("-----------")
        # print(conductance[0], conductance[1])
        # print("-----------")
        
        ax.plot(result['tran']['T'], conductance, **condct_time_style)
        # ax.plot(x_interval, conductance, **input_voltage_style)
        # ax.plot(x_interval, y, **condct_time_style)


    if plot_type=="volt_dep":

        ax.plot(result['tran']['VN1'], conductance, **input_voltage_style)

    ax.grid(ls=':')
    # ax.set_ylim([0.01, 0.03])
    # ax.set_xlabel(r't[ms]', fontsize = axis_fontsize)

    ax.tick_params('y', labelsize=size_ticks)
    ax.tick_params('x', labelsize=size_ticks)
    ax.legend()


def plot_current(ax, result):
    
    
    ax.plot(result['tran']['T'], result['tran']['I(V1)'], **current_res_style)

    data = np.loadtxt(f"../data/current_memr.txt", unpack=True)
    current = data[1]
    time=data[0]

    ax.plot(time, current, **current_written_style)

    ax.grid(ls=':')

    ax.tick_params('y', labelsize=size_ticks)
    ax.tick_params('x', labelsize=size_ticks)
    ax.legend()
