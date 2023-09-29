import matplotlib.pyplot as plt
from numpy import arange

def draw_properties(data, gasa, gasb):
    gasa_str = gasa.replace('2', r'$_{2}$')
    if gasb == 'CH4':
        gasb_str = gasb.replace('4', r'$_{4}$')
    else:
        gasb_str = gasb.replace('2', r'$_{2}$')
    data['alpha_'+gasa+'_'+gasb] = data[gasa] - data[gasb]
    if gasb == 'CH4':
        data['k_'+gasa+'_'+gasb] = 0.3794 * data[gasa] + data['alpha_'+gasa+'_'+gasb]
    else:
        data['k_'+gasa+'_'+gasb] = 0.2933 * data[gasa] + data['alpha_'+gasa+'_'+gasb]
    fig, ax = plt.subplots(figsize=(4.8,3.33))   #4.12, 3.33
    plt.rc('font',family='Arial', weight='medium')
    ax.set_facecolor('white')
    kuang = ['top','bottom','left','right']
    for i in kuang:
        ax.spines[i].set_color('black')
        ax.spines[i].set_linewidth(1.5)

    im = ax.scatter(data['CO2'], data['alpha_'+gasa+'_'+gasb],s=0.3, c='b', alpha =1)
    ax.set_xlabel('log$_{10}$(P)', fontsize=12)
    ax.set_ylabel(r'log$_{10}$($\alpha$)', fontsize=12)
    ax.set_xticks(arange(0,5.01))
    ax.set_yticks(arange(0,5.01))
    ax.set_xlim(-1,5)
    ax.set_ylim(-1,5)
    plt.tick_params(labelsize=12)

    ax.text(-0.8,4.5, gasa_str+'/'+gasb_str, fontsize=12)

    ax.axline((0,2.5531), slope=-0.3794,ls="--", c=".3", label='2008 upper bound')
    ax.axline((0,2.7898), slope=-0.3794,ls=":", c=".3",  label='2019 upper bound')
    leg = ax.legend(fontsize=12, framealpha=0.1)
    leg.get_frame().set_linewidth(0.0)

    return plt