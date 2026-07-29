import numpy as np
import matplotlib.pyplot as plt
import math

def plot_delta_e2(output_folder):
    plt.clf() #Clear the current figure (prevents multiple labels)

    labelfont = {
            'family' : 'sans-serif',  # (cursive, fantasy, monospace, serif)
            'color'  : 'black',       # html hex or colour name
            'weight' : 'normal',      # (normal, bold, bolder, lighter)
            'size'   : 14,            # default value:12
            }

    titlefont = {
            'family' : 'sans-serif',
            'color'  : 'black',
            'weight' : 'normal',
            'size'   : 18,
            }

    pi = np.pi
    Ek = 1.33
    Ek = 0.662
    x = np.linspace(0.0, Ek, 1000)
    # f1 = np.arccos(1 - 0.511*(1.0/(1.17 - x) - 1.0/1.17)) - np.arccos(1 - 0.511*(1.0/(1.17 - x - 0.4) - 1.0/1.17))
    # f2 = np.arccos(1 - 0.511*(1.0/(1.17 - x) - 1.0/1.17)) - np.arccos(1 - 0.511*(1.0/(1.17 - x - 0.2) - 1.0/1.17))


    compton_edge = Ek * (2.0 * Ek) / (0.511 + 2.0 * Ek)
    for E1 in np.linspace(0.05, 0.95*compton_edge, 6):
    # for E1 in [0.20, 0.40, 0.6, 0.8]:
        f = np.arccos(1.0 - 0.511*(1.0/(Ek - E1 - x) - 1.0/(Ek - x))) - np.arccos(1 - 0.511*(1.0/(Ek - E1) - 1.0/Ek))
        # f = np.arccos(1.0 - 0.511*(1.0/(Ek - E1 - x) - 1.0/(Ek - x)))
        # f = np.arccos(1 - 0.511*(1.0/(Ek - E1) - 1.0/Ek)) + 0*x
        # f *= 180.0/np.pi
        f = np.sin(f)
        notnans = np.logical_not(np.isnan(f))
        f = f[notnans]
        x = x[notnans]

        i_max = np.argmax(f, axis=None, out=None)
        print "max i:", E1, Ek, compton_edge, i_max, np.max(f), len(f)

        plt.plot(x[:i_max], f[:i_max],
                 linestyle='-',                    # line style
                 linewidth=3,                       # line width
                 label='E$_1 =$ %0.2f MeV' % (E1))      # plot label

    axes = plt.gca()
    axes.set_xlim([0.0, 0.8])            # x-axis bounds
    # axes.set_ylim([0.0, 90.0])              # y-axis bounds

    legend = plt.legend(loc='lower right', shadow=True, fontsize='small', title="First Scatter Energy [MeV]")

    plt.title('Cone Angular Displacement ($\gamma=1.33$ MeV)', fontdict=titlefont)
    plt.xlabel('$\delta$E$_2$ [MeV]', fontdict=labelfont)
    plt.ylabel('sin($\delta\\theta$)', fontdict=labelfont)

    plt.subplots_adjust(left=0.15)        # prevents overlapping of the y label

    plt.savefig('%s/delta_e2.png' % output_folder, bbox_inches='tight')


def plot_E1_dependence(output_folder):
    plt.clf() #Clear the current figure (prevents multiple labels)

    labelfont = {
            'family' : 'sans-serif',  # (cursive, fantasy, monospace, serif)
            'color'  : 'black',       # html hex or colour name
            'weight' : 'normal',      # (normal, bold, bolder, lighter)
            'size'   : 14,            # default value:12
            }

    titlefont = {
            'family' : 'sans-serif',
            'color'  : 'black',
            'weight' : 'normal',
            'size'   : 18,
            }

    pi = np.pi
    Ek = 1.33
    # Ek = 4.44
    Ek = 0.662
    compton_edge = Ek * (2.0 * Ek) / (0.511 + 2.0 * Ek)
    x = np.linspace(0.001, 0.95*compton_edge, 10000)
    x = x[x < 3.0]
    # f1 = np.arccos(1 - 0.511*(1.0/(1.17 - x) - 1.0/1.17)) - np.arccos(1 - 0.511*(1.0/(1.17 - x - 0.4) - 1.0/1.17))
    # f2 = np.arccos(1 - 0.511*(1.0/(1.17 - x) - 1.0/1.17)) - np.arccos(1 - 0.511*(1.0/(1.17 - x - 0.2) - 1.0/1.17))

    for dE2 in np.linspace(0.005, 0.10*compton_edge, 6):
    # for dE2 in [0.01, 0.05, 0.10, 0.20, 0.30, 0.50]:

        f = np.arccos(1.0 - 0.511*(1.0/(Ek - x - dE2) - 1.0/(Ek - dE2))) - np.arccos(1 - 0.511*(1.0/(Ek - x) - 1.0/Ek))
        # f = np.arccos(1.0 - 0.511*(1.0/(Ek - E1 - x) - 1.0/(Ek - x)))
        # f = np.arccos(1 - 0.511*(1.0/(Ek - E1) - 1.0/Ek)) + 0*x
        # f *= 180.0/np.pi
        f = np.sin(f)
        notnans = np.logical_not(np.isnan(f))
        f = f[notnans]
        x = x[notnans]

        i_max = np.argmax(f, axis=None, out=None)
        print "max i:", dE2, i_max, np.max(f), len(f)

        plt.plot(x[:i_max], f[:i_max],
                 linestyle='-',
                 linewidth=3,
                 label='%0.2f MeV' % (dE2))

    axes = plt.gca()
    axes.set_xlim([0.0, 0.85])
    # axes.set_ylim([0.0, 90.0])

    legend = plt.legend(loc='upper left', shadow=True, fontsize='small', title="$\delta$E")

    plt.title('Cone Angular Displacement [$\gamma=$%.3f MeV]' % Ek, fontdict=titlefont)
    plt.xlabel('E$_1$ [MeV]', fontdict=labelfont)
    plt.ylabel('sin($\delta\\theta$)', fontdict=labelfont)

    plt.subplots_adjust(left=0.15)        # prevents overlapping of the y label

    plt.savefig('%s/E1_dependence.png' % output_folder, bbox_inches='tight')
    plt.savefig('%s/E1_dependence.eps' % output_folder, bbox_inches='tight')
    plt.savefig('/y_drive/projects/PromptGamma/Images/E1_dependence.pdf' , bbox_inches='tight')


def plot_wrong_label_penalty(output_folder):
    plt.clf() #Clear the current figure (prevents multiple labels)

    labelfont = {
            'family' : 'sans-serif',  # (cursive, fantasy, monospace, serif)
            'color'  : 'black',       # html hex or colour name
            'weight' : 'normal',      # (normal, bold, bolder, lighter)
            'size'   : 14,            # default value:12
            }

    titlefont = {
            'family' : 'sans-serif',
            'color'  : 'black',
            'weight' : 'bold',
            'size'   : 16,
            }

    pi = np.pi
    Ek = 1.17
    Ek = [1.17, 1.33]
    # Ek = [4.44, 6.13]
    x = np.linspace(0.0, 3.0, 1000)
    # f1 = np.arccos(1 - 0.511*(1.0/(1.17 - x) - 1.0/1.17)) - np.arccos(1 - 0.511*(1.0/(1.17 - x - 0.4) - 1.0/1.17))
    # f2 = np.arccos(1 - 0.511*(1.0/(1.17 - x) - 1.0/1.17)) - np.arccos(1 - 0.511*(1.0/(1.17 - x - 0.2) - 1.0/1.17))

    # f = np.arccos(1.0 - 0.511*(1.0/(Ek[0] - x) - 1.0/(Ek[0]))) - np.arccos(1.0 - 0.511*(1.0/(Ek[1] - x) - 1.0/(Ek[1])))

    f1 = np.arccos(1.0 - 0.511*(1.0/(Ek[1] - x) - 1.0/(Ek[1])))
    sinf1 = np.sin(f1)

    f2 = np.arccos(1.0 - 0.511*(1.0/(Ek[0] - x) - 1.0/(Ek[0])))
    sinf2 = np.sin(f2)

    f = np.arccos(1.0 - 0.511*(1.0/(Ek[0] - x) - 1.0/(Ek[0]))) - np.arccos(1.0 - 0.511*(1.0/(Ek[1] - x) - 1.0/(Ek[1])))
    sinf = np.sin(f)
    # f *= 180.0/np.pi

    #
    # plt.plot(x, f,
    #          linestyle='-',                    # line style
    #          linewidth=5,                       # line width
    #          label='$\Delta\\theta$')      # plot label


    plt.plot(x, f1*180/math.pi,
             linestyle='-',                    # line style
             linewidth=1,                       # line width
             label='$\\theta$(1.17 MeV)')      # plot label
    plt.plot(x, f2*180/math.pi,
             linestyle='-   ',                    # line style
             linewidth=1,                       # line width
             label='$\\theta$(1.33 MeV)')      # plot label
    plt.plot(x, f*180/math.pi,
             linestyle='-',                    # line style
             linewidth=3,                       # line width
             label='$\\theta$(1.33 MeV) - $\\theta$(1.17 MeV)')      # plot label
    axes = plt.gca()
    axes.set_xlim([0.0, 0.9])
    axes.set_ylim([0.0, 120.0])

    legend = plt.legend(loc='upper left', shadow=True, fontsize='small')

    plt.title(' ', fontdict=titlefont)
    plt.xlabel('E$_1$ [MeV]', fontdict=labelfont)
    plt.ylabel('$\\delta\\theta$ [degrees]', fontdict=labelfont)

    plt.subplots_adjust(left=0.15)        # prevents overlapping of the y label

    plt.savefig('%s/wronglabel.png' % output_folder, bbox_inches='tight')


#------------- MAIN -------------------------------------------------------
if __name__ == "__main__":

    # cProfile.run("main()", sort=2)
    # cProfile.run("main()", sort=2)
    OUTPUT_FOLDER = "/home/dsmackin/public_html/cc"
    # plot_wrong_label_penalty(OUTPUT_FOLDER)
    plot_E1_dependence(OUTPUT_FOLDER)
    plot_delta_e2(OUTPUT_FOLDER)