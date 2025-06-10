import numpy as np
from datetime import datetime
import warnings
warnings.filterwarnings("ignore")


if __name__ == '__main__':
    #constantsFile = PATH TO TXT FILE WITH CONSTANTS FROM LEOS 
    
    cycles = []
    with open(constantsFile, 'r') as f:
        lines = f.readlines()
        for l in lines:
            if l.startswith('New cycle:'):

                date_str = l.split(': ')[-1]
                date_str = date_str.strip()
                date_object = datetime.strptime(date_str, '%Y%m%d')

                cycles.append(date_object)

    print(cycles)

    cycleSTA = cycles[0]
    cycleEND = cycles[-1]

    #saveto = PATH TO TXT FILE WHERE YOU WANT TO SAVE THE AVERAGED CONSTANTS
    alpha = []
    beta = []
    T0 = []
    S0 = []
    dens0 = []
    
    with open(constantsFile, 'r') as f:
        lines = f.readlines()
        for l in lines:
            if l.startswith('Alpha:'):
                alpha.append(float(l.split(':')[-1]))
            elif l.startswith('Beta:'):
                beta.append(float(l.split(':')[-1]))
            elif l.startswith('T0:'):
                T0.append(float(l.split(':')[-1]))
            elif l.startswith('S0:'):
                S0.append(float(l.split(':')[-1]))
            elif l.startswith('dens0:'):
                dens0.append(float(l.split(':')[-1]))

    alpham = np.mean(np.array(alpha))
    betam = np.mean(np.array(beta))
    T0m = np.mean(np.array(T0))
    S0m = np.mean(np.array(S0))
    dens0m = np.mean(np.array(dens0))

    with open(saveto, 'w') as f:
        f.write('Alphas: ' + str(alpha) + '\n')
        f.write('Betas: ' + str(beta) + '\n')
        f.write('T0s: ' + str(T0) + '\n')
        f.write('S0s: ' + str(S0) + '\n')
        f.write('dens0s: ' + str(dens0) + '\n\n')

        f.write('Final constants: \n')
        f.write('Alpha:' + str(alpham) + '\n')
        f.write('Beta:' + str(betam) + '\n')
        f.write('T0:' + str(T0m) + '\n')
        f.write('S0:' + str(S0m) + '\n')
        f.write('dens0:' + str(dens0m) + '\n')

  