
# coding: utf-8

# Calculating achievable depth and quantum volume as a function of experimental paramaters for a trapped ion design

import numpy as np
import matplotlib.pyplot as plt


# Shuttling dependence from ion routing model. See paper for more details
def shuttlesFromQubit(qubit):
    """ Number of shuttling operations (defined as shuttling between two
        adjacent X-Junctions) as a function of qubit number, N
    Parameters
    ----------
    totalError: float

    qubit: int
        Number of qubits, N

    Returns
    ----------
    total shuttling time scaling: float
        tau = 1.3 * sqrt(N) + 2
    """
    return 1.3 * np.sqrt(qubit) + 2

# Number of times an ion on average passes through the centre of an X-Junction, as a function of qubit number


def xPassFromQubit(qubit):
    """ Number of times an ion on average passes through the centre of an
        X-Junction, as a function of qubit number, N
    Parameters
    ----------
    totalError: float

    qubit: int
        Number of qubits, N

    Returns
    ----------
    X_count (paper): float
        X_count = [0.6/sqrt(2) * sqrt(N)] + 2.4
    """
    return (0.6/np.sqrt(2)) * np.sqrt(qubit) + 2.4

# Depth overhead for swapping on a superconducting square grid using CQC's publically available tket


def scDepthOverhead(qubit):
    return max(1, ((2.77*(qubit**0.5)-4.53)))


# Equations to calculate QV_native for the various architectures

def fid2Error(gatefid):
    """ Convert a gate fidelity (%) into an error
    Parameters
    ----------
    gatefid: float
        fidelity of gate
    Returns
    ----------
    gate error: float
        100 - fidelity / 100
    """
    return (100 - gatefid)/100


def depth(totalError, qubit):
    """ Calculate the achievable depth, D, as a function of the total effective
        error, epsilon_eff, and qubit number, N.
    Parameters
    ----------
    totalError: float

    qubit: int
        Number of qubits, N

    Returns
    ----------
    depth: float
        D = 1 / (N * epsilon_eff)
    """
    return 1/(qubit*(totalError))

# Calculate trapped ion total effective error


def ionEffectiveError(qubit, gatefid, shuttleTime, coherenceTime, ionLossRate):
    """ Calculate trapped ion total effective error
    Parameters
    ----------
    qubit: int
        Number of qubits, N

    gatefid: float
        fidelity of gate

    total shuttling time scaling: float
        tau = 1.3 * sqrt(N) + 2

    coherenceTim : float
        time in seconds (c in paper)

    ionLossRate : float
        rate of ions lost per shuttle rate (X_loss in paper)

    Returns
    ----------
    effective ion error: float
        Equation 2 in paper:

        epsilon_eff = epsilon_gate + (1 - e^{-t / c}) + (X_count * X_loss)

        t = time spent shuttling and on separating and merging [microseconds]
        X_count = mean number of ions passing over X-junctions
    """

    t = (shuttleTime*shuttlesFromQubit(qubit)+2*separationMerge)*10**(-6)

    # Error from decoherence 1-e^(t/c)
    errorFromShuttling = 1 - np.exp(-t / coherenceTime)

    # Ion Loss error (X_count * X_loss)
    ionLoss = xPassFromQubit(qubit) * ionLossRate

    # Assuming gate requirement is the native two qubit gate
    errorFromGates = fid2Error(gatefid)

    errorFromConnectivity = errorFromShuttling + ionLoss

    return errorFromGates + errorFromConnectivity

# QV native as a function of qubit number


def QVnative(architecture, qubit, gatefid, square, coherenceTime):
    """ Quantum Vplume (QV) native as a function of qubit number
    Parameters
    ----------
    architecture: int
        1 = Ion Trap, 2 = All-to-All Connected, 3 = Superconducting

    qubit: int
        Number of qubits, N

    gatefid: float
        fidelity of gate

    total shuttling time scaling: float
        tau = 1.3 * sqrt(N) + 2

    square: int
        square QV if square == 1, else plot log_2(QV)
    Returns
    ----------
    effective ion error: float
        Equation 2 in paper:

        epsilon_eff = epsilon_gate + (1 - e^{-t / c}) + (X_count * X_loss)

        t = time spent shuttling and on separating and merging [microseconds]
        X_count = mean number of ions passing over X-junctions
    """
    # Ion Trap Archirecture
    if architecture == 1:
        totalError = ionEffectiveError(
            qubit, gatefid, shuttleSpeed, coherenceTime, ionLossRate)
    
    # All-to-All Connected Archirecture
    elif architecture == 2:
        totalError = fid2Error(gatefid)

    # Superconducting where depth overhead scales as 2.77N^0.5 -4.53 (adjust above)
    elif architecture == 3:
        totalError = scDepthOverhead(qubit) * fid2Error(gatefid)

    circuitDepth = depth(totalError, qubit)

    if square == 1:
        QV = min(circuitDepth, qubit)**2

    else:
        QV = min(circuitDepth, qubit)

    return QV

# The peak value of QV as a function of qubit number
def maxQVnative(architecture, qubits, gatefid, square, coherenceTime):
    """ Quantum Vplume (QV) native as a function of qubit number
    Parameters
    ----------
    architecture: int
        1 = Ion Trap, 2 = All-to-All Connected, 3 = Superconducting

    qubits: list 
        qubit numbers from a defined minimum to maximum in steps of 2

    gatefid: float
        fidelity of gate
        
    square: int
        square QV if square == 1, else plot log_2(QV)

    Returns
    ----------
    maxQV: float
        Equation 2 in paper:

        epsilon_eff = epsilon_gate + (1 - e^{-t / c}) + (X_count * X_loss)

        t = time spent shuttling and on separating and merging [microseconds]
        X_count = mean number of ions passing over X-junctions
    """
    maxQV = 0
    
    for qubit in qubits:    
        newQV = QVnative(architecture, qubit, gatefid, square, coherenceTime)
        if newQV < maxQV:
            return maxQV
        maxQV = newQV
   
    return maxQV


# Qubit number range investigated
qubit_min = 2
qubit_max = 120 #25 #120
qubit_step = 2
qubits = np.arange(qubit_min, qubit_max + qubit_step, qubit_step)

# Gate fidelity range investigated
gatefidmin = 99 #10
gatefidmax = 99.99 #99.9
dataPoints = 10000
step = (gatefidmax - gatefidmin)/dataPoints
gatefids = np.arange(gatefidmin, gatefidmax + step, step)

# Convert gate fidelity % to inverse error
two_qubit_errors = 1 / ((100 - gatefids) / 100)

# Ion trapping experimental paramaters

# Coherence time (seconds)
coherenceTime = 2.13
coherenceTimeLong = coherenceTime * 10
# Likely hood of loss per X-junction travel
ionLossRate = 10**-6
# Time to shuttle between two adjacent X-junctions (microseconds)
shuttleSpeed = 114
# Time to perform a separation or merge (microseconds)
separationMerge = 80

# square QV if square==1, else plot log_2(QV)
square = 0

#Calculate and plot
ion_5 = []
ion_6 = []
sc = []
all2all = []

for gatefid in gatefids:
    ion_5.append(maxQVnative(1, qubits, gatefid, square, coherenceTime))
    ion_6.append(maxQVnative(1, qubits, gatefid, square, coherenceTimeLong))
    all2all.append(maxQVnative(2, qubits, gatefid, square, coherenceTime))
    sc.append(maxQVnative(3, qubits, gatefid, square, coherenceTime))

plt.rcParams.update({'font.size': 14})
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

plt.plot(two_qubit_errors, all2all, 'r-', label="All-to-All")
plt.plot(two_qubit_errors, ion_5, 'g-', label="Ion "+r'$t_{s}/t_{c}\ \approx\ 5 \times 10^{-5}$')
plt.plot(two_qubit_errors, ion_6, 'b-', label="Ion "+r'$t_{s}/t_{c}\ \approx\ 5 \times 10^{-6}$')
plt.plot(two_qubit_errors, sc, 'y-', label="Superconducting")
label = 'Ion'
plt.legend(loc="best")
plt.xlabel('1/'+r'$\epsilon$')
plt.ylabel('$log_2$ ($QV_{native}$)')
plt.xlim(min(two_qubit_errors), max(two_qubit_errors))
plt.grid()
ax.set_xscale('log')
plt.savefig(label+'.png', dpi = 300)
 