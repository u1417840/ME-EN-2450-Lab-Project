import numpy as np

def gaussian_plume_dep(X, Y, WindSpeed, WindDir, dep_area, Q):
    # Parameters for simple spread rate model (following Miller et al., 2018 perpendicular spread fits)
    sigy0 = 1.0
    my = 1.0
    set_vel = 0.1  # Accounts for gravity and deposition to leaves

    # Calculate the 'outgoing' wind direction in radians
    WindDir = np.radians(WindDir)
    dep_length = np.sqrt(dep_area / 10000)  # Convert to m^2

    # Initialize deposition array
    C = np.full(len(X), np.nan)
    X = np.array(X)
    Y = np.array(Y)
    # Calculate distance in rotated coords
    x = X * np.cos(WindDir) + Y * np.sin(WindDir)
    PosX = x > 0  # Skip if it's invalid (Gauss Plumes are invalid in the -dir)

    y = -X * np.sin(WindDir) + Y * np.cos(WindDir)
    sigy = np.sqrt(sigy0**2 + my**2 * x[PosX]**2)

    C[PosX] = (1. / (2 * np.pi * sigy * WindSpeed)) * np.exp(-((y[PosX])**2) / (2 * sigy**2))
    C[PosX] = np.exp(-0.05 * np.abs(x[PosX])) * C[PosX]  # Miller et al (2018) removal rate
    C[PosX] = np.real(set_vel * C[PosX] * dep_length * Q)

    # Note that this C is the deposition in amount/second (a rate)
    return C