import numpy as np


def generate_fsd_matrix(FI):
    # function to generate the FSD matrix (generates a fragemntation matrix based on the selected fragmentation style determined by FI)
    # Initialize a 5x5 matrix with zeros
    matrix = np.zeros((5, 5))
    c1 = 0.2
    c2 = 0.15
    c3 = 0.1

    matrix[1, 0] = 1
    matrix[2, 0] = 1 - FI
    matrix[2, 1] = FI
    if FI <= 0.5:
        matrix[3, 1] = FI * 2 * c1
        matrix[4, 1] = FI * 2 * c2
        matrix[4, 2] = FI * 2 * c3
    else:
        matrix[3, 1] = c1 - ((FI - 0.5) * 2 * c1)
        matrix[4, 1] = c2 - ((FI - 0.5) * 2 * c2)
        matrix[4, 2] = c3 - ((FI - 0.5) * 2 * c3)
    matrix[3, 0] = matrix[2, 0] + (0.5 * matrix[3, 1])
    matrix[3, 2] = 1 - matrix[3, 0] - matrix[3, 1]
    matrix[4, 0] = matrix[3, 0] + (0.5 * matrix[4, 1]) + (0.25 * matrix[4, 2])
    matrix[4, 3] = 1 - matrix[4, 0] - matrix[4, 1] - matrix[4, 2]

    return matrix


# function to convert mass to number
def mass_to_num(mass_g, volume_m3, density_kg_m3):
    number = mass_g / 1000 / density_kg_m3 / volume_m3
    # number of particles has always to be integer?
    return number


# function to convert number to mass
def num_to_mass(number, volume_m3, density_kg_m3):
    mass_g = number * volume_m3 * density_kg_m3 * 1000
    return mass_g
