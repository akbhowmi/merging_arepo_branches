""" Idealized test, AMR Kelvin-Helmholtz (static) in 2D. """
import numpy as np
import sys

sys.path.append('./')  # test.py
sys.path.append('./examples/')  # test.sh
import utils

Lx = 1.0
Ly = 1.0

#numTasksMPI = 2 # ngbtree error with any


def create_ics(path, filename='ics.hdf5'):
    # config
    N_x = 64
    N_y = 64

    P = 2.5
    omega0 = 0.1
    sigma = 0.05 / np.sqrt(2)
    GAMMA = 5.0 / 3.0

    # derived
    dx = Lx / N_x
    dy = Ly / N_y

    N = int(N_x * N_y / 2 + N_x * N_y * 4 / 2)
    ind4 = int(N_x * N_y / 4)
    ind2 = int(N_x * N_y / 2)
    ind34 = int(N_x * N_y * 3 / 4)

    # allocate
    pos = np.zeros((N, 3), dtype='float64')
    vel = np.zeros((N, 3), dtype='float64')
    u = np.zeros(N, dtype='float64')
    mass = np.zeros(N, dtype='float64')

    # set gas cell propeties
    mesh = np.meshgrid(np.arange(N_x), np.arange(N_y))

    posx = (mesh[0] * dx).reshape(N_x * N_y) + 0.5 * dx
    posy = (mesh[1] * dy).reshape(N_x * N_y) + 0.5 * dy

    mesh = np.meshgrid(np.arange(N_x), np.arange(N_y))
    posx2 = (mesh[0] * dx).reshape(N_x * N_y) + 0.5 * dx
    posy2 = (mesh[1] * dy).reshape(N_x * N_y) + 0.5 * dy

    pos[:ind4, 0] = posx[:ind4]
    pos[:ind4, 1] = posy[:ind4]
    mass[:ind4] = 1. * dx * dy
    u[:ind4] = P / ((GAMMA - 1) * 1.)
    vel[:ind4, 0] = -0.5

    pos[ind4:ind4 + ind2, 0] = posx2[ind4:ind34]
    pos[ind4:ind4 + ind2, 1] = posy2[ind4:ind34]
    mass[ind4:ind4 + ind2] = 2. * dx * dy
    u[ind4:ind4 + ind2] = P / ((GAMMA - 1) * 2.)
    vel[ind4:ind4 + ind2, 0] = +0.5

    pos[-ind4:, 0] = posx[-ind4:]
    pos[-ind4:, 1] = posy[-ind4:]
    mass[-ind4:] = 1. * dx * dy
    u[-ind4:] = P / ((GAMMA - 1) * 1.)
    vel[-ind4:, 0] = -0.5

    # perturbation
    vel[:, 1] = omega0 * np.sin(4 * np.pi * pos[:, 0]) * (
        np.exp(-(pos[:, 1] - 0.25)**2 * 0.5 /
               (sigma**2)) + np.exp(-(pos[:, 1] - 0.75)**2 * 0.5 / (sigma**2)))

    ids = np.arange(N) + 1

    # write
    pt0 = {
        'Coordinates': pos,
        'Velocities': vel,
        'Masses': mass,
        'InternalEnergy': u,
        'ParticleIDs': ids
    }
    h = {'Flag_DoublePrecision': 1}

    utils.write_ic_file(path + filename, {'PartType0': pt0},
                        boxSize=Lx,
                        headerAttrs=h)


if __name__ == '__main__':

    try:
        create_ics(path=sys.argv[1])
    except:
        sys.exit(1)

    sys.exit(0)  # normal exit
