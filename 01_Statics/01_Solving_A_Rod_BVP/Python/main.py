import math
import numpy as np
import scipy.optimize
from scipy import integrate
import matplotlib.pyplot as plt


def main():
    # Independent Parameters
    E = 200e9  # Young's modulus
    G = 80e9  # Shear modulus
    r = 0.001  # Cross-sectional radius
    rho = 8000  # Density
    g = np.array([9.81, 0, 0]).T  # Gravitational acceleration
    L = 0.5  # Length(before strain)

    # Dependent Parameters
    A = math.pi * r ** 2  # Cross-sectional area
    I = math.pi * r ** 4 / 4  # Area moment of inertia
    J = 2 * I  # Polar moment of inertia

    Kse = np.diag([G * A, G * A, E * A])  # Stiffness matrices
    Kbt = np.diag([E * I, E * I, G * J])

    # Boundary Conditions
    p0 = np.array([0, 0, 0]).T
    R0 = np.eye(3)
    pL = np.array([0, -0.1*L, 0.8*L]).T
    RL = np.eye(3)

    # Sub-functions
    def obj_fun(guess):  # Optimization objective function
        # Update guessed initial conditions
        n0 = guess[0:3]
        m0 = guess[3:6]
        y0 = np.concatenate([p0, np.reshape(R0, 9), n0, m0])

        # Numerically solve the IVP
        nonlocal Y
        Y = integrate.solve_ivp(rod_ode, (0, L), y0, max_step=0.01).y

        # Calculate distal constraint violation
        pL_shot = Y[0:3, -1]
        RL_shot = np.reshape(Y[3:12, -1], (3, 3))
        position_error = pL_shot - pL
        rotation_error = inv_hat(RL_shot.T @ RL - RL_shot @ RL.T)
        return np.concatenate([position_error, rotation_error])

    def rod_ode(s, y):  # State vector derivative function
        del s  # Integration variable unused in autonomous ODE
        # Unpack state vector
        R = np.reshape(y[3:12], (3, 3))
        n = y[12:15].T
        m = y[15:18].T

        # Constitutive equation
        v = np.linalg.inv(Kse) @ R.T @ n + np.array([0, 0, 1]).T
        u = np.linalg.inv(Kbt) @ R.T @ m

        # Static Cosserat rod equations - system of nonlinear ODEs
        ps = R @ v
        Rs = R @ hat(u)
        ns = -rho * A * g
        ms = -np.cross(ps, n)

        # Pack state vector derivative
        return np.concatenate([ps, np.reshape(Rs, 9), ns, ms]).T

    def hat(y):
        return np.array([[0, -y[2], y[1]],
                         [y[2], 0, -y[0]],
                         [-y[1], y[0], 0]])

    def inv_hat(skew):
        return np.array([skew[2, 1], skew[0, 2], skew[1, 0]])

    # Numerical Integration
    init_guess = np.zeros(6)
    Y = None
    sol = scipy.optimize.root(obj_fun, init_guess, method="lm").x

    # Visualization
    ax = plt.axes(projection='3d')
    ax.plot3D(Y[0, :], Y[1, :], Y[2, :])
    ax.set_xlim([-L / 2, L / 2])
    ax.set_ylim([-L / 2, L / 2])
    ax.set_zlim([0, L])
    plt.title('Rod IVP Solution')
    ax.set_xlabel('x (m)')
    ax.set_ylabel('y (m)')
    ax.set_zlabel('z (m)')
    plt.show()


if __name__ == '__main__':
    main()
