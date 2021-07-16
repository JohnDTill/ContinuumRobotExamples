import math
import numpy as np
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

    # Measured base force and moment
    n0 = np.array([0, 1, 0]).T
    m0 = np.array([0, 0, 0]).T

    # Arbitrary base frame assignment
    p0 = np.array([0, 0, 0]).T
    R0 = np.eye(3)

    # Sub-functions
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

    # Numerical Integration
    y0 = np.concatenate([p0.T, np.reshape(R0, 9), n0.T, m0.T]).T  # Combine states into single state vector
    Y = integrate.solve_ivp(rod_ode, (0, L), y0, max_step=0.01).y  # Solve IVP with numerical integration

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
