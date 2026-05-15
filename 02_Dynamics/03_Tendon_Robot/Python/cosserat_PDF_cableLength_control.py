import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from numba import njit
import time
import os

# =====================================================================
# Core low-level calculation functions compiled with Numba JIT (Extreme Acceleration)
# =====================================================================

# =====================================================================
# Numba JIT Core functions (Fixed memory contiguity issue)
# =====================================================================

@njit(cache=True)
def hat(y):
    res = np.zeros((3, 3))
    res[0, 1] = -y[2]; res[0, 2] = y[1]
    res[1, 0] = y[2];  res[1, 2] = -y[0]
    res[2, 0] = -y[1]; res[2, 1] = y[0]
    return res

@njit(cache=True)
def staticODE_njit(y, tau_guess, r, Kse, Kbt, rhoAg, num_tendons):
    # [Fix]: Force conversion to a contiguous array before reshape
    R = np.ascontiguousarray(y[3:12]).reshape((3, 3)) 
    v = y[12:15]
    u = y[15:18]
    pib_s_norm = np.zeros(num_tendons)
    
    a = np.zeros(3); b = np.zeros(3)
    A = np.zeros((3, 3)); G = np.zeros((3, 3)); H = np.zeros((3, 3))

    for i in range(num_tendons):
        pb_si = np.cross(u, r[i]) + v
        pib_s_norm[i] = np.linalg.norm(pb_si)
        hat_pb = hat(pb_si)
        A_i = - (hat_pb @ hat_pb) * (tau_guess[i] / (pib_s_norm[i]**3 + 1e-12))
        G_i = -A_i @ hat(r[i])
        a_i = A_i @ np.cross(u, pb_si)
        
        a += a_i; b += np.cross(r[i], a_i)
        A += A_i; G += G_i; H += hat(r[i]) @ G_i

    nb = Kse @ (v - np.array([0.0, 0.0, 1.0]))
    mb = Kbt @ u

    rhs1 = -np.cross(u, nb) - R.T @ rhoAg - a
    rhs2 = -np.cross(u, mb) - np.cross(v, nb) - b
    rhs = np.empty(6)
    rhs[0:3] = rhs1
    rhs[3:6] = rhs2

    K = np.zeros((6, 6))
    K[0:3, 0:3] = A + Kse; K[0:3, 3:6] = G
    K[3:6, 0:3] = G.T;     K[3:6, 3:6] = H + Kbt

    vs_and_us = np.linalg.solve(K, rhs)

    p_s = R @ v
    R_s = R @ hat(u)

    y_s = np.zeros(24 + num_tendons)
    y_s[0:3] = p_s
    y_s[3:12] = R_s.flatten() 
    y_s[12:18] = vs_and_us
    y_s[24:24+num_tendons] = pib_s_norm
    return y_s

@njit(cache=True)
def integrate_static(Y_init, Z_init, ds, N, tau_guess, r, Kse, Kbt, rhoAg, num_tendons):
    Y = Y_init.copy()
    Z = Z_init.copy()
    for j in range(N - 1):
        y_s = staticODE_njit(Y[:, j], tau_guess, r, Kse, Kbt, rhoAg, num_tendons)
        Y[:, j+1] = Y[:, j] + ds * y_s
        Z[0:6, j] = Y[12:18, j]
        Z[12:18, j] = y_s[12:18]
    return Y, Z

@njit(cache=True)
def dynamicODE_njit(y, z_h, tau_guess, r, Kse, Kbt, Bse, Bbt, rhoA, rhoAg, rho, J, C_mat, c0, num_tendons):
    # [Fix]: Force conversion to a contiguous array before reshape
    R = np.ascontiguousarray(y[3:12]).reshape((3, 3))
    v, u, q, w = y[12:15], y[15:18], y[18:21], y[21:24]
    v_h, u_h, q_h, w_h, v_sh, u_sh = z_h[0:3], z_h[3:6], z_h[6:9], z_h[9:12], z_h[12:15], z_h[15:18]

    a = np.zeros(3); b = np.zeros(3)
    A = np.zeros((3, 3)); G = np.zeros((3, 3)); H = np.zeros((3, 3))
    pib_s_norm = np.zeros(num_tendons)

    for i in range(num_tendons):
        pb_si = np.cross(u, r[i]) + v
        pib_s_norm[i] = np.linalg.norm(pb_si)
        hat_pb = hat(pb_si)
        A_i = - (hat_pb @ hat_pb) * (tau_guess[i] / (pib_s_norm[i]**3 + 1e-12))
        G_i = -A_i @ hat(r[i])
        a_i = A_i @ np.cross(u, pb_si)
        
        a += a_i; b += np.cross(r[i], a_i)
        A += A_i; G += G_i; H += hat(r[i]) @ G_i

    v_t = c0 * v + v_h; u_t = c0 * u + u_h
    q_t = c0 * q + q_h; w_t = c0 * w + w_h

    nb = Kse @ (v - np.array([0.0, 0.0, 1.0])) + Bse @ v_t
    mb = Kbt @ u + Bbt @ u_t

    rhs1 = -a + rhoA * (np.cross(w, q) + q_t) + (C_mat @ q) * np.linalg.norm(q) - R.T @ rhoAg - np.cross(u, nb) - Bse @ v_sh
    rhs2 = -b + rho * np.cross(w, J @ w) + rho * (J @ w_t) - np.cross(v, nb) - np.cross(u, mb) - Bbt @ u_sh

    rhs = np.empty(6)
    rhs[0:3] = rhs1; rhs[3:6] = rhs2

    K = np.zeros((6, 6))
    K[0:3, 0:3] = A + Kse + c0 * Bse; K[0:3, 3:6] = G
    K[3:6, 0:3] = G.T;                K[3:6, 3:6] = H + Kbt + c0 * Bbt

    vs_and_us = np.linalg.solve(K, rhs)

    p_s = R @ v; R_s = R @ hat(u)
    q_s = v_t - hat(u) @ q + hat(w) @ v; w_s = u_t - hat(u) @ w

    y_s = np.zeros(24 + num_tendons)
    y_s[0:3] = p_s; y_s[3:12] = R_s.flatten()
    y_s[12:15] = vs_and_us[0:3]; y_s[15:18] = vs_and_us[3:6]
    y_s[18:21] = q_s; y_s[21:24] = w_s
    y_s[24:24+num_tendons] = pib_s_norm

    z = np.zeros(18)
    z[0:3] = v; z[3:6] = u; z[6:9] = q; z[9:12] = w
    z[12:15] = vs_and_us[0:3]; z[15:18] = vs_and_us[3:6]

    return y_s, z

@njit(cache=True)
def integrate_dynamic(Y_init, Z_init, Z_h, ds, N, tau_guess, r, Kse, Kbt, Bse, Bbt, rhoA, rhoAg, rho, J, C_mat, c0, num_tendons):
    Y = Y_init.copy()
    Z = Z_init.copy()
    for j in range(N - 1):
        y_s, z = dynamicODE_njit(Y[:, j], Z_h[:, j], tau_guess, r, Kse, Kbt, Bse, Bbt, rhoA, rhoAg, rho, J, C_mat, c0, num_tendons)
        Y[:, j+1] = Y[:, j] + ds * y_s
        Z[:, j] = z
    return Y, Z


# =====================================================================
# Main Robot Control Class
# =====================================================================

class TendonRobotPDE:
    def __init__(self, enable_plot=False):
        self.enable_plot = enable_plot 
        
        # --- Physical and System Parameters ---
        self.L = 0.5            
        self.N = 200            
        self.E = 207e9          
        self.shear_modulus = self.E / (2 * (1 + 0.3)) 
        self.radiu = 0.001      
        self.total_mass = 0.034 
        
        self.g = np.array([-9.81, 0, 0], dtype=np.float64) 
        self.Bse = np.zeros((3, 3), dtype=np.float64)      
        self.Bbt = np.diag(np.array([5e-4, 5e-4, 5e-4])) 
        self.C = np.diag(np.array([1e-4, 1e-4, 1e-4]))   
        
        self.num_tendons = 1
        self.compliance = np.array([1e-4])
        self.num_disks = 9
        self.tendon_offset = 0.01506
        self.r_arr = np.array([self.tendon_offset * np.array([np.cos(0), np.sin(0), 0])], dtype=np.float64) 
        
        self.area = np.pi * self.radiu**2
        self.J = np.diag(np.array([np.pi * self.radiu**4 / 4, np.pi * self.radiu**4 / 4, np.pi * self.radiu**4 / 2]))
        
        self.Kse = np.diag(np.array([self.shear_modulus * self.area, self.shear_modulus * self.area, self.E * self.area]))
        self.Kbt = np.diag(np.array([self.E * self.J[0, 0], self.E * self.J[1, 1], self.shear_modulus * self.J[2, 2]]))
        
        self.rho = self.total_mass / (self.L * self.area)
        self.rhoA = self.rho * self.area
        self.rhoAg = self.rhoA * self.g
        
        self.q0 = np.zeros(3); self.w0 = np.zeros(3)
        self.R0 = np.eye(3)
        
        self.T = 10.0
        self.dt = 1e-2
        self.alpha = 0.0
        self.t = 0.0
        self.STEPS = int(self.T / self.dt)
        self.ds = self.L / (self.N - 1)
        
        self.c0 = (1.5 + self.alpha) / (self.dt * (1 + self.alpha))
        self.c1 = -2 / self.dt
        self.c2 = (0.5 + self.alpha) / (self.dt * (1 + self.alpha))
        
        self.zA = 0.0; self.zB = -0.01619
        self.t1, self.t2, self.t3, self.t4 = 1.652, 1.967, 5.63, 5.94
        
        y_size = 24 + self.num_tendons
        self.Y = np.zeros((y_size, self.N))
        self.Y[2, :] = np.linspace(0, self.L, self.N) 
        R0_flat = self.R0.flatten() 
        for i in range(self.N):
            self.Y[3:12, i] = R0_flat
            
        self.Z = np.zeros((18, self.N))
        self.Z[2, :] = 1.0 
        
        self.Y_prev = np.copy(self.Y); self.Z_prev = np.copy(self.Z)
        self.Z_h = np.zeros_like(self.Z)
        self._Y_temp_latest = np.zeros_like(self.Y)
        self._Z_temp_latest = np.zeros_like(self.Z)
        self.running = True

        # --- Data Recording Containers ---
        self.history_time = []
        self.history_tip_pos = []
        self.history_backbone = []
        self.history_step_time = []

    def Z_t(self, t):
        if self.t1 < t <= self.t2:
            return self.zA + (self.zB - self.zA) * (t - self.t1) / (self.t2 - self.t1)
        elif self.t2 < t <= self.t3: return self.zB
        elif self.t3 < t <= self.t4:
            return self.zB + (self.zA - self.zB) * (t - self.t3) / (self.t4 - self.t3)
        return self.zA

    def staticBVP(self, guess):
        v0 = np.linalg.solve(self.Kse, guess[0:3]) + np.array([0.0, 0.0, 1.0])
        tau_val = guess[6:6+self.num_tendons]
        tau = np.maximum(tau_val, 0)
        slack = -np.minimum(tau_val, 0)

        Y_temp = np.copy(self.Y); Z_temp = np.copy(self.Z)
        Y_temp[:, 0] = np.concatenate([Y_temp[0:12, 0], v0, guess[3:6], self.q0, self.w0, -self.Z_t(self.t) * np.ones(self.num_tendons)])
        
        Y_temp, Z_temp = integrate_static(Y_temp, Z_temp, self.ds, self.N, tau, self.r_arr, self.Kse, self.Kbt, self.rhoAg, self.num_tendons)

        self._Y_temp_latest = Y_temp; self._Z_temp_latest = Z_temp
        vL = Y_temp[12:15, -1]; uL = Y_temp[15:18, -1]
        
        force_error = -(self.Kse @ (vL - np.array([0.0, 0.0, 1.0])))
        moment_error = -(self.Kbt @ uL)

        for idx in range(self.num_tendons):
            pb_si = np.cross(uL, self.r_arr[idx]) + vL
            Fb_i = -tau[idx] * pb_si / np.linalg.norm(pb_si)
            force_error += Fb_i
            moment_error += np.cross(self.r_arr[idx], Fb_i)

        length_error = (Y_temp[24:24+self.num_tendons, -1] + slack) - (self.L + self.L * (self.compliance * tau))
        return np.concatenate([force_error, moment_error, length_error])

    def dynamicBVP(self, guess):
        v0 = np.linalg.solve(self.Kse, guess[0:3]) + np.array([0.0, 0.0, 1.0])
        tau_val = guess[6:6+self.num_tendons]
        tau = np.maximum(tau_val, 0)
        slack = -np.minimum(tau_val, 0)

        Y_temp = np.copy(self.Y); Z_temp = np.copy(self.Z)
        Y_temp[:, 0] = np.concatenate([Y_temp[0:12, 0], v0, guess[3:6], self.q0, self.w0, -self.Z_t(self.t) * np.ones(self.num_tendons)])

        Y_temp, Z_temp = integrate_dynamic(Y_temp, Z_temp, self.Z_h, self.ds, self.N, tau, self.r_arr, 
                                           self.Kse, self.Kbt, self.Bse, self.Bbt, self.rhoA, self.rhoAg, 
                                           self.rho, self.J, self.C, self.c0, self.num_tendons)

        self._Y_temp_latest = Y_temp; self._Z_temp_latest = Z_temp
        vL = Y_temp[12:15, -1]; uL = Y_temp[15:18, -1]
        vL_t = self.c0 * vL + self.Z_h[0:3, -1]; uL_t = self.c0 * uL + self.Z_h[3:6, -1]

        force_error = -(self.Kse @ (vL - np.array([0.0, 0.0, 1.0])) + self.Bse @ vL_t)
        moment_error = -(self.Kbt @ uL + self.Bbt @ uL_t)

        for idx in range(self.num_tendons):
            pb_si = np.cross(uL, self.r_arr[idx]) + vL
            Fb_i = -tau[idx] * pb_si / np.linalg.norm(pb_si)
            force_error += Fb_i
            moment_error += np.cross(self.r_arr[idx], Fb_i)

        length_error = Y_temp[24:24+self.num_tendons, -1] + slack - (self.L + self.L * (self.compliance * tau))
        return np.concatenate([force_error, moment_error, length_error])

    # ------------------ Plot Export System ------------------
    def export_data(self):
        """Plot recorded data and export to the output folder"""
        output_dir = "output"
        os.makedirs(output_dir, exist_ok=True)
        print(f"\n>> Computation complete. Exporting charts to the '{output_dir}/' folder...")

        time_arr = np.array(self.history_time)
        tip_arr = np.array(self.history_tip_pos)

        # Plot 1: Tip trajectory over time
        fig, axs = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
        fig.suptitle('Robot Tip Trajectory over Time', fontsize=16)
        
        axs[0].plot(time_arr, tip_arr[:, 0], 'r-', linewidth=2)
        axs[0].set_ylabel('Tip X (m)'); axs[0].grid(True)
        
        axs[1].plot(time_arr, tip_arr[:, 1], 'g-', linewidth=2)
        axs[1].set_ylabel('Tip Y (m)'); axs[1].grid(True)
        
        axs[2].plot(time_arr, tip_arr[:, 2], 'b-', linewidth=2)
        axs[2].set_ylabel('Tip Z (m)'); axs[2].set_xlabel('Time (s)'); axs[2].grid(True)
        
        fig.tight_layout()
        traj_path = os.path.join(output_dir, "tip_trajectory.png")
        plt.savefig(traj_path, dpi=300)
        plt.close(fig)

        # Plot 2: 3D shape evolution sampled at key frames
        fig_3d = plt.figure(figsize=(10, 8))
        ax_3d = fig_3d.add_subplot(111, projection='3d')
        ax_3d.set_title('Robot Backbone Shape Evolution')
        
        # Uniformly sample 5-10 frames to draw on the same canvas
        num_frames = min(10, len(self.history_backbone))
        indices = np.linspace(0, len(self.history_backbone)-1, num_frames, dtype=int)
        
        colors = plt.cm.viridis(np.linspace(0, 1, num_frames))
        for idx, color in zip(indices, colors):
            backbone = self.history_backbone[idx]
            t_val = self.history_time[idx]
            ax_3d.plot(backbone[0, :], backbone[1, :], backbone[2, :], 
                       color=color, linewidth=2, alpha=0.8, label=f't={t_val:.1f}s')

        ax_3d.set_xlabel('X (m)'); ax_3d.set_ylabel('Y (m)'); ax_3d.set_zlabel('Z (m)')
        ax_3d.set_xlim([-0.05, 0.05]); ax_3d.set_ylim([-0.05, 0.05]); ax_3d.set_zlim([0, self.L + 0.05])
        ax_3d.legend(loc='upper left', bbox_to_anchor=(1.05, 1))
        
        fig_3d.tight_layout()
        shape_path = os.path.join(output_dir, "backbone_shapes.png")
        plt.savefig(shape_path, dpi=300)
        plt.close(fig_3d)

        # Plot 3: Step computation time
        fig_time, ax_time = plt.subplots(figsize=(10, 5))
        ax_time.plot(time_arr[1:], self.history_step_time, 'm-', linewidth=1.5)
        ax_time.set_title('Computation Time per Simulation Step', fontsize=14)
        ax_time.set_xlabel('Simulation Time (s)')
        ax_time.set_ylabel('Computation Time (s)')
        ax_time.grid(True)
        
        fig_time.tight_layout()
        time_path = os.path.join(output_dir, "step_computation_time.png")
        plt.savefig(time_path, dpi=300)
        plt.close(fig_time)

        print(f">> Charts exported successfully:\n  - {traj_path}\n  - {shape_path}\n  - {time_path}")

    def run(self):
        print(">> Starting initial static boundary value problem solving...")
        t_start_static = time.time()
        guessPre = fsolve(self.staticBVP, np.zeros(6 + self.num_tendons))
        self.Y = np.copy(self._Y_temp_latest);   self.Y_prev = np.copy(self.Y)
        self.Z = np.copy(self._Z_temp_latest);   self.Z_prev = np.copy(self.Z)
        print(f">> Static solving complete. Time elapsed: {time.time() - t_start_static:.4f}s")
        
        # Record initial state
        self.history_time.append(self.t)
        self.history_tip_pos.append(self.Y[0:3, -1].copy())
        self.history_backbone.append(self.Y[0:3, :].copy())

        print(">> Starting dynamic simulation...")
        if not self.enable_plot:
            print(">> [INFO] Plotting disabled. Running full-speed background computation. Will output step duration...")
            
        for tStep in range(self.STEPS):
            step_t_start = time.time() # Start step timer

            self.Z_h = self.c1 * self.Z + self.c2 * self.Z_prev
            self.Y_prev = np.copy(self.Y)
            self.Z_prev = np.copy(self.Z)
            
            guessPre = fsolve(self.dynamicBVP, guessPre)
            
            self.Y = np.copy(self._Y_temp_latest)
            self.Z = np.copy(self._Z_temp_latest)
            
            step_time = time.time() - step_t_start # End step timer

            self.t += self.dt
            
            # Record data
            self.history_time.append(self.t)
            self.history_tip_pos.append(self.Y[0:3, -1].copy())
            self.history_backbone.append(self.Y[0:3, :].copy())
            self.history_step_time.append(step_time)

            # Print step time and progress
            progress = (tStep + 1) / self.STEPS * 100
            print(f"Step {tStep+1:03d}/{self.STEPS} | Progress: {progress:5.1f}% | t = {self.t:.2f}s | Elapsed: {step_time:.4f}s")
            
        # Loop complete, export charts
        self.export_data()

if __name__ == "__main__":
    # Default to plot disabled, running full-speed background calculation and exporting images
    robot = TendonRobotPDE(enable_plot=False)
    robot.run()