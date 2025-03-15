import numpy as np
import control as cnt
import VTOLParam as P

class ctrlStateFeedback:
    # dirty derivatives to estimate thetadot
    def __init__(self):
        # Longitudinal dynamics ---------------------------------------------
        tr_lon = 1.5
        zeta_lon  = 0.707

        wn_lon = 2.2 / tr_lon  # natural frequency
        des_char_poly_lon = [1, 2 * zeta_lon * wn_lon, wn_lon**2]
        des_poles_lon = np.roots(des_char_poly_lon)

        # Compute the gains if the system is controllable
        if np.linalg.matrix_rank(cnt.ctrb(P.A_lon, P.B_lon)) != 2:
            print("The system is not controllable")
        else:
            self.K_lon = (cnt.place(P.A_lon, P.B_lon, des_poles_lon))
            self.kr_lon = -1.0 / (P.C_lon @ np.linalg.inv(P.A_lon - P.B_lon @ self.K_lon) @ P.B_lon)
        print('K: ', self.K_lon)
        print('kr: ', self.kr_lon)
        print(des_poles_lon)

        # Lateral dynamics ---------------------------------------------
        #  tuning parameters

        # "Inner Loop"
        tr_lat_th = .2
        wn_lat_th = 2.2 / tr_lat_th  # natural frequency
        zeta_lat_th = 0.707

        # "Outer Loop"
        tr_lat_z = 2
        wn_lat_z = 2.2 / tr_lat_z
        zeta_lat_z = 0.707

        des_char_poly_lat = np.convolve(
            [1, 2 * zeta_lat_z * wn_lat_z, wn_lat_z**2],
            [1, 2 * zeta_lat_th * wn_lat_th, wn_lat_th**2])
        des_poles_lat = np.roots(des_char_poly_lat)

        # Compute the gains if the system is controllable
        if np.linalg.matrix_rank(cnt.ctrb(P.A_lat, P.B_lat)) != 4:
            print("The system is not controllable")
        else:
            self.K_lat = (cnt.place(P.A_lat, P.B_lat, des_poles_lat))
            self.kr_lat = -1.0 / (P.C_lat @ np.linalg.inv(P.A_lat - P.B_lat @ self.K_lat) @ P.B_lat)
        print('K: ', self.K_lat)
        print('kr: ', self.kr_lat)
        print(des_poles_lat)


    def update(self, r, x):
        Href = r[0]
        Zref = r[1]
        # compute force
        x_lon = x[[1,4]]
        Fe = P.g * (P.mc + 2 * P.mr)  # equilibrium force
        F_tilde = -self.K_lon @ x_lon + self.kr_lat * Href
        F = saturate(Fe + F_tilde[0][0], P.F_max)

        # compute torque
        x_lat = x[[0, 2, 3, 5]]
        tau_tilde = -self.K_lat @ x_lat + self.kr_lat * Zref
        tau = saturate(tau_tilde[0][0], P.Tau_max)        

        # compute fr and fl
        fr = (F + tau * P.d) / 2.0
        fl = (F - tau * P.d) / 2.0
        U = np.array([[fr], [fl]])

        return U


def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u


