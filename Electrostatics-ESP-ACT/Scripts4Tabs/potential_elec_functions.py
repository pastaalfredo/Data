import numpy as np
import math

one_4pi_eps0 = 1389.3544561


def slater2_charge(distance, alpha):
    if 0 == distance:
         return alpha/2
    else:
         return 1/distance - (6 + 9*distance*alpha + 6*(alpha*distance)**2 + 2*(distance*alpha)**3) / (6*math.exp(2*distance*alpha)*distance)



def slater_charge(distance, alpha):
    if 0 == distance:
        return alpha
    else:
        return (1/distance) - (1 + distance * alpha) / (math.exp(2 * distance * alpha) * distance)



def gaussian(distance, alpha):
    if 0 == distance:
        return 2 * alpha / np.sqrt(np.pi)
    else:
        return math.erf(distance * alpha) / distance




##############################################################################

#2S-2S- shell shell

def calculate_energy_shellshell2(r, xi, xj):
    rxi = r * xi
    rxj = r * xj

    S = (1 / r) * (
             (6 * np.exp(2 * (rxi + rxj)) * (np.power(np.power(rxi, 2) - np.power(rxj, 2), 7))) -
             (np.exp(2 * rxi) * np.power(rxi, 6) *
             (21 * np.power(rxi, 4) * np.power(rxj, 4) * (6 + 11 * rxj + 2 * np.power(rxj, 2)) -
              2 * np.power(rxj, 8) * (90 + 54 * rxj + 12 * np.power(rxj, 2) + np.power(rxj, 3)) +
              np.power(rxi, 8) * (6 + 9 * rxj + 6 * np.power(rxj, 2) + 2 * np.power(rxj, 3)) +
              np.power(rxi, 2) * np.power(rxj, 6) * (-390 - 69 * rxj + 18 * np.power(rxj, 2) + 4 * np.power(rxj, 3)) -
              np.power(rxi, 6) * np.power(rxj, 2) * (42 + 63 * rxj + 42 * np.power(rxj, 2) + 4 * np.power(rxj, 3)))) +
              (np.exp(2 * rxj) * np.power(rxj, 6) *
              (-24 * np.power(rxi, 10) - 2 * np.power(rxi, 11) - 69 * np.power(rxi, 7) * np.power(rxj, 2) +
              6 * np.power(rxj, 8) + 9 * rxi * np.power(rxj, 8) +
              4 * np.power(rxi, 9) * (-27 + np.power(rxj, 2)) +
              18 * np.power(rxi, 8) * (-10 + np.power(rxj, 2)) +
              6 * np.power(rxi, 2) * np.power(rxj, 6) * (-7 + np.power(rxj, 2)) -
              42 * np.power(rxi, 4) * np.power(rxj, 4) * (-3 + np.power(rxj, 2)) +
              np.power(rxi, 3) * np.power(rxj, 6) * (-63 + 2 * np.power(rxj, 2)) +
              6 * np.power(rxi, 6) * np.power(rxj, 2) * (-65 + 7 * np.power(rxj, 2)) +
              np.power(rxi, 5) * (231 * np.power(rxj, 4) - 4 * np.power(rxj, 6))))) / \
              (6 * np.exp(2 * (rxi + rxj)) * np.power(rxi - rxj, 7) * np.power(rxi + rxj, 7))

    return S


#core-2S
def slater2_core_shell_potential(r, xi):
    S = 1 / r - (6 + 9 * r * xi + 6 * pow(r, 2) * pow(xi, 2) + 2 * pow(r, 3) * pow(xi, 3)) / \
        (6 * math.exp(2 * r * xi) * r)

    return S



#core-1S
def slater_core_shell_potential(distance, alpha):
    return 1/distance - (1 + distance * alpha) / (math.exp(2 * distance * alpha) * distance)



#1S-1S- shell shell

def calculate_energy_shellshell(r, rxi, rxj):
    rxi2 = rxi * rxi
    rxj2 = rxj * rxj
    rxij = rxi + rxj
    exp2rxij = np.exp(2 * rxij)

    n = (exp2rxij * pow((rxi2 - rxj2), 3) +
                 np.exp(2 * rxj) * pow(rxj, 4) *
                 (-3 * rxi2 - pow(rxi, 3) + rxj2 + rxi * rxj2) -
                 np.exp(2 * rxi) * pow(rxi, 4) *
                 (rxi2 * (1 + rxj) - rxj2 * (3 + rxj)))

    d = exp2rxij * pow((rxi - rxj), 3) * pow(rxij, 3)
    energy = (1 / r) * (n / d)

    return energy


#1S-2S shell-shell


def double_Slater_1S_2S(r, xi, xj):
    rxi = r * xi
    rxj = r * xj
    n= (6 * np.exp(2*(rxi + rxj)) * np.power(np.power(rxi, 2) - np.power(rxj, 2), 5) +
         6 * np.exp(2*rxj) * np.power(rxj, 6) *
         (-4*np.power(rxi, 4) - np.power(rxi, 5) - 5*np.power(rxi, 2)*np.power(rxj, 2) +
          np.power(rxj, 4) + rxi*np.power(rxj, 4)) -
         np.exp(2*rxi) * np.power(rxi, 4) *
         (np.power(rxi, 6)*(6 + 9*rxj + 6*np.power(rxj, 2) + 2*np.power(rxj, 3)) -
          3*np.power(rxi, 4)*np.power(rxj, 2) *
          (10 + 15*rxj + 10*np.power(rxj, 2) + 2*np.power(rxj, 3)) +
          3*np.power(rxi, 2)*np.power(rxj, 4) *
          (20 + 33*rxj + 14*np.power(rxj, 2) + 2*np.power(rxj, 3)) -
          np.power(rxj, 6)*(84 + 63*rxj + 18*np.power(rxj, 2) + 2*np.power(rxj, 3))))


    d=  (6 * np.exp(2*(rxi + rxj)) * np.power(rxi - rxj, 5) * np.power(rxi + rxj, 5))

    S = (1 / r) *  (n /d )


    return S


def gaussian(distance:float, zeta:float)->float:
    if distance == 0:
        return 2.0*zeta/math.sqrt(math.pi)
    else:
        return math.erf(distance*zeta)/distance

def gaussian_pointshell(distance, q_c, q_s,z2):
    erf2 = math.erf(distance * z2)
    energy=(q_c *q_s * erf2 / distance)
    return energy



def gaussian_shellshell(distance, q_c, q_s, z1, z2):
    erf = math.erf(distance * z1 * z2 / (z1**2 + z2**2)**0.5)
    energy= (q_c * q_s*erf / distance)
    return energy


#"point core 1 slater shell charge"

def Point_core_1slater_shell(distances, q_c_na, q_s_na, q_c_cl, q_s_cl, z_na, z_cl):
    energy = []
    for i in range(len(distances)):

        energy.append(
            one_4pi_eps0 * (
                (q_c_na * q_c_cl / distances[i] ) +
                (q_s_na*q_s_cl*(calculate_energy_shellshell(distances[i], z_cl*distances[i], z_na*distances[i])) ) +
                (q_s_cl*q_c_na * slater_core_shell_potential(distances[i], z_cl)) +
                (q_s_na*q_c_cl * slater_core_shell_potential(distances[i], z_na))
            )

        )

    return energy


def Point_core_1slater_2slater_shell(distances,
                                     q_c_na, q_s1_na, q_s2_na,
                                     q_c_cl, q_s1_cl, q_s2_cl,
                                     z_c_na, z_c_cl, z_s_na, z_s_cl):
    energy = []
    for i in range(len(distances)):

        energy.append(
            one_4pi_eps0 * (
                 (q_c_na * q_c_cl / distances[i] ) +

                 (q_s1_cl*q_c_na * slater_core_shell_potential(distances[i], z_c_cl)) +
                 (q_s1_na*q_c_cl * slater_core_shell_potential(distances[i], z_c_na))+

                 (q_s2_cl*q_c_na * slater2_core_shell_potential(distances[i], z_s_cl)) +
                 (q_s2_na*q_c_cl * slater2_core_shell_potential(distances[i], z_s_na))+

                 (q_s1_na*q_s1_cl * (calculate_energy_shellshell(distances[i], z_c_cl*distances[i], z_c_na*distances[i])) ) +
                 (q_s2_na*q_s2_cl * (calculate_energy_shellshell2(distances[i], z_s_cl, z_s_na)) )  +
                 (q_s1_na*q_s2_cl * (double_Slater_1S_2S(distances[i], z_c_na, z_s_cl)))+
                 (q_s1_cl*q_s2_na * (double_Slater_1S_2S(distances[i], z_c_cl, z_s_na)))

            )

        )

    return energy



#"point core gaussian gaussian shell charge"
def Point_core_2gaussian_shell(distances, q_c_na, q_s1_na, q_s2_na, q_c_cl, q_s1_cl, q_s2_cl, z1_na, z1_cl, z2_na, z2_cl):
    energy = []
    for i in range(len(distances)):

        energy.append(
             (one_4pi_eps0 *
                ((q_c_na * q_c_cl / distances[i])  +
                (gaussian_shellshell(distances[i], q_s1_cl, q_s1_na, z1_cl, z1_na)) +
                (gaussian_shellshell(distances[i], q_s2_cl, q_s2_na, z2_cl, z2_na)) +
                (gaussian_shellshell(distances[i], q_s1_cl, q_s2_na, z1_cl, z2_na)) +
                (gaussian_shellshell(distances[i], q_s2_cl, q_s1_na, z2_cl, z1_na)) +
                (gaussian_pointshell(distances[i], q_c_cl, q_s1_na, z1_na)) +
                (gaussian_pointshell(distances[i], q_c_cl, q_s2_na, z2_na)) +
                (gaussian_pointshell(distances[i], q_c_na, q_s1_cl, z1_cl)) +
                (gaussian_pointshell(distances[i], q_c_na, q_s2_cl, z2_cl)) )
            )

        )

    return energy


#"point core + gaussian shell charge"

def Point_core_gaussian_shell(distances, q_c_na, q_s1_na, q_c_cl, q_s1_cl, z1_na, z1_cl):
    energy = []
    for i in range(len(distances)):

        energy.append(
             one_4pi_eps0 *
                (  (q_c_na * q_c_cl / distances[i]) +
                (gaussian_shellshell(distances[i], q_s1_cl, q_s1_na, z1_cl, z1_na)) +

                (gaussian_pointshell(distances[i], q_c_cl, q_s1_na, z1_na)) +

                (gaussian_pointshell(distances[i], q_c_na, q_s1_cl, z1_cl))
                )
            )

    return energy
