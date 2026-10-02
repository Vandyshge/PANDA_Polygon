import numpy as np
from numpy import pi, cos, sin, sqrt, abs

BIG_NUM = np.inf


# ============================================================
# Common geometry of the square pore
#
# Dimensionless convention:
#     square side a = 1
#
# z = 0 is the centre of the square cross-section.
#
# Effective wall-layer correction:
#     a_eff = 1 - 2*delta
#     z_bottom = -a_eff/2
#     z_top    =  a_eff/2
# ============================================================


def a_eff_square(delta):
    """Effective side a' = 1 - 2*delta."""
    return 1.0 - 2.0 * delta


def a_prime_square(delta):
    """Alias for a_eff_square."""
    return a_eff_square(delta)


def z_bottom_square(delta):
    """Bottom accessible z coordinate: -a'/2."""
    return -a_eff_square(delta) / 2.0


def z_top_square(delta):
    """Top accessible z coordinate: a'/2."""
    return a_eff_square(delta) / 2.0
# ============================================================
# Density profiles
# ============================================================


""" filament_w0c0 """


def r_filament_w0c0_square(l, phi, th, delta):
    """
    Radius of the axial cylindrical phase.

    pi r^2 = phi   (в сечении)
    """
    return sqrt(phi  / pi)


def rho_filament_w0c0_square(z, l, ly, phi, th, delta):
    """
    Density profile for filament_w0c0.

    Цилиндр по оси поры, в сечении круг радиуса r.
    z = 0 — центр поры.
    """
    z = np.asarray(z, dtype=float)

    r = r_filament_w0c0_square(l, phi, th, delta)

    mask = abs(z) <= r
    y = sqrt(np.maximum(0.0, r ** 2 - z ** 2))

    return 2.0 * y * mask / ly


""" filament_w0c0_straight """


def b_filament_w0c0_straight_square(l, phi, th, delta):
    """
    Side of the square cross-section of the layer.

    b^2 = phi 
    """
    return sqrt(phi)


def rho_filament_w0c0_straight_square(z, l, ly, phi, th, delta):
    """
    Density profile for filament_w0c0_straight.

    Квадратный шнур в сечении, сторона b.
    z = 0 — центр поры.
    """
    z = np.asarray(z, dtype=float)

    b = b_filament_w0c0_straight_square(l, phi, th, delta)

    mask = abs(z) <= b / 2.0

    return b * mask / ly


""" filament_w1c0 """


def r_filament_w1c0_square(l, phi, th, delta):
    """
    Radius of the circular arc for filament_w1c0.

    S_seg = r^2 (theta - sin(theta) cos(theta))
    V = l * S_seg = phi 
    """
    denom = th - sin(th) * cos(th)
    return sqrt(phi/ denom)


def rho_filament_w1c0_square(z, l, ly, phi, th, delta):
    """
    Density profile for filament_w1c0.

    Цилиндрический сегмент, опирающийся на нижнюю стенку.
    Центр окружности на высоте z_center = z_bottom + r*cos(theta').
    """
    z = np.asarray(z, dtype=float)

    r = r_filament_w1c0_square(l, phi, th, delta)
    th_p = pi - th

    z_bottom = z_bottom_square(delta)
    z_center = z_bottom + r * cos(th_p)

    mask = (z >= z_bottom) & (z <= z_center + r)

    r_sq = np.maximum(0.0, r ** 2 - (z - z_center) ** 2)
    y = sqrt(r_sq)

    return 2.0 * y * mask / ly
# ============================================================
# Auxiliary quantities for empty-corner worm configurations
# ============================================================


def A_square(theta):
    """
    Фактор площади одного пустого угла A(theta).

    Из tex:
        S_ang = r^2 (sqrt(2) sin(theta + pi/4) cos(theta) + 3pi/4 - theta)
              ≡ A(theta) r^2
    """
    return (
        sqrt(2.0) * sin(theta + pi / 4.0) * cos(theta)
        + 3.0 * pi / 4.0
        - theta
    )





def L1_square(r, theta):
    """
    Расстояние L_1 от вершины угла до ближней точки контакта.

    Из tex:
        L_1 = -sqrt(2) r sin(theta + pi/4)
    """
    return -sqrt(2.0) * r * sin(theta + pi / 4.0)


def L2_square(r, theta):
  
    return sqrt(2.0) * r * sin(theta - pi / 4.0)




""" filament_w4c0 """


def r_filament_w4c0_square(l, phi, th, delta):

    a_eff = a_eff_square(delta)
    return sqrt((a_eff ** 2 - phi) / (4.0 * A_square(th)))


def rho_filament_w4c0_square(z, l, ly, phi, th, delta):

    z = np.asarray(z, dtype=float)

    a_eff = a_eff_square(delta)
    r = r_filament_w4c0_square(l, phi, th, delta)
    th_p = pi - th

    L1 = L1_square(r, th)

    z_bottom = z_bottom_square(delta)
    z_top = z_top_square(delta)

    z_center_low = z_bottom + r * cos(th_p)
    z_center_up = z_top - r * cos(th_p)

    z1 = z_bottom
    z2 = z_bottom + L1
    z3 = z_top - L1
    z4 = z_top
    rho = np.zeros_like(z, dtype=float)

    # нижняя закруглённая часть
    mask1 = (z >= z1) & (z <= z2)
    rho[mask1] = (
        a_eff
        - 2.0 * sqrt(np.maximum(0.0, r ** 2 - (z[mask1] - z_center_low) ** 2))
    ) / ly

    # центральная полоса
    mask2 = (z > z2) & (z <= z3)
    rho[mask2] = a_eff / ly

    # верхняя закруглённая часть
    mask3 = (z > z3) & (z <= z4)
    rho[mask3] = (
        a_eff
        - 2.0 * sqrt(np.maximum(0.0, r ** 2 - (z[mask3] - z_center_up) ** 2))
    ) / ly

    return rho


""" filament_w4c1 """


def r_filament_w4c1_square(l, phi, th, delta):
    """


    """
    a_eff = a_eff_square(delta)
    return sqrt((a_eff ** 2 - phi) / (3.0 * A_square(th)))


def rho_filament_w4c1_square(z, l, ly, phi, th, delta):
    """
    Density profile for filament_w4c1.

    """
    z = np.asarray(z, dtype=float)

    a_eff = a_eff_square(delta)
    r = r_filament_w4c1_square(l, phi, th, delta)
    th_p = pi - th

    L1 = L1_square(r, th)

    z_bottom = z_bottom_square(delta)
    z_top = z_top_square(delta)

    z_center_low = z_bottom + r * cos(th_p)
    z_center_up = z_top - r * cos(th_p)

    z1 = z_bottom
    z2 = z_bottom + L1
    z3 = z_top - L1
    z4 = z_top
    rho = np.zeros_like(z, dtype=float)

    # нижняя часть: один пустой угол
    mask1 = (z >= z1) & (z <= z2)
    rho[mask1] = (
        a_eff
        - sqrt(np.maximum(0.0, r ** 2 - (z[mask1] - z_center_low) ** 2))
    ) / ly

    # центральная полоса
    mask2 = (z > z2) & (z <= z3)
    rho[mask2] = a_eff / ly

    # верхняя часть: два пустых угла
    mask3 = (z > z3) & (z <= z4)
    rho[mask3] = (
        a_eff
        - 2.0 * sqrt(np.maximum(0.0, r ** 2 - (z[mask3] - z_center_up) ** 2))
    ) / ly

    return rho


""" filament_w4c2 """


def r_filament_w4c2_square(l, phi, th, delta):
    """

    """
    a_eff = a_eff_square(delta)
    return sqrt((a_eff ** 2 - phi) / (2.0 * A_square(th)))


def rho_filament_w4c2_square(z, l, ly, phi, th, delta):
    """
    Density profile for filament_w4c2.

    """
    z = np.asarray(z, dtype=float)

    a_eff = a_eff_square(delta)
    r = r_filament_w4c2_square(l, phi, th, delta)
    th_p = pi - th

    L1 = L1_square(r, th)

    z_bottom = z_bottom_square(delta)
    z_top = z_top_square(delta)


    z_center_up = z_top - r * cos(th_p)

    z1 = z_bottom
    z2 = z_bottom + L1
    z3 = z_top - L1
  

    rho = np.zeros_like(z, dtype=float)

    # центральная полоса
    mask1 = (z >= z1) & (z <= z2)
    rho[mask1] = a_eff / ly

    # верхняя часть: два пустых угла
    mask2 = (z > z2) & (z <= z3)
    rho[mask2] = (
        a_eff
        - 2.0 * sqrt(np.maximum(0.0, r ** 2 - (z[mask2] - z_center_up) ** 2))
    ) / ly

    return rho
""" filament_w4c3 """


def r_filament_w4c3_square(l, phi, th, delta):
    """

    """
    a_eff = a_eff_square(delta)
    return sqrt((a_eff ** 2 - phi) / A_square(th))


def rho_filament_w4c3_square(z, l, ly, phi, th, delta):
    """
 
    """
    z = np.asarray(z, dtype=float)

    a_eff = a_eff_square(delta)
    r = r_filament_w4c3_square(l, phi, th, delta)
    th_p = pi - th

    L1 = L1_square(r, th)

    z_bottom = z_bottom_square(delta)
    z_top = z_top_square(delta)

    z_center_up = z_top - r * cos(th_p)

    z1 = z_bottom
    z2 = z_top - L1
    z3 = z_top

    rho = np.zeros_like(z, dtype=float)

    # центральная полоса
    mask1 = (z >= z1) & (z <= z2)
    rho[mask1] = a_eff / ly

    # верхняя часть: один пустой угол
    mask2 = (z > z2) & (z <= z3)
    rho[mask2] = (
        a_eff
        - sqrt(np.maximum(0.0, r ** 2 - (z[mask2] - z_center_up) ** 2))
    ) / ly

    return rho

""" filament_w3c2 """
def alpha_square(theta):
    """
    Угол alpha, опирающийся на пустой угол квадрата.

    Из tex:
        alpha = pi/2 - 2*theta',  где theta' = pi - theta
    """
    theta_p = pi - theta
    return pi  - 2.0 * theta_p

def r_filament_w3c2_square(l, phi, th, delta):
    """
    Радиус сегмента окружности.

    Из tex:
        r = (1/2 - delta) / cos(theta')
    """
    th_p = pi - th
    return (0.5 - delta) / cos(th_p)


def h_filament_w3c2_square(l, phi, th, delta):
    """
    Расстояние между хордой сектора и дном поры.

  
        => h = (phi - 0.5 r^2 (alpha - sin alpha)) / (1 - 2 delta)
    """
    a_eff = a_eff_square(delta)
    r = r_filament_w3c2_square(l, phi, th, delta)
    alpha = alpha_square(th)
    return (phi - 0.5 * (r ** 2) * (alpha - sin(alpha))) / a_eff


def rho_filament_w3c2_square(z, l, ly, phi, th, delta):
    """
    Density profile for filament_w3c2.

    """
    z = np.asarray(z, dtype=float)

    a_eff = a_eff_square(delta)
    r = r_filament_w3c2_square(l, phi, th, delta)
    alpha = alpha_square(th)
    h = h_filament_w3c2_square(l, phi, th, delta)

    z_bottom = z_bottom_square(delta)
    z_center_low = z_bottom + h - r * cos(alpha / 2.0)

    z1 = z_bottom
    z2 = z_bottom + h
    z3 = z_center_low + r

    rho = np.zeros_like(z, dtype=float)

    # плоская часть
    mask1 = (z >= z1) & (z <= z2)
    rho[mask1] = a_eff / ly

    # закруглённая часть
    mask2 = (z > z2) & (z <= z3)
    rho[mask2] = (
        2.0 * sqrt(np.maximum(0.0, r ** 2 - (z[mask2] - z_center_low) ** 2))
    ) / ly

    return rho


""" filament_w2c0_opposite """


def r_filament_w2c0_opposite_square(l, phi, th, delta):
    """
    Радиус сегментов окружности.

    Из tex:
        r = (1/2 - delta) / cos(theta')
    """
    th_p = pi - th
    return (0.5 - delta) / cos(th_p)


def h_filament_w2c0_opposite_square(l, phi, th, delta):
    """
    Расстояние между хордами секторов.


    """
    a_eff = a_eff_square(delta)
    r = r_filament_w2c0_opposite_square(l, phi, th, delta)
    alpha = alpha_square(th)
    return (phi - (r ** 2) * (alpha - sin(alpha))) / a_eff


def rho_filament_w2c0_opposite_square(z, l, ly, phi, th, delta):
    """
    Density profile for filament_w2c0_opposite.

    """
    z = np.asarray(z, dtype=float)

    a_eff = a_eff_square(delta)
    r = r_filament_w2c0_opposite_square(l, phi, th, delta)
    alpha = alpha_square(th)
    h = h_filament_w2c0_opposite_square(l, phi, th, delta)

    z_center_low = h/2.0 - r * cos(alpha / 2.0)
    z_center_high = -h / 2.0 + r * cos(alpha / 2.0)

    z1 = z_center_high - r
    z2 = -h / 2.0
    z3 = h / 2.0
    z4 = z_center_low + r

    rho = np.zeros_like(z, dtype=float)

    # верхняя закруглённая часть
    mask1 = (z >= z1) & (z <= z2)
    rho[mask1] = (
        2.0 * sqrt(np.maximum(0.0, r ** 2 - (z[mask1] - z_center_high) ** 2))
    ) / ly

    # центральная полоса
    mask2 = (z > z2) & (z <= z3)
    rho[mask2] = a_eff / ly

    # нижняя закруглённая часть
    mask3 = (z > z3) & (z <= z4)
    rho[mask3] = (
        2.0 * sqrt(np.maximum(0.0, r ** 2 - (z[mask3] - z_center_low) ** 2))
    ) / ly

    return rho

""" filament_w2c1 """


def r_filament_w2c1_square(l, phi, th, delta):
    """
    Радиус дуги для filament_w2c1.

    Из tex:
        S_liquid = r^2 (theta - pi/4 + cos^2(theta) - sin(theta) cos(theta))
        V = S_liquid L = phi L
        r = sqrt( phi / (theta - pi/4 + cos^2 theta - sin theta cos theta) )
    """
    denom = th - pi / 4.0 + cos(th) ** 2 - sin(th) * cos(th)
    return sqrt(phi / denom)


def rho_filament_w2c1_square(z, l, ly, phi, th, delta):
    """
    Density profile for filament_w2c1.

    Из tex (density_profiles):
        z_c = r cos(theta') - 1/2 + delta

        [z_bottom, z_bottom + L2]:
            ( r cos(theta') + sqrt(r^2 - (z - z_c)^2) ) / ly
        [z_bottom + L2, z_bottom + r(1 + cos(theta'))]:
            2 sqrt(r^2 - (z - z_c)^2) / ly
        иначе 0.

    Замечание: в tex в формуле для z_c написано
        z_c = r cos(theta') - 1/2 + delta,
    что совпадает с z_bottom + r cos(theta').
    """
    z = np.asarray(z, dtype=float)

    r = r_filament_w2c1_square(l, phi, th, delta)
    th_p = pi - th

    z_bottom = z_bottom_square(delta)

    z_c = z_bottom + r * cos(th_p)

    L2 = L2_square(r, th)

    z1 = z_bottom
    z2 = z_bottom + L2
    z3 = z_bottom + r * (1.0 + cos(th_p))

    rho = np.zeros_like(z, dtype=float)

    # нижняя часть: одна стенка смочена
    mask1 = (z >= z1) & (z <= z2)
    rho[mask1] = (
        r * cos(th_p)
        + sqrt(np.maximum(0.0, r ** 2 - (z[mask1] - z_c) ** 2))
    ) / ly

    # верхняя часть: симметричный сегмент
    mask2 = (z > z2) & (z <= z3)
    rho[mask2] = (
        2.0 * sqrt(np.maximum(0.0, r ** 2 - (z[mask2] - z_c) ** 2))
    ) / ly

    return rho
""" droplet_w1c0 """


def r_droplet_w1c0_square(l, phi, th, delta):
    """
    Радиус сферы-сегмента для droplet_w1c0.


    """
    vol_factor = 2.0 - 3.0 * cos(th) + cos(th) ** 3
    return (3.0 * phi * l / (pi * vol_factor)) ** (1.0 / 3.0)


def rho_droplet_w1c0_square(z, l, ly, phi, th, delta):
    """
    Density profile for droplet_w1c0.

  
    """
    z = np.asarray(z, dtype=float)

    r = r_droplet_w1c0_square(l, phi, th, delta)

    z_bottom = z_bottom_square(delta)
    z_c = z_bottom - r * cos(th)

    z_upper = z_c + r

    mask = (z >= z_bottom) & (z <= z_upper)

    r_sq = np.maximum(0.0, r ** 2 - (z - z_c) ** 2)

    return pi * r_sq * mask / (l * ly)


""" droplet_w2c0 """


def r_droplet_w2c0_square(l, phi, th, delta):
    """
    Радиус сферы для droplet_w2c0.

 
    """
    return (
        3.0 * phi * l
        / (2.0 * pi * (-cos(th)) * (3.0 - cos(th) ** 2))
    ) ** (1.0 / 3.0)


def rho_droplet_w2c0_square(z, l, ly, phi, th, delta):
    """

    """
    z = np.asarray(z, dtype=float)

    r = r_droplet_w2c0_square(l, phi, th, delta)
    th_p = pi - th

    z_bottom = z_bottom_square(delta)

    z_c = z_bottom + r * cos(th_p)
    d = r * cos(th_p)

    L1 = L1_square(r, th)
    L2 = L2_square(r, th)

    z1 = z_bottom + L1
    z2 = z_bottom + L2

    rho = np.zeros_like(z, dtype=float)

    upper = z_c + r
    mask_any = (z >= z_bottom) & (z <= upper)

    if not np.any(mask_any):
        return rho

    z_m = z[mask_any]

    r_sq = np.maximum(0.0, r ** 2 - (z_m - z_c) ** 2)
    r_z = sqrt(r_sq)

    base = pi * r_sq

    in_cut = (z_m > z1) & (z_m < z2)

    cut = np.zeros_like(z_m)
    if np.any(in_cut):
        rz = r_z[in_cut]
        ratio = np.divide(d, rz, out=np.ones_like(rz), where=rz > 0)
        ratio = np.clip(ratio, -1.0, 1.0)
        cut[in_cut] = (
            rz ** 2 * np.arccos(ratio)
            - d * sqrt(np.maximum(0.0, rz ** 2 - d ** 2))
        )

    S = base.copy()
    S[in_cut] = base[in_cut] - cut[in_cut]

    rho[mask_any] = S / (l * ly)

    return rho

""" droplet_w2c1 """


def r_droplet_w2c1_square(l, phi, th, delta):
    """
    Радиус сферы для droplet_w2c1.

    Из tex:
        V_droplet = (2 R^3 / 3) vol_factor(theta) = phi L
        R = ( 3 phi L / (2 vol_factor(theta)) )^(1/3)

        vol_factor(theta) = arccos(cot^2 theta)
                          - cos(theta) (2 + sin^2 theta) arccos(cot theta)
                          + cos^2 theta sqrt(-cos(2 theta))
    """
    vol_factor = vol_factor_droplet_w2c1_square(th)
    return (3.0 * phi * l / (2.0 * vol_factor)) ** (1.0 / 3.0)


def vol_factor_droplet_w2c1_square(th):
    """
    Из tex:
        vol_factor(theta) = arccos(cot^2 theta)
                          - cos(theta) (2 + sin^2 theta) arccos(cot theta)
                          + cos^2 theta sqrt(-cos(2 theta))
    """
    cot_t = cos(th) / sin(th)
    return (
        np.arccos(cot_t ** 2)
        - cos(th) * (2.0 + sin(th) ** 2) * np.arccos(cot_t)
        + (cos(th) ** 2)* sqrt(-cos(2.0 * th))
    )


def rho_droplet_w2c1_square(z, l, ly, phi, th, delta):
    """
    Density profile for droplet_w2c1.

    Формула профиля та же, что у droplet_w2c0, отличается только R.
    """
    z = np.asarray(z, dtype=float)

    r = r_droplet_w2c1_square(l, phi, th, delta)
    th_p = pi - th

    z_bottom = z_bottom_square(delta)

    z_c = z_bottom + r * cos(th_p)
    d = r * cos(th_p)

    L1 = L1_square(r, th)
    L2 = L2_square(r, th)

    z1 = z_bottom + L1
    z2 = z_bottom + L2

    rho = np.zeros_like(z, dtype=float)

    upper = z_c + r
    mask_any = (z >= z_bottom) & (z <= upper)

    if not np.any(mask_any):
        return rho

    z_m = z[mask_any]

    r_sq = np.maximum(0.0, r ** 2 - (z_m - z_c) ** 2)
    r_z = sqrt(r_sq)

    base = pi * r_sq

    in_cut = (z_m > z1) & (z_m < z2)

    cut = np.zeros_like(z_m)
    if np.any(in_cut):
        rz = r_z[in_cut]
        ratio = np.divide(d, rz, out=np.ones_like(rz), where=rz > 0)
        ratio = np.clip(ratio, -1.0, 1.0)
        cut[in_cut] = (
            rz ** 2 * np.arccos(ratio)
            - d * sqrt(np.maximum(0.0, rz ** 2 - d ** 2))
        )

    S = base.copy()
    S[in_cut] = base[in_cut] - cut[in_cut]

    rho[mask_any] = S / (l * ly)

    return rho

""" capsule_w4c0 """


def r_capsule_w4c0_square(l, phi, th, delta):
    """
    Радиус кривизны мениска для capsule_w4c0.

    Из tex:
        r = -a' / (2 cos(theta)) = (1 - 2 delta) / (2 cos(theta'))
    """
    a_eff = a_eff_square(delta)
    return a_eff / (2.0 * cos(pi - th))


def l_tilde_capsule_w4c0_square(l, phi, th, delta):
    """
    Длина цилиндрической части capsule_w4c0.

    Из tex:
        l_tilde = ( F L + (4/3) pi r^3 - 2 pi a' r^2 + (pi/6) a'^3 )
                  / ( r^2 (pi - 4 theta' + 2 sin(2 theta')) )
    """
    a_eff = a_eff_square(delta)
    r = r_capsule_w4c0_square(l, phi, th, delta)
    th_p = pi - th

    numerator = (
        phi * l
        + (4.0 / 3.0) * pi * (r ** 3)
        - 2.0 * pi * a_eff *( r ** 2)
        + (pi / 6.0) * a_eff ** 3
    )
    denominator = r ** 2 * (pi - 4.0 * th_p + 2.0 * sin(2.0 * th_p))

    return numerator / denominator


def S_cross_capsule_w4c0_square(l, phi, th, delta):
    """
    Площадь поперечного сечения цилиндрической части capsule_w4c0.

    Из tex:
        S_cross = r^2 (pi - 4 theta' + 2 sin(2 theta'))
                = r^2 (4 theta - 3 pi - 2 sin(2 theta))
    """
    r = r_capsule_w4c0_square(l, phi, th, delta)
    th_p = pi - th
    return r ** 2 * (pi - 4.0 * th_p + 2.0 * sin(2.0 * th_p))


def rho_capsule_w4c0_square(z, l, ly, phi, th, delta):
    """
    Density profile for capsule_w4c0.

    """
    z = np.asarray(z, dtype=float)

 
    r = r_capsule_w4c0_square(l, phi, th, delta)
    l_tilde = l_tilde_capsule_w4c0_square(l, phi, th, delta)

    d = 0.5 - delta

    z_bottom = z_bottom_square(delta)
    z_top = z_top_square(delta)

    rho = np.zeros_like(z, dtype=float)

    mask = (z >= z_bottom) & (z <= z_top)

    if not np.any(mask):
        return rho

    z_m = z[mask]

    r_sq = np.maximum(0.0, r ** 2 - z_m ** 2)
    r_z = sqrt(r_sq)

    S = np.zeros_like(z_m)

    inside = r_z <= d
    clipped = ~inside

    S[inside] = 2.0 * l_tilde * r_z[inside] + pi * r_sq[inside]

    if np.any(clipped):
        rz = r_z[clipped]
        ratio = np.divide(d, rz, out=np.ones_like(rz), where=rz > 0)
        ratio = np.clip(ratio, -1.0, 1.0)
        S_seg = rz ** 2 * np.arccos(ratio) - d * sqrt(np.maximum(0.0, rz ** 2 - d ** 2))
        S[clipped] = 2.0 * l_tilde * d + pi * rz ** 2 - 2.0 * S_seg

    rho[mask] = S / (l * ly)

    return rho


""" capsule_w4c4 """


def r_capsule_w4c4_square(l, phi, th, delta):
    """
    Радиус кривизны мениска для capsule_w4c4.

    Из tex:
        r = -a' / (2 cos(theta)) = (1 - 2 delta) / (-2 cos(theta))
    """
    a_eff = a_eff_square(delta)
    return a_eff / (-2.0 * cos(th))


def Omega_capsule_w4c4_square(th):
    """
    Телесный угол мениска capsule_w4c4.

    Из tex:
        Omega = 4 arccos(-cot^2 theta) - 2 pi
    """
    arg = -cos(th) ** 2 / sin(th) ** 2
    arg = np.clip(arg, -1.0, 1.0)
    return 4.0 * np.arccos(arg) - 2.0 * pi


def h_capsule_w4c4_square(l, phi, th, delta):
    """
    Величина h для capsule_w4c4 (высота среза мениска).

    Из tex (perforation / capsule_w4c4 в square_pore_delta):
        h = sqrt(r^2 - (a')^2 / 2)
    """
    a_eff = a_eff_square(delta)
    r = r_capsule_w4c4_square(l, phi, th, delta)
    return sqrt(np.maximum(0.0, r ** 2 - a_eff ** 2 / 2.0))


def l_tilde_capsule_w4c4_square(l, phi, th, delta):
    """
    Длина цилиндрической перемычки capsule_w4c4.

        l_tilde = F L / a'^2 - (4/3) h - (2 r^3 / (3 a'^2)) Omega

    (в tex конечная формула содержит опечатку: первый член записан как F L,
    но из баланса V = a'^2 (l_tilde + 2h) + 2 V_caps = F L следует F L / a'^2)
    """
    a_eff = a_eff_square(delta)
    r = r_capsule_w4c4_square(l, phi, th, delta)
    Omega = Omega_capsule_w4c4_square(th)
    h = h_capsule_w4c4_square(l, phi, th, delta)

    return (
        phi * l / a_eff ** 2            # <-- добавили / a_eff**2
        - (4.0 / 3.0) * h
        - (2.0 * r ** 3 / (3.0 * a_eff ** 2)) * Omega
    )



def rho_capsule_w4c4_square(z, l, ly, phi, th, delta):
    """
    Density profile for capsule_w4c4.

    Из tex (density_profiles):
        r(z) = sqrt(r^2 - z^2)
        d = 1/2 - delta
        S_seg(z) = r^2(z) arccos(d/r(z)) - d sqrt(r^2(z) - d^2)

        [z_bottom, z_top]:
            2 l_tilde d + pi r^2(z) - 2 S_seg(z)
        иначе 0.
    """
    z = np.asarray(z, dtype=float)

   
    r = r_capsule_w4c4_square(l, phi, th, delta)
    l_tilde = l_tilde_capsule_w4c4_square(l, phi, th, delta)

    d = 0.5 - delta

    z_bottom = z_bottom_square(delta)
    z_top = z_top_square(delta)

    rho = np.zeros_like(z, dtype=float)

    mask = (z >= z_bottom) & (z <= z_top)

    if not np.any(mask):
        return rho

    z_m = z[mask]

    r_sq = np.maximum(0.0, r ** 2 - z_m ** 2)
    r_z = sqrt(r_sq)

    ratio = np.divide(d, r_z, out=np.ones_like(r_z), where=r_z > 0)
    ratio = np.clip(ratio, -1.0, 1.0)

    S_seg = r_z ** 2 * np.arccos(ratio) - d * sqrt(np.maximum(0.0, r_z ** 2 - d ** 2))

    S = 2.0 * l_tilde * d + pi * r_sq - 2.0 * S_seg

    rho[mask] = S / (l * ly)

    return rho
""" perforation """


def r_perforation_square(l, phi, th, delta):
    """
    Радиус кривизны перфорации.

    Из tex:
        r = (1/2 - delta) / cos(theta')
    """
    th_p = pi - th
    return (0.5 - delta) / cos(th_p)


def S_u_perforation_square(th):
    """
    Величина S_u из tex.

    Из tex:
        S_u = 2u + sin(2u) = pi - 2t + sin(2t),  где t = theta' = pi - theta
    """
    t = pi - th
    return pi - 2.0 * t + sin(2.0 * t)


def D_isc_perforation_square(l, phi, th, delta):
    """
    Дискриминант D_isc из tex.

    Из tex:
        D_isc = r^2 ( S_u^2 - 16 c^2 + (16/3) c^4 )
                + 8 phi_gas L c / (pi r)

    где c = cos(t) = cos(theta') = -cos(theta),
        phi_gas = 1 - phi.
    """
    r = r_perforation_square(l, phi, th, delta)
    S_u = S_u_perforation_square(th)
    c = cos(pi - th)
    phi_gas = 1.0 - phi

    return (
        r ** 2 * (S_u ** 2 - 16.0 * c ** 2 + (16.0 / 3.0) * c ** 4)
        + 8.0 * phi_gas * l * c / (pi * r)
    )


def A_plus_perforation_square(l, phi, th, delta):
    """
    Величина A_+ из tex.

    Из tex:
        A_+ = (r S_u + sqrt(D_isc)) / (2 c)
    """
    r = r_perforation_square(l, phi, th, delta)
    S_u = S_u_perforation_square(th)
    c = cos(pi - th)
    D_isc = D_isc_perforation_square(l, phi, th, delta)

    return (r * S_u + sqrt(np.maximum(0.0, D_isc))) / (2.0 * c)


def D0_perforation_square(l, phi, th, delta):
    """
    Диаметр перфорации в центре поры.

    Из tex:
        D(0) = A_+ - 2r
    """
    r = r_perforation_square(l, phi, th, delta)
    A_plus = A_plus_perforation_square(l, phi, th, delta)
    return A_plus - 2.0 * r


def D_perforation_square(z, l, phi, th, delta):
    """
    Диаметр перфорации D(z).

    Из tex:
        D(z) = D(0) + 2r - 2 sqrt(r^2 - z^2)

    Эквивалентно D(beta) = A - 2r cos(beta), A = D(0) + 2r.
    """
    r = r_perforation_square(l, phi, th, delta)
    D0 = D0_perforation_square(l, phi, th, delta)
    return D0 + 2.0 * r - 2.0 * sqrt(np.maximum(0.0, r ** 2 - z ** 2))


def rho_perforation_square(z, l, ly, phi, th, delta):
    """
    Density profile for perforation.

    Из tex (density_profiles):
        rho(z) = ((1 - 2 delta) L - pi D^2(z) / 4) / (L Ly),
                 z in [z_bottom, z_top]
        иначе 0.

        D(z) = D(0) + 2r - 2 sqrt(r^2 - z^2)
    """
    z = np.asarray(z, dtype=float)

    a_eff = a_eff_square(delta)

    z_bottom = z_bottom_square(delta)
    z_top = z_top_square(delta)

    rho = np.zeros_like(z, dtype=float)

    mask = (z >= z_bottom) & (z <= z_top)

    if not np.any(mask):
        return rho

    D_z = D_perforation_square(z[mask], l, phi, th, delta)

    rho[mask] = (a_eff * l - pi * D_z ** 2 / 4.0) / (l * ly)

    return rho

# ============================================================
# Analytical surface areas / capillary functional
# ============================================================
#
# These expressions are transferred from the stability-diagram
# notebook.  The volume fraction phi is always defined relative
# to the ORIGINAL unit-side square pore:
#
#     V = phi * l
#
# The wall layer only changes the accessible square side:
#
#     a' = 1 - 2*delta
#
# flag = False:
#     return the capillary functional
#         S_LV - cos(theta) * S_SL
#
# flag = True:
#     return the purely geometrical sum
#         S_LV + S_SL
#
# All functions accept either scalars or broadcastable NumPy arrays.
# ============================================================


def a_prime_square(delta):
    """Accessible square side after a parallel inward wall offset delta."""
    return 1.0 - 2.0 * delta


def _square_broadcast(l, phi):
    """Broadcast l and phi and remember whether the result was scalar."""
    l_arr, phi_arr = np.broadcast_arrays(
        np.asarray(l, dtype=float),
        np.asarray(phi, dtype=float),
    )
    scalar = l_arr.ndim == 0
    return l_arr, phi_arr, scalar


def _square_restore(value, scalar):
    arr = np.asarray(value)
    return float(arr) if scalar else arr


def _square_inf(l, phi):
    l_arr, phi_arr, scalar = _square_broadcast(l, phi)
    out = np.full(l_arr.shape, np.inf, dtype=float)
    return l_arr, phi_arr, out, scalar


def _square_energy(S_LV, S_SL, th, flag):
    if flag:
        return S_LV + S_SL
    return S_LV - cos(th) * S_SL

# ------------------------------------------------------------
# filament_w0c0
# ------------------------------------------------------------

def cond_filament_w0c0_square(l, phi, th, delta):
    """
    Условие существования центрального цилиндра в квадратной поре.

    Из tex / ipynb:
        phi_max = (pi/4) * a'^2
        cond: 0 < phi < phi_max
    """
    del th
    l_arr, phi_arr, scalar = _square_broadcast(l, phi)
    a_prime = a_prime_square(delta)

    if a_prime <= 0:
        result = np.zeros(l_arr.shape, dtype=bool)
    else:
        phi_max = (pi / 4.0) * (a_prime ** 2)
        result = (phi_arr > 0.0) & (phi_arr < phi_max)

    return bool(result) if scalar else result


def area_filament_w0c0_square(l, phi, th, delta, flag=False):
    """
    Площадь поверхности / капиллярный функционал для filament_w0c0.

    Из tex:
        S_LV = 2 pi r L,   r = sqrt(phi / pi)
        S_SL = 0
    """
    del th, flag
    l_arr, phi_arr, S, scalar = _square_inf(l, phi)
    cond = np.asarray(
        cond_filament_w0c0_square(l_arr, phi_arr, 0.0, delta),
        dtype=bool,
    )

    r = r_filament_w0c0_square(l_arr, phi_arr, 0.0, delta)
    values = 2.0 * pi * r * l_arr
    S = np.where(cond, values, S)

    return _square_restore(S, scalar)


# ------------------------------------------------------------
# filament_w0c0_straight  (filament_w0c0_straight_square)
# ------------------------------------------------------------

def cond_filament_w0c0_straight_square(l, phi, th, delta):
    """
    Условие существования соосного квадратного шнура.

    Из tex / ipynb:
        cond: 0 < phi < a'^2
    """
    del th
    l_arr, phi_arr, scalar = _square_broadcast(l, phi)
    a_prime = a_prime_square(delta)

    if a_prime <= 0:
        result = np.zeros(l_arr.shape, dtype=bool)
    else:
        result = (phi_arr > 0.0) & (phi_arr < a_prime ** 2)

    return bool(result) if scalar else result


def area_filament_w0c0_straight_square(l, phi, th, delta, flag=False):
    """
    Площадь поверхности / капиллярный функционал для filament_w0c0_straight.

    Из tex:
        S_LV = 4 sqrt(phi) L,   b = sqrt(phi)
        S_SL = 0
    """
    del th, flag
    l_arr, phi_arr, S, scalar = _square_inf(l, phi)
    cond = np.asarray(
        cond_filament_w0c0_straight_square(l_arr, phi_arr, 0.0, delta),
        dtype=bool,
    )

    b = b_filament_w0c0_straight_square(l_arr, phi_arr, 0.0, delta)
    values = 4.0 * b * l_arr
    S = np.where(cond, values, S)

    return _square_restore(S, scalar)


# ------------------------------------------------------------
# filament_w1c0
# ------------------------------------------------------------

def cond_filament_w1c0_square(l, phi, th, delta):
    """
    Условие существования filament_w1c0.

    Из tex:
        denom = theta - sin(theta) cos(theta)
        r = sqrt(phi / denom)
        cond: 2 r < a'
    """
    del l
    a_prime = a_prime_square(delta)

    phi_arr = np.asarray(phi, dtype=float)
    scalar = phi_arr.ndim == 0

    if a_prime <= 0:
        result = np.zeros_like(phi_arr, dtype=bool)
        return bool(result) if scalar else result

    denom = th - sin(th) * cos(th)
    if denom <= 0:
        result = np.zeros_like(phi_arr, dtype=bool)
        return bool(result) if scalar else result

    r = sqrt(phi_arr / denom)
    cond_wall = 2.0 * r < a_prime

    result = (phi_arr > 0.0) & cond_wall
    return bool(result) if scalar else result


def area_filament_w1c0_square(l, phi, th, delta, flag=False):
    """
    Площадь поверхности / капиллярный функционал для filament_w1c0.

    Из tex:
        P_LV = 2 r theta
        P_SL = 2 r sin(theta)
        S_LV = P_LV L
        S_SL = P_SL L
    """
    l_arr, phi_arr, S, scalar = _square_inf(l, phi)
    cond = np.asarray(
        cond_filament_w1c0_square(l_arr, phi_arr, th, delta),
        dtype=bool,
    )

    r = r_filament_w1c0_square(l_arr, phi_arr, th, delta)

    P_LV = 2.0 * r * th
    P_SL = 2.0 * r * sin(th)

    values = l_arr * _square_energy(P_LV, P_SL, th, flag)
    S = np.where(cond, values, S)

    return _square_restore(S, scalar)

# ------------------------------------------------------------
# filament_w4c0
# ------------------------------------------------------------

def cond_filament_w4c0_square(l, phi, th, delta):
    """
    Условие существования filament_w4c0.

    Из tex:
        cond1: theta > 3 pi / 4
        cond2: L1 < 0.5 - delta
    """
    del l
    a_prime = a_prime_square(delta)

    phi_arr = np.asarray(phi, dtype=float)
    scalar = phi_arr.ndim == 0

    if a_prime <= 0 or th <= 3.0 * pi / 4.0:
        result = np.zeros_like(phi_arr, dtype=bool)
        return bool(result) if scalar else result

    num = a_prime ** 2 - phi_arr
    valid_num = num >= 0
    num_safe = np.maximum(num, 0.0)

    r = r_filament_w4c0_square(l, num_safe, th, delta)
    L1 = L1_square(r, th)

    result = (phi_arr > 0.0) & valid_num & (L1 < 0.5 - delta)
    return bool(result) if scalar else result


def area_filament_w4c0_square(l, phi, th, delta, flag=False):
    """
    Площадь поверхности / капиллярный функционал для filament_w4c0.

    Из tex:
        S_SL = 4 L (1 - 2 delta - 2 L1)
        S_LV = 4 alpha r L
    """
    l_arr, phi_arr, S, scalar = _square_inf(l, phi)
    cond = np.asarray(
        cond_filament_w4c0_square(l_arr, phi_arr, th, delta),
        dtype=bool,
    )

    if th <= 3.0 * pi / 4.0:
        return _square_restore(S, scalar)

    a_prime = a_prime_square(delta)
    alpha = np.pi/2 - 2*(np.pi - th)

    phi_safe = np.maximum(a_prime ** 2 - phi_arr, 0.0)
    r = r_filament_w4c0_square(l_arr, phi_safe, th, delta)
    L1 = L1_square(r, th)

    S_SL = 4.0 * l_arr * (a_prime - 2.0 * L1)
    S_LV = 4.0 * alpha * r * l_arr

    values = _square_energy(S_LV, S_SL, th, flag)
    S = np.where(cond, values, S)

    return _square_restore(S, scalar)


# ------------------------------------------------------------
# filament_w4c1
# ------------------------------------------------------------

def cond_filament_w4c1_square(l, phi, th, delta):
    """
    Условие существования filament_w4c1.

    Из tex:
        cond1: theta > 3 pi / 4
        cond2: L1 < 0.5 - delta
    """
    del l
    a_prime = a_prime_square(delta)

    phi_arr = np.asarray(phi, dtype=float)
    scalar = phi_arr.ndim == 0

    if a_prime <= 0 or th <= 3.0 * pi / 4.0:
        result = np.zeros_like(phi_arr, dtype=bool)
        return bool(result) if scalar else result

    num = a_prime ** 2 - phi_arr
    valid_num = num >= 0
    num_safe = np.maximum(num, 0.0)

    r = r_filament_w4c1_square(l, num_safe, th, delta)
    L1 = L1_square(r, th)

    result = (phi_arr > 0.0) & valid_num & (L1 < 0.5 - delta)
    return bool(result) if scalar else result


def area_filament_w4c1_square(l, phi, th, delta, flag=False):
    """
    Площадь поверхности / капиллярный функционал для filament_w4c1.

    Из tex:
        S_SL = 2 L (1 - 2 delta - 2 L1) + 2 L (1 - 2 delta - L1)
        S_LV = 3 alpha r L
    """
    l_arr, phi_arr, S, scalar = _square_inf(l, phi)
    cond = np.asarray(
        cond_filament_w4c1_square(l_arr, phi_arr, th, delta),
        dtype=bool,
    )

    if th <= 3.0 * pi / 4.0:
        return _square_restore(S, scalar)

    a_prime = a_prime_square(delta)
    alpha = np.pi/2 - 2*(np.pi - th)

    phi_safe = np.maximum(a_prime ** 2 - phi_arr, 0.0)
    r = r_filament_w4c1_square(l_arr, phi_safe, th, delta)
    L1 = L1_square(r, th)

    S_SL = (
        2.0 * l_arr * (a_prime - 2.0 * L1)
        + 2.0 * l_arr * (a_prime - L1)
    )
    S_LV = 3.0 * alpha * r * l_arr

    values = _square_energy(S_LV, S_SL, th, flag)
    S = np.where(cond, values, S)

    return _square_restore(S, scalar)


# ------------------------------------------------------------
# filament_w4c2
# ------------------------------------------------------------

def cond_filament_w4c2_square(l, phi, th, delta):
    """
    Условие существования filament_w4c2.

    Из tex:
        cond1: theta > 3 pi / 4
        cond2: L1 < 0.5 - delta
    """
    del l
    a_prime = a_prime_square(delta)

    phi_arr = np.asarray(phi, dtype=float)
    scalar = phi_arr.ndim == 0

    if a_prime <= 0 or th <= 3.0 * pi / 4.0:
        result = np.zeros_like(phi_arr, dtype=bool)
        return bool(result) if scalar else result

    num = a_prime ** 2 - phi_arr
    valid_num = num >= 0
    num_safe = np.maximum(num, 0.0)

    r = r_filament_w4c2_square(l, num_safe, th, delta)
    L1 = L1_square(r, th)

    result = (phi_arr > 0.0) & valid_num & (L1 < 0.5 - delta)
    return bool(result) if scalar else result


def area_filament_w4c2_square(l, phi, th, delta, flag=False):
    """
    Площадь поверхности / капиллярный функционал для filament_w4c2.

    Из tex:
        S_SL = 4 L (1 - 2 delta - L1)
        S_LV = 2 alpha r L
    """
    l_arr, phi_arr, S, scalar = _square_inf(l, phi)
    cond = np.asarray(
        cond_filament_w4c2_square(l_arr, phi_arr, th, delta),
        dtype=bool,
    )

    if th <= 3.0 * pi / 4.0:
        return _square_restore(S, scalar)

    a_prime = a_prime_square(delta)
    alpha = np.pi/2 - 2*(np.pi - th)

    phi_safe = np.maximum(a_prime ** 2 - phi_arr, 0.0)
    r = r_filament_w4c2_square(l_arr, phi_safe, th, delta)
    L1 = L1_square(r, th)

    S_SL = 4.0 * l_arr * (a_prime - L1)
    S_LV = 2.0 * alpha * r * l_arr

    values = _square_energy(S_LV, S_SL, th, flag)
    S = np.where(cond, values, S)

    return _square_restore(S, scalar)
# ------------------------------------------------------------
# filament_w4c3
# ------------------------------------------------------------

def cond_filament_w4c3_square(l, phi, th, delta):
    """
    Условие существования filament_w4c3.

    Из tex:
        cond1: theta > 3 pi / 4
        cond2: L1 < 1 - 2 delta
    """
    del l
    a_prime = a_prime_square(delta)

    phi_arr = np.asarray(phi, dtype=float)
    scalar = phi_arr.ndim == 0

    if a_prime <= 0 or th <= 3.0 * pi / 4.0:
        result = np.zeros_like(phi_arr, dtype=bool)
        return bool(result) if scalar else result

    num = a_prime ** 2 - phi_arr
    valid_num = num >= 0
    num_safe = np.maximum(num, 0.0)

    r = r_filament_w4c3_square(l, num_safe, th, delta)
    L1 = L1_square(r, th)

    result = (phi_arr > 0.0) & valid_num & (L1 < a_prime)
    return bool(result) if scalar else result


def area_filament_w4c3_square(l, phi, th, delta, flag=False):
    """
    Площадь поверхности / капиллярный функционал для filament_w4c3.

    Из tex:
        S_SL = 2 L (1 - 2 delta - L1) + 2 L (1 - 2 delta)
        S_LV = alpha r L
    """
    l_arr, phi_arr, S, scalar = _square_inf(l, phi)
    cond = np.asarray(
        cond_filament_w4c3_square(l_arr, phi_arr, th, delta),
        dtype=bool,
    )

    if th <= 3.0 * pi / 4.0:
        return _square_restore(S, scalar)

    a_prime = a_prime_square(delta)
    alpha = np.pi/2 - 2*(np.pi - th)

    phi_safe = np.maximum(a_prime ** 2 - phi_arr, 0.0)
    r = r_filament_w4c3_square(l_arr, phi_safe, th, delta)
    L1 = L1_square(r, th)

    S_SL = (
        2.0 * l_arr * (a_prime - L1)
        + 2.0 * l_arr * a_prime
    )
    S_LV = alpha * r * l_arr

    values = _square_energy(S_LV, S_SL, th, flag)
    S = np.where(cond, values, S)

    return _square_restore(S, scalar)
# ------------------------------------------------------------
# filament_w3c2
# ------------------------------------------------------------

def cond_filament_w3c2_square(l, phi, th, delta):
    del l
    a_prime = a_prime_square(delta)

    phi_arr = np.asarray(phi, dtype=float)
    scalar = phi_arr.ndim == 0

    if a_prime <= 0 or th <= pi / 2.0:
        result = np.zeros_like(phi_arr, dtype=bool)
        return bool(result) if scalar else result

    r = r_filament_w3c2_square(l, phi_arr, th, delta)
    alpha = alpha_square(th)

    h = (phi_arr - 0.5 * r ** 2 * (alpha - sin(alpha))) / a_prime

    cond_h = h > 0
    cond_top = (r * (1.0 - cos(alpha / 2.0)) + h) < a_prime

    result = (phi_arr > 0.0) & cond_h & cond_top
    return bool(result) if scalar else result

def area_filament_w3c2_square(l, phi, th, delta, flag=False):
    """
    Площадь поверхности / капиллярный функционал для filament_w3c2.

    Из tex:
        S_SL = L (2h + (1 - 2 delta))
        S_LV = L r alpha
    """
    l_arr, phi_arr, S, scalar = _square_inf(l, phi)
    cond = np.asarray(
        cond_filament_w3c2_square(l_arr, phi_arr, th, delta),
        dtype=bool,
    )

    a_prime = a_prime_square(delta)

    r = r_filament_w3c2_square(l_arr, phi_arr, th, delta)
    alpha = alpha_square(th)

    h = (phi_arr - 0.5 * r ** 2 * (alpha - sin(alpha))) / a_prime

    S_SL = l_arr * (2.0 * h + a_prime)
    S_LV = l_arr * r * alpha

    values = _square_energy(S_LV, S_SL, th, flag)
    S = np.where(cond, values, S)

    return _square_restore(S, scalar)


# ------------------------------------------------------------
# filament_w2c0_opposite
# ------------------------------------------------------------
def cond_filament_w2c0_opposite_square(l, phi, th, delta):
    """
    Условие существования filament_w2c0_opposite.

    Из tex:
        cond1: h > 0
        cond2: 2r(1 - cos(alpha/2)) + h < 1 + 2 delta
    """
    del l
    a_prime = a_prime_square(delta)

    phi_arr = np.asarray(phi, dtype=float)
    scalar = phi_arr.ndim == 0

    if a_prime <= 0 or th <= pi / 2.0 or th >= pi:
        result = np.zeros_like(phi_arr, dtype=bool)
        return bool(result) if scalar else result

    r = r_filament_w2c0_opposite_square(l, phi_arr, th, delta)
    alpha = alpha_square(th)

    h = (phi_arr - r ** 2 * (alpha - sin(alpha))) / a_prime

    cond_h = h > 0
    cond_top = (2.0 * r * (1.0 - cos(alpha / 2.0)) + h) < a_prime

    result = (phi_arr > 0.0) & cond_h & cond_top
    return bool(result) if scalar else result

def area_filament_w2c0_opposite_square(l, phi, th, delta, flag=False):
    """
    Площадь поверхности / капиллярный функционал для filament_w2c0_opposite.

    Из tex:
        S_SL = L 2h
        S_LV = 2 L r alpha
    """
    l_arr, phi_arr, S, scalar = _square_inf(l, phi)
    cond = np.asarray(
        cond_filament_w2c0_opposite_square(l_arr, phi_arr, th, delta),
        dtype=bool,
    )

    a_prime = a_prime_square(delta)

    r = r_filament_w2c0_opposite_square(l_arr, phi_arr, th, delta)
    alpha = alpha_square(th)

    h = (phi_arr - r ** 2 * (alpha - sin(alpha))) / a_prime

    S_SL = l_arr * 2.0 * h
    S_LV = 2.0 * l_arr * r * alpha

    values = _square_energy(S_LV, S_SL, th, flag)
    S = np.where(cond, values, S)

    return _square_restore(S, scalar)
# ------------------------------------------------------------
# filament_w2c1
# ------------------------------------------------------------

def cond_filament_w2c1_square(l, phi, th, delta):
    """
    Условие существования filament_w2c1.

    Из tex:
        cond1: L_2 <= 1 - 2 delta
        cond2: r (1 - cos(theta)) < 1 - 2 delta
    """
    del l
    a_prime = a_prime_square(delta)

    phi_arr = np.asarray(phi, dtype=float)
    scalar = phi_arr.ndim == 0

    if a_prime <= 0:
        result = np.zeros_like(phi_arr, dtype=bool)
        return bool(result) if scalar else result

    denom = th - pi / 4.0 + cos(th) ** 2 - sin(th) * cos(th)
    if denom <= 0:
        result = np.zeros_like(phi_arr, dtype=bool)
        return bool(result) if scalar else result

    r = r_filament_w2c1_square(l, phi_arr, th, delta)

    cond_contact = L2_square(r, th) <= a_prime
    cond_opp = r * (1.0 - cos(th)) < a_prime

    result = (phi_arr > 0.0) & cond_contact & cond_opp
    return bool(result) if scalar else result

def area_filament_w2c1_square(l, phi, th, delta, flag=False):
    """
    Площадь поверхности / капиллярный функционал для filament_w2c1.

    Из tex:
        P_SL = 2 L_2
        P_LV = r (2 theta - pi/2)
        S_SL = P_SL L
        S_LV = P_LV L
    """
    l_arr, phi_arr, S, scalar = _square_inf(l, phi)
    cond = np.asarray(
        cond_filament_w2c1_square(l_arr, phi_arr, th, delta),
        dtype=bool,
    )

    denom = th - pi / 4.0 + cos(th) ** 2 - sin(th) * cos(th)
    if denom <= 0:
        return _square_restore(S, scalar)

    r = r_filament_w2c1_square(l_arr, phi_arr, th, delta)

    P_SL = 2.0 * L2_square(r, th)
    P_LV = r * (2.0 * th - pi / 2.0)

    S_SL = P_SL * l_arr
    S_LV = P_LV * l_arr

    values = _square_energy(S_LV, S_SL, th, flag)
    S = np.where(cond, values, S)

    return _square_restore(S, scalar)
# ------------------------------------------------------------
# droplet_w1c0
# ------------------------------------------------------------

def cond_droplet_w1c0_square(l, phi, th, delta):
    """
    Условие существования droplet_w1c0.

    Из tex:
        cond1: theta >= pi/2
        cond2: 2R < L
        cond3: 2R <= 1 - 2 delta
    """
    l_arr, phi_arr, scalar = _square_broadcast(l, phi)
    a_prime = a_prime_square(delta)

    if a_prime <= 0 or th < pi / 2.0:
        result = np.zeros(l_arr.shape, dtype=bool)
        return bool(result) if scalar else result

    vol_factor = (pi / 3.0) * (2.0 - 3.0 * cos(th) + cos(th) ** 3)
    if vol_factor <= 0:
        result = np.zeros(l_arr.shape, dtype=bool)
        return bool(result) if scalar else result

    r = r_droplet_w1c0_square(l_arr, phi_arr, th, delta)

    cond_length = 2.0 * r < l_arr
    cond_width = 2.0 * r <= a_prime

    result = (phi_arr > 0.0) & cond_length & cond_width
    return bool(result) if scalar else result


def area_droplet_w1c0_square(l, phi, th, delta, flag=False):
    """
    Площадь поверхности / капиллярный функционал для droplet_w1c0.

    Из tex:
        S_LV = 2 pi r^2 (1 - cos(theta))
        S_SL = pi r^2 sin^2(theta)
    """
    l_arr, phi_arr, S, scalar = _square_inf(l, phi)
    cond = np.asarray(
        cond_droplet_w1c0_square(l_arr, phi_arr, th, delta),
        dtype=bool,
    )

    if th < pi / 2.0:
        return _square_restore(S, scalar)

    r = r_droplet_w1c0_square(l_arr, phi_arr, th, delta)

    S_LV = 2.0 * pi * r ** 2 * (1.0 - cos(th))
    S_SL = pi * r ** 2 * sin(th) ** 2

    values = _square_energy(S_LV, S_SL, th, flag)
    S = np.where(cond, values, S)

    return _square_restore(S, scalar)


# ------------------------------------------------------------
# droplet_w2c0
# ------------------------------------------------------------

def cond_droplet_w2c0_square(l, phi, th, delta):
    """
    Условие существования droplet_w2c0.

    Из tex:
        cond1: theta > 3 pi / 4
        cond2: 2R <= L
        cond3: R (1 - cos(theta)) <= 1 - 2 delta
    """
    l_arr, phi_arr, scalar = _square_broadcast(l, phi)
    a_prime = a_prime_square(delta)

    if a_prime <= 0 or th <= 3.0 * pi / 4.0:
        result = np.zeros(l_arr.shape, dtype=bool)
        return bool(result) if scalar else result

    vol_factor = (2.0 * pi / 3.0) * (-cos(th)) * (3.0 - cos(th) ** 2)
    if vol_factor <= 0:
        result = np.zeros(l_arr.shape, dtype=bool)
        return bool(result) if scalar else result

    r = r_droplet_w2c0_square(l_arr, phi_arr, th, delta)

    cond_length = 2.0 * r <= l_arr
    cond_wall = r * (1.0 - cos(th)) <= a_prime

    result = (phi_arr > 0.0) & cond_length & cond_wall
    return bool(result) if scalar else result


def area_droplet_w2c0_square(l, phi, th, delta, flag=False):
    """
    Площадь поверхности / капиллярный функционал для droplet_w2c0.

    Из tex:
        S_LV = -4 pi r^2 cos(theta)
        S_SL = 2 pi r^2 sin^2(theta)
    """
    l_arr, phi_arr, S, scalar = _square_inf(l, phi)
    cond = np.asarray(
        cond_droplet_w2c0_square(l_arr, phi_arr, th, delta),
        dtype=bool,
    )

    if th <= 3.0 * pi / 4.0 or th > pi:
        return _square_restore(S, scalar)

    r = r_droplet_w2c0_square(l_arr, phi_arr, th, delta)

    S_LV = -4.0 * pi * (r ** 2) * cos(th)
    S_SL = 2.0 * pi * (r ** 2) * (sin(th) ** 2)

    values = _square_energy(S_LV, S_SL, th, flag)
    S = np.where(cond, values, S)

    return _square_restore(S, scalar)


# ------------------------------------------------------------
# droplet_w2c1
# ------------------------------------------------------------

def cond_droplet_w2c1_square(l, phi, th, delta):
    """
    Условие существования droplet_w2c1.

    Из tex:
        cond1: pi/2 < theta < 3 pi / 4
        cond2: 2R <= L
        cond3: R (1 - cos(theta)) <= 1 - 2 delta
    """
    l_arr, phi_arr, scalar = _square_broadcast(l, phi)
    a_prime = a_prime_square(delta)

    if a_prime <= 0 or th <= pi / 2.0 or th > 3.0 * pi / 4.0:
        result = np.zeros(l_arr.shape, dtype=bool)
        return bool(result) if scalar else result

    vol_factor = vol_factor_droplet_w2c1_square(th)
    if vol_factor <= 0:
        result = np.zeros(l_arr.shape, dtype=bool)
        return bool(result) if scalar else result

    r = r_droplet_w2c1_square(l_arr, phi_arr, th, delta)

    cond_length = 2.0 * r <= l_arr
    cond_wall = r * (1.0 - cos(th)) <= a_prime

    result = (phi_arr > 0.0) & cond_length & cond_wall
    return bool(result) if scalar else result


def area_droplet_w2c1_square(l, phi, th, delta, flag=False):
    """
    Площадь поверхности / капиллярный функционал для droplet_w2c1.

    Из tex:
        S_LV = 2 r^2 ( arccos(cot^2 theta) - 2 cos(theta) arccos(cot theta) )
        S_SL = 2 r^2 ( sin^2 theta arccos(cot theta)
                       - cos(theta) sqrt(-cos(2 theta)) )
    """
    l_arr, phi_arr, S, scalar = _square_inf(l, phi)
    cond = np.asarray(
        cond_droplet_w2c1_square(l_arr, phi_arr, th, delta),
        dtype=bool,
    )

    if th <= pi / 2.0 or th >= 3.0 * pi / 4.0:
        return _square_restore(S, scalar)

    vol_factor = vol_factor_droplet_w2c1_square(th)
    if vol_factor <= 0:
        return _square_restore(S, scalar)

    r = r_droplet_w2c1_square(l_arr, phi_arr, th, delta)

    cot_t = cos(th) / sin(th)

    S_LV = 2.0 * r ** 2 * (
        np.arccos(cot_t ** 2)
        - 2.0 * cos(th) * np.arccos(cot_t)
    )
    S_SL = 2.0 * r ** 2 * (
        sin(th) ** 2 * np.arccos(cot_t)
        - cos(th) * sqrt(-cos(2.0 * th))
    )

    values = _square_energy(S_LV, S_SL, th, flag)
    S = np.where(cond, values, S)

    return _square_restore(S, scalar)

# ------------------------------------------------------------
# capsule_w4c0
# ------------------------------------------------------------

def cond_capsule_w4c0_square(l, phi, th, delta):
    """
    Условие существования capsule_w4c0.

    Из tex:
        cond1: 3 pi / 4 <= theta <= pi
        cond2: l_tilde >= 0
        cond3: 2r + l_tilde <= L
    """
    l_arr, phi_arr, scalar = _square_broadcast(l, phi)
    a_prime = a_prime_square(delta)

    if a_prime <= 0 or th < 3.0 * pi / 4.0 or th > pi:
        result = np.zeros(l_arr.shape, dtype=bool)
        return bool(result) if scalar else result

    S_cross = S_cross_capsule_w4c0_square(l_arr, phi_arr, th, delta)
    if np.any(S_cross <= 0):
        result = np.zeros(l_arr.shape, dtype=bool)
        return bool(result) if scalar else result

    r = r_capsule_w4c0_square(l_arr, phi_arr, th, delta)
    l_tilde = l_tilde_capsule_w4c0_square(l_arr, phi_arr, th, delta)

    cond_l = l_tilde >= 0.0
    cond_len = (2.0 * r + l_tilde) <= l_arr

    result = (phi_arr > 0.0) & cond_l & cond_len
    return bool(result) if scalar else result


def area_capsule_w4c0_square(l, phi, th, delta, flag=False):
    """
    Площадь поверхности / капиллярный функционал для capsule_w4c0.

    Из tex:
        P_LV = 2 r (pi - 4 theta') = 8 r (theta - 3 pi / 4)
        P_SL = 8 r sin(theta') = -8 r sin(theta)

        S_LV = 4 pi r (a' - r) + l_tilde P_LV
        S_SL = 4 pi r^2 - pi a'^2 + l_tilde P_SL
    """
    l_arr, phi_arr, S, scalar = _square_inf(l, phi)
    cond = np.asarray(
        cond_capsule_w4c0_square(l_arr, phi_arr, th, delta),
        dtype=bool,
    )

    if th < 3.0 * pi / 4.0 or th > pi:
        return _square_restore(S, scalar)

    a_prime = a_prime_square(delta)
    th_p = pi - th

    r = r_capsule_w4c0_square(l_arr, phi_arr, th, delta)
    l_tilde = l_tilde_capsule_w4c0_square(l_arr, phi_arr, th, delta)

    P_LV = 2.0 * r * (pi - 4.0 * th_p)
    P_SL = 8.0 * r * sin(th_p)

    S_LV = 4.0 * pi * r * (a_prime - r) + l_tilde * P_LV
    S_SL = 4.0 * pi * r ** 2 - pi * a_prime ** 2 + l_tilde * P_SL

    values = _square_energy(S_LV, S_SL, th, flag)
    S = np.where(cond, values, S)

    return _square_restore(S, scalar)


# ------------------------------------------------------------
# capsule_w4c4
# ------------------------------------------------------------

def cond_capsule_w4c4_square(l, phi, th, delta):
    """
    Условие существования capsule_w4c4.

    Из tex:
        cond1: pi / 2 <= theta <= 3 pi / 4
        cond2: l_tilde >= 0
        cond3: 2r + l_tilde <= L
    """
    l_arr, phi_arr, scalar = _square_broadcast(l, phi)
    a_prime = a_prime_square(delta)

    if a_prime <= 0 or th < pi / 2.0 or th > 3.0 * pi / 4.0:
        result = np.zeros(l_arr.shape, dtype=bool)
        return bool(result) if scalar else result

    r = r_capsule_w4c4_square(l_arr, phi_arr, th, delta)
    l_tilde = l_tilde_capsule_w4c4_square(l_arr, phi_arr, th, delta)

    cond_l = l_tilde >= 0.0
    cond_len = (2.0 * r + l_tilde) <= l_arr

    result = (phi_arr > 0.0) & cond_l & cond_len
    return bool(result) if scalar else result


def area_capsule_w4c4_square(l, phi, th, delta, flag=False):
    """
    Площадь поверхности / капиллярный функционал для capsule_w4c4.

    Из tex (строгий вывод):
        S_LV = 2 r^2 Omega
        S_SL = 4 a' l_tilde + S_spot

        S_spot = 4 S_strip
        S_strip = a' x + 2 rho^2 arcsin(a' / (2 rho))
        rho = sqrt(r^2 - a'^2 / 4)
        x   = sqrt(r^2 - a'^2 / 2)
    """
    l_arr, phi_arr, S, scalar = _square_inf(l, phi)
    cond = np.asarray(
        cond_capsule_w4c4_square(l_arr, phi_arr, th, delta),
        dtype=bool,
    )

    if th < pi / 2.0 or th > 3.0 * pi / 4.0:
        return _square_restore(S, scalar)

    a_prime = a_prime_square(delta)

    r = r_capsule_w4c4_square(l_arr, phi_arr, th, delta)
    l_tilde = l_tilde_capsule_w4c4_square(l_arr, phi_arr, th, delta)
    Omega = Omega_capsule_w4c4_square(th)

    rho_sq = np.maximum(0.0, r ** 2 - a_prime ** 2 / 4.0)
    rho = sqrt(rho_sq)

    x_sq = np.maximum(0.0, r ** 2 - a_prime ** 2 / 2.0)
    x = sqrt(x_sq)

    arg_arcsin = np.divide(
        a_prime / 2.0,
        rho,
        out=np.zeros_like(rho),
        where=rho > 0,
    )
    arg_arcsin = np.clip(arg_arcsin, -1.0, 1.0)

    S_strip = a_prime * x + 2.0 * rho_sq * np.arcsin(arg_arcsin)
    S_spot = 4.0 * S_strip

    S_LV = 2.0 * r ** 2 * Omega
    S_SL = 4.0 * a_prime * l_tilde + S_spot

    values = _square_energy(S_LV, S_SL, th, flag)
    S = np.where(cond, values, S)

    return _square_restore(S, scalar)

# ------------------------------------------------------------
# perforation
# ------------------------------------------------------------

def cond_perforation_square(l, phi, th, delta):
    """
    Условие существования perforation.

    Из tex:
        cond1: D_isc >= 0
        cond2: D(0) > 0
        cond3: A_+ - 2r cos(u) < 1 - 2 delta
        cond4: A_+ - 2r cos(u) < L

    phi — относительный объём жидкости F = V_liq / V_pore.
    phi_gas = 1 - F.
    """
    l_arr, phi_arr, scalar = _square_broadcast(l, phi)
    a_prime = a_prime_square(delta)

    if a_prime <= 0 or th < pi / 2.0 or th > pi:
        result = np.zeros(l_arr.shape, dtype=bool)
        return bool(result) if scalar else result

    c = cos(pi - th)
    if abs(c) < 1e-12:
        result = np.zeros(l_arr.shape, dtype=bool)
        return bool(result) if scalar else result

    r = r_perforation_square(l_arr, phi_arr, th, delta)
    S_u = S_u_perforation_square(th)
    phi_gas = 1.0 - phi_arr

    D_isc = (
        r ** 2 * (S_u ** 2 - 16.0 * c ** 2 + (16.0 / 3.0) * c ** 4)
        + 8.0 * phi_gas * l_arr * c / (pi * r)
    )

    cond_disc = D_isc >= 0.0
    sqrt_D = sqrt(np.maximum(D_isc, 0.0))

    A_plus = (r * S_u + sqrt_D) / (2.0 * c)
    D0 = A_plus - 2.0 * r

    cond_D0 = D0 > 0.0

    cos_u = sin(th)
    D_u = A_plus - 2.0 * r * cos_u
    cond_Du_pore = D_u < a_prime
    cond_Du_L = D_u < l_arr

    result = (
        (phi_arr > 0.0)
        & (phi_arr < 1.0)
        & cond_disc
        & cond_D0
        & cond_Du_pore
        & cond_Du_L
    )
    return bool(result) if scalar else result


def area_perforation_square(l, phi, th, delta, flag=False):
    """
    Площадь поверхности / капиллярный функционал для perforation.

    Из tex:
        S_LG = pi sqrt(D_isc)
        S_SL = 4 (1 - 2 delta) L - pi (A_+ - 2r cos(u))^2 / 2
    """
    l_arr, phi_arr, S, scalar = _square_inf(l, phi)
    cond = np.asarray(
        cond_perforation_square(l_arr, phi_arr, th, delta),
        dtype=bool,
    )

    if th < pi / 2.0 or th > pi:
        return _square_restore(S, scalar)

    a_prime = a_prime_square(delta)
    c = cos(pi - th)
    if abs(c) < 1e-12:
        return _square_restore(S, scalar)

    r = r_perforation_square(l_arr, phi_arr, th, delta)
    S_u = S_u_perforation_square(th)
    phi_gas = 1.0 - phi_arr

    D_isc = (
        r ** 2 * (S_u ** 2 - 16.0 * c ** 2 + (16.0 / 3.0) * c ** 4)
        + 8.0 * phi_gas * l_arr * c / (pi * r)
    )
    sqrt_D = sqrt(np.maximum(D_isc, 0.0))

    A_plus = (r * S_u + sqrt_D) / (2.0 * c)

    cos_u = sin(th)
    D_u = A_plus - 2.0 * r * cos_u

    S_LG = pi * sqrt_D
    S_SL = 4.0 * a_prime * l_arr - 0.5 * pi * D_u ** 2

    values = _square_energy(S_LG, S_SL, th, flag)
    S = np.where(cond, values, S)

    return _square_restore(S, scalar)

# ------------------------------------------------------------
# PANDA-style S_* aliases
# ------------------------------------------------------------

S_filament_w0c0_square = area_filament_w0c0_square
S_filament_w0c0_straight_square = area_filament_w0c0_straight_square
S_filament_w1c0_square = area_filament_w1c0_square
S_filament_w4c0_square = area_filament_w4c0_square
S_filament_w4c1_square = area_filament_w4c1_square
S_filament_w4c2_square = area_filament_w4c2_square
S_filament_w4c3_square = area_filament_w4c3_square
S_filament_w3c2_square = area_filament_w3c2_square
S_filament_w2c0_opposite_square = area_filament_w2c0_opposite_square
S_filament_w2c1_square = area_filament_w2c1_square
S_droplet_w1c0_square = area_droplet_w1c0_square
S_droplet_w2c0_square = area_droplet_w2c0_square
S_droplet_w2c1_square = area_droplet_w2c1_square
S_capsule_w4c0_square = area_capsule_w4c0_square
S_capsule_w4c4_square = area_capsule_w4c4_square
S_perforation_square = area_perforation_square


# ------------------------------------------------------------
# Registries
# ------------------------------------------------------------

AREA_SQUARE = {
    "filament_w0c0": area_filament_w0c0_square,
    "filament_w0c0_straight": area_filament_w0c0_straight_square,
    "filament_w1c0": area_filament_w1c0_square,
    "filament_w4c0": area_filament_w4c0_square,
    "filament_w4c1": area_filament_w4c1_square,
    "filament_w4c2": area_filament_w4c2_square,
    "filament_w4c3": area_filament_w4c3_square,
    "filament_w3c2": area_filament_w3c2_square,
    "filament_w2c0_opposite": area_filament_w2c0_opposite_square,
    "filament_w2c1": area_filament_w2c1_square,
    "droplet_w1c0": area_droplet_w1c0_square,
    "droplet_w2c0": area_droplet_w2c0_square,
    "droplet_w2c1": area_droplet_w2c1_square,
    "capsule_w4c0": area_capsule_w4c0_square,
    "capsule_w4c4": area_capsule_w4c4_square,
    "perforation": area_perforation_square,
}

COND_SQUARE = {
    "filament_w0c0": cond_filament_w0c0_square,
    "filament_w0c0_straight": cond_filament_w0c0_straight_square,
    "filament_w1c0": cond_filament_w1c0_square,
    "filament_w4c0": cond_filament_w4c0_square,
    "filament_w4c1": cond_filament_w4c1_square,
    "filament_w4c2": cond_filament_w4c2_square,
    "filament_w4c3": cond_filament_w4c3_square,
    "filament_w3c2": cond_filament_w3c2_square,
    "filament_w2c0_opposite": cond_filament_w2c0_opposite_square,
    "filament_w2c1": cond_filament_w2c1_square,
    "droplet_w1c0": cond_droplet_w1c0_square,
    "droplet_w2c0": cond_droplet_w2c0_square,
    "droplet_w2c1": cond_droplet_w2c1_square,
    "capsule_w4c0": cond_capsule_w4c0_square,
    "capsule_w4c4": cond_capsule_w4c4_square,
    "perforation": cond_perforation_square,
}

RHO_SQUARE = {
    "filament_w0c0": rho_filament_w0c0_square,
    "filament_w0c0_straight": rho_filament_w0c0_straight_square,
    "filament_w1c0": rho_filament_w1c0_square,
    "filament_w4c0": rho_filament_w4c0_square,
    "filament_w4c1": rho_filament_w4c1_square,
    "filament_w4c2": rho_filament_w4c2_square,
    "filament_w4c3": rho_filament_w4c3_square,
    "filament_w3c2": rho_filament_w3c2_square,
    "filament_w2c0_opposite": rho_filament_w2c0_opposite_square,
    "filament_w2c1": rho_filament_w2c1_square,
    "droplet_w1c0": rho_droplet_w1c0_square,
    "droplet_w2c0": rho_droplet_w2c0_square,
    "droplet_w2c1": rho_droplet_w2c1_square,
    "capsule_w4c0": rho_capsule_w4c0_square,
    "capsule_w4c4": rho_capsule_w4c4_square,
    "perforation": rho_perforation_square,
}


# ------------------------------------------------------------
# get-functions
# ------------------------------------------------------------

def get_area_square(name):
    return AREA_SQUARE[name]


def get_cond_square(name):
    return COND_SQUARE[name]


def get_rho_square(name):
    return RHO_SQUARE[name]


# ------------------------------------------------------------
# __all__
# ------------------------------------------------------------

__all__ = [
    # common geometry
    "BIG_NUM",
    "a_eff_square",
    "a_prime_square",
    "z_bottom_square",
    "z_top_square",

    # auxiliary quantities
    "A_square",
    "alpha_square",
    "L1_square",
    "L2_square",
    "z_center_low_square",
    "z_center_up_square",
    "z1_square",
    "z2_square",
    "z3_square",

    # density profiles
    "r_filament_w0c0_square",
    "rho_filament_w0c0_square",
    "b_filament_w0c0_straight_square",
    "rho_filament_w0c0_straight_square",
    "r_filament_w1c0_square",
    "rho_filament_w1c0_square",
    "r_filament_w4c0_square",
    "rho_filament_w4c0_square",
    "r_filament_w4c1_square",
    "rho_filament_w4c1_square",
    "r_filament_w4c2_square",
    "rho_filament_w4c2_square",
    "r_filament_w4c3_square",
    "rho_filament_w4c3_square",
    "r_filament_w3c2_square",
    "h_filament_w3c2_square",
    "rho_filament_w3c2_square",
    "r_filament_w2c0_opposite_square",
    "h_filament_w2c0_opposite_square",
    "rho_filament_w2c0_opposite_square",
    "r_filament_w2c1_square",
    "rho_filament_w2c1_square",
    "R_droplet_w1c0_square",
    "rho_droplet_w1c0_square",
    "R_droplet_w2c0_square",
    "rho_droplet_w2c0_square",
    "vol_factor_droplet_w2c1_square",
    "R_droplet_w2c1_square",
    "rho_droplet_w2c1_square",
    "r_capsule_w4c0_square",
    "S_cross_capsule_w4c0_square",
    "V_caps_capsule_w4c0_square",
    "l_tilde_capsule_w4c0_square",
    "rho_capsule_w4c0_square",
    "r_capsule_w4c4_square",
    "Omega_capsule_w4c4_square",
    "h_capsule_w4c4_square",
    "l_tilde_capsule_w4c4_square",
    "rho_capsule_w4c4_square",
    "r_perforation_square",
    "S_u_perforation_square",
    "D_isc_perforation_square",
    "A_plus_perforation_square",
    "D0_perforation_square",
    "D_perforation_square",
    "rho_perforation_square",

    # area/cond registries
    "AREA_SQUARE",
    "COND_SQUARE",
    "RHO_SQUARE",
    "get_area_square",
    "get_cond_square",
    "get_rho_square",
]

__all__ += [
    name
    for name in globals()
    if name.startswith(("area_", "S_", "cond_"))
    and name.endswith("_square")
]