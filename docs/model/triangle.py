import numpy as np
from numpy import pi, cos, sin, sqrt, abs

BIG_NUM = np.inf


# ============================================================
# Common geometry of the equilateral triangular pore
#
# Dimensionless convention:
#     original triangle side a = 1
#
# z = 0 is the centroid / intersection of angle bisectors.
#
# Effective wall-layer correction:
#     a_eff = 1 - 2 * sqrt(3) * delta
#     H_eff = sqrt(3) * a_eff / 2
#
# Hence
#     z_bottom = -H_eff / 3
#     z_top    =  2 H_eff / 3
# ============================================================


def a_eff_triangle(delta):
    """Effective side a' = 1 - 2*sqrt(3)*delta."""
    return 1.0 - 2.0 * sqrt(3.0) * delta


def H_eff_triangle(delta):
    return 0.5 * sqrt(3.0) * a_eff_triangle(delta)


def z_bottom_triangle(delta):
    return -H_eff_triangle(delta) / 3.0


def z_top_triangle(delta):
    return 2.0 * H_eff_triangle(delta) / 3.0


""" filament_w0c0 """


def r_filament_w0c0_triangle(l, phi, th, delta):
    """
    Radius of the axial cylindrical phase.

    pi r^2 l = phi * (sqrt(3)/4) * l
    """
    return sqrt(phi * sqrt(3.0) / (4.0 * pi))


def rho_filament_w0c0_triangle(z, l, ly, phi, th, delta):
    """
    Density profile for filament_w0c0.

    The cylinder axis coincides with the pore axis, therefore its
    centre is z = 0.
    """
    z = np.asarray(z, dtype=float)

    r = r_filament_w0c0_triangle(l, phi, th, delta)

    mask = abs(z) <= r
    y = sqrt(np.maximum(0.0, r**2 - z**2))

    return 2.0 * y * mask / ly


""" droplet_w1c0 """


def r_droplet_w1c0_triangle(l, phi, th, delta):
    """
    Radius of the spherical cap.

    V = pi R^3 / 3 * (2 - 3 cos(theta) + cos(theta)^3)
      = phi * l * sqrt(3) / 4
    """
    vol_factor = 2.0 - 3.0 * cos(th) + cos(th)**3

    return (
        3.0 * sqrt(3.0) * phi * l
        / (4.0 * pi * vol_factor)
    ) ** (1.0 / 3.0)


def rho_droplet_w1c0_triangle(z, l, ly, phi, th, delta):
    """
    Density profile for droplet_w1c0.

    The droplet lies on the lower effective wall.
    """
    z = np.asarray(z, dtype=float)

    R = r_droplet_w1c0_triangle(l, phi, th, delta)

    z_bottom = z_bottom_triangle(delta)

    # Centre of the sphere:
    # z_c = z_bottom - R cos(theta)
    z_c = z_bottom - R * cos(th)

    z_upper = z_c + R

    mask = (z >= z_bottom) & (z <= z_upper)

    r_sq = np.maximum(
        0.0,
        R**2 - (z - z_c)**2,
    )

    # Circular xy cross-section divided by the full xy layer area.
    return pi * r_sq * mask / (l * ly)


""" filament_w3c0 """


def r_filament_w3c0_triangle(l, phi, th, delta):
    """
    Radius of the three identical rounded empty corners.

    The liquid occupies the effective triangle minus three
    identical empty corner regions.

    phi is referred to the original unit-side pore volume.
    """
    theta_prime = pi - th
    a_eff = a_eff_triangle(delta)

    corner_factor = (
        theta_prime
        - pi / 3.0
        - sin(2.0 * theta_prime - pi / 3.0)
        + sqrt(3.0) / 2.0
    )

    numerator = (
        (a_eff**2 - phi)
        * sqrt(3.0)
        / 12.0
    )

    return sqrt(numerator / corner_factor)


def rho_filament_w3c0_triangle(z, l, ly, phi, th, delta):
    """
    Analytical density profile for filament_w3c0.

    The configuration fills the central part of the triangular pore
    and leaves all three corners empty.

    z = 0 is the triangle centroid.
    """
    z = np.asarray(z, dtype=float)

    sqrt3 = sqrt(3.0)
    theta_prime = pi - th

    a_eff = a_eff_triangle(delta)
    H_eff = H_eff_triangle(delta)

    r = r_filament_w3c0_triangle(
        l, phi, th, delta
    )

    z_bottom = -H_eff / 3.0
    z_top = 2.0 * H_eff / 3.0

    # Centres of circles constructing the lower and upper rounded corners.
    zc_low = (
        z_bottom
        + r * cos(theta_prime)
    )

    zc_up = (
        z_top
        - 2.0 * r * cos(theta_prime)
    )

    # Distance from a vertex to the contact point along a wall.
    L1 = (
        2.0
        * r
        * sin(pi / 3.0 - theta_prime)
    )

    # Horizontal liquid width in the lower rounded part.
    L2 = (
        a_eff
        - 2.0 * sqrt3 * r * cos(theta_prime)
    )

    # Boundaries between the three analytical branches.
    z1 = z_bottom
    z2 = z_bottom + L1 * sqrt3 / 2.0
    z3 = z_top - L1 * sqrt3 / 2.0
    z4 = zc_up + r

    rho = np.zeros_like(z, dtype=float)

    # 1. Lower rounded corners.
    mask1 = (z >= z1) & (z <= z2)

    rho[mask1] = (
        L2
        + 2.0
        * sqrt(
            np.maximum(
                0.0,
                r**2 - (z[mask1] - zc_low)**2,
            )
        )
    ) / ly

    # 2. Full-width triangular section.
    mask2 = (z > z2) & (z <= z3)

    rho[mask2] = (
        a_eff
        - 2.0 / sqrt3
        * (z[mask2] - z_bottom)
    ) / ly

    # 3. Upper rounded corner.
    mask3 = (z > z3) & (z <= z4)

    rho[mask3] = (
        2.0
        * sqrt(
            np.maximum(
                0.0,
                r**2 - (z[mask3] - zc_up)**2,
            )
        )
        / ly
    )

    return rho


""" capsule_w3c0 """


def r_capsule_w3c0_triangle(l, phi, th, delta):
    """
    Meniscus radius for capsule_w3c0.
    """
    a_eff = a_eff_triangle(delta)

    return (
        -a_eff
        / (2.0 * sqrt(3.0) * cos(th))
    )


def d_capsule_w3c0_triangle(l, phi, th, delta):
    """
    Length of the cylindrical insertion between the two end menisci.

    This is the analytical expression used in the notebook,
    rewritten for the effective equilateral triangle.
    """
    sqrt3 = sqrt(3.0)

    a_eff = a_eff_triangle(delta)
    r = r_capsule_w3c0_triangle(
        l, phi, th, delta
    )

    numerator = (
        phi * l * sqrt3 / 4.0
        - pi * sqrt3 / 2.0 * r**2 * a_eff
        + 2.0 * pi * r**3 / 3.0
        + pi * sqrt3 / 72.0 * a_eff**3
    )

    arg = (
        a_eff
        / (2.0 * sqrt3 * r)
    )
    arg = np.clip(arg, -1.0, 1.0)

    root = sqrt(
        np.maximum(
            0.0,
            r**2 - a_eff**2 / 12.0,
        )
    )

    denominator = (
        r**2
        * (
            pi
            - 3.0 * np.arccos(arg)
        )
        + sqrt3 / 2.0
        * a_eff
        * root
    )

    return numerator / denominator


def rho_capsule_w3c0_triangle(z, l, ly, phi, th, delta):
    """
    Analytical density profile for capsule_w3c0.

    At a given z the xy cross-section is a stadium/capsule.
    When its circular part crosses the inclined triangle walls,
    the corresponding circular segments are subtracted analytically.
    """
    z = np.asarray(z, dtype=float)

    sqrt3 = sqrt(3.0)

    a_eff = a_eff_triangle(delta)
    H_eff = H_eff_triangle(delta)

    r = r_capsule_w3c0_triangle(
        l, phi, th, delta
    )

    d = d_capsule_w3c0_triangle(
        l, phi, th, delta
    )

    # Pore axis / centroid.
    z_c = 0.0

    z_bottom = -H_eff / 3.0
    z_top = 2.0 * H_eff / 3.0

    # Highest point of the spherical end-cap.
    z_cap_top = z_c + r

    rho = np.zeros_like(z, dtype=float)

    mask = (
        (z >= z_bottom)
        & (z <= z_cap_top)
    )

    if not np.any(mask):
        return rho

    z_m = z[mask]

    # Radius of the circular section of the capsule at height z.
    r_z = sqrt(
        np.maximum(
            0.0,
            r**2 - (z_m - z_c)**2,
        )
    )

    # Half-width of the equilateral pore at height z.
    pore_half_width = (
        z_top - z_m
    ) / sqrt3

    S = np.zeros_like(z_m)

    # --------------------------------------------------------
    # Case 1:
    # the capsule cross-section lies completely inside the pore
    # --------------------------------------------------------
    inside = r_z <= pore_half_width

    S[inside] = (
        2.0 * d * r_z[inside]
        + pi * r_z[inside]**2
    )

    # --------------------------------------------------------
    # Case 2:
    # the inclined pore walls cut the capsule cross-section
    # --------------------------------------------------------
    clipped = ~inside

    if np.any(clipped):
        rz = r_z[clipped]
        h = pore_half_width[clipped]

        ratio = np.divide(
            h,
            rz,
            out=np.ones_like(h),
            where=rz > 0,
        )
        ratio = np.clip(
            ratio,
            -1.0,
            1.0,
        )

        segment = (
            rz**2 * np.arccos(ratio)
            - h
            * sqrt(
                np.maximum(
                    0.0,
                    rz**2 - h**2,
                )
            )
        )

        S[clipped] = (
            2.0 * d * h
            + pi * rz**2
            - 2.0 * segment
        )

    rho[mask] = S / (l * ly)

    return rho


# ============================================================
# Analytical surface areas / capillary functional
# ============================================================
#
# These expressions are transferred from the stability-diagram
# notebook.  The volume fraction phi is always defined relative
# to the ORIGINAL unit-side triangular pore:
#
#     V = phi * (sqrt(3)/4) * l
#
# The wall layer only changes the accessible triangular side:
#
#     a' = 1 - 2*sqrt(3)*delta
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


def a_prime_triangle(delta):
    """Accessible triangle side after a parallel inward wall offset delta."""
    return 1.0 - 2.0 * sqrt(3.0) * delta


def _triangle_broadcast(l, phi):
    """Broadcast l and phi and remember whether the result was scalar."""
    l_arr, phi_arr = np.broadcast_arrays(
        np.asarray(l, dtype=float),
        np.asarray(phi, dtype=float),
    )
    scalar = l_arr.ndim == 0
    return l_arr, phi_arr, scalar


def _triangle_restore(value, scalar):
    arr = np.asarray(value)
    return float(arr) if scalar else arr


def _triangle_inf(l, phi):
    l_arr, phi_arr, scalar = _triangle_broadcast(l, phi)
    out = np.full(l_arr.shape, np.inf, dtype=float)
    return l_arr, phi_arr, out, scalar


def _triangle_energy(S_LV, S_SL, th, flag):
    if flag:
        return S_LV + S_SL
    return S_LV - cos(th) * S_SL


# ------------------------------------------------------------
# filament_w0c0
# ------------------------------------------------------------

def cond_filament_w0c0_triangle(l, phi, th, delta):
    del th
    l_arr, phi_arr, scalar = _triangle_broadcast(l, phi)
    a_prime = a_prime_triangle(delta)

    if a_prime <= 0:
        result = np.zeros(l_arr.shape, dtype=bool)
    else:
        phi_max = pi * a_prime**2 / (3.0 * sqrt(3.0))
        result = (phi_arr > 0.0) & (phi_arr <= phi_max)

    return bool(result) if scalar else result


def area_filament_w0c0_triangle(l, phi, th, delta, flag=False):
    del th, flag
    l_arr, phi_arr, S, scalar = _triangle_inf(l, phi)
    cond = np.asarray(
        cond_filament_w0c0_triangle(l_arr, phi_arr, 0.0, delta),
        dtype=bool,
    )

    r = sqrt(sqrt(3.0) * phi_arr / (4.0 * pi))
    values = 2.0 * pi * r * l_arr
    S = np.where(cond, values, S)

    return _triangle_restore(S, scalar)


# ------------------------------------------------------------
# filament_w1c0
# ------------------------------------------------------------

def cond_filament_w1c0_triangle(l, phi, th, delta):
    l_arr, phi_arr, scalar = _triangle_broadcast(l, phi)

    denom = th - sin(th) * cos(th)
    a_prime = a_prime_triangle(delta)

    if denom <= 0 or a_prime <= 0:
        result = np.zeros(l_arr.shape, dtype=bool)
        return bool(result) if scalar else result

    r = sqrt(phi_arr * sqrt(3.0) / (4.0 * denom))

    cond_base = 2.0 * r * sin(th) <= a_prime

    aa = -4.0 / 3.0
    bb = a_prime / sqrt(3.0) - 2.0 * r * cos(th)
    cc = (r * sin(th))**2 - 0.25 * a_prime**2

    z_max = r * (1.0 - cos(th))

    E_0 = cc
    E_zmax = aa * z_max**2 + bb * z_max + cc

    z_vert = -bb / (2.0 * aa)
    E_vert = aa * z_vert**2 + bb * z_vert + cc

    max_E = np.maximum(E_0, E_zmax)
    in_bounds = (z_vert >= 0.0) & (z_vert <= z_max)
    max_E = np.where(
        in_bounds,
        np.maximum(max_E, E_vert),
        max_E,
    )

    cond_walls = max_E < 0.0
    result = cond_base & cond_walls & (phi_arr > 0.0)

    return bool(result) if scalar else result


def area_filament_w1c0_triangle(l, phi, th, delta, flag=False):
    l_arr, phi_arr, S, scalar = _triangle_inf(l, phi)
    cond = np.asarray(
        cond_filament_w1c0_triangle(l_arr, phi_arr, th, delta),
        dtype=bool,
    )

    denom = th - sin(th) * cos(th)
    if denom <= 0:
        return _triangle_restore(S, scalar)

    r = sqrt(phi_arr * sqrt(3.0) / (4.0 * denom))

    P_LV = 2.0 * r * th
    P_SL = 2.0 * r * sin(th)

    values = l_arr * _triangle_energy(
        P_LV, P_SL, th, flag
    )
    S = np.where(cond, values, S)

    return _triangle_restore(S, scalar)


# ------------------------------------------------------------
# filament_w3c0
# ------------------------------------------------------------

def cond_filament_w3c0_triangle(l, phi, th, delta):
    l_arr, phi_arr, scalar = _triangle_broadcast(l, phi)
    a_prime = a_prime_triangle(delta)

    if a_prime <= 0 or th <= 2.0 * pi / 3.0:
        result = np.zeros(l_arr.shape, dtype=bool)
        return bool(result) if scalar else result

    th_p = pi - th
    denom = (
        th_p
        - pi / 3.0
        - sin(2.0 * th_p - pi / 3.0)
        + sqrt(3.0) / 2.0
    )

    if denom <= 0:
        result = np.zeros(l_arr.shape, dtype=bool)
        return bool(result) if scalar else result

    r_sq = (
        np.maximum(0.0, a_prime**2 - phi_arr)
        * (sqrt(3.0) / 12.0)
        / denom
    )
    r = sqrt(r_sq)
    x = 2.0 * r * sin(pi / 3.0 - th_p)

    cond_vol = (
        (phi_arr > 0.0)
        & (phi_arr < a_prime**2)
        & (r > 0.0)
    )
    cond_x = x <= a_prime / 2.0

    result = cond_vol & cond_x
    return bool(result) if scalar else result


def area_filament_w3c0_triangle(l, phi, th, delta, flag=False):
    l_arr, phi_arr, S, scalar = _triangle_inf(l, phi)
    cond = np.asarray(
        cond_filament_w3c0_triangle(l_arr, phi_arr, th, delta),
        dtype=bool,
    )

    a_prime = a_prime_triangle(delta)
    if a_prime <= 0 or th <= 2.0 * pi / 3.0:
        return _triangle_restore(S, scalar)

    th_p = pi - th
    denom = (
        th_p
        - pi / 3.0
        - sin(2.0 * th_p - pi / 3.0)
        + sqrt(3.0) / 2.0
    )
    if denom <= 0:
        return _triangle_restore(S, scalar)

    r = sqrt(
        np.maximum(0.0, a_prime**2 - phi_arr)
        * (sqrt(3.0) / 12.0)
        / denom
    )
    x = 2.0 * r * sin(pi / 3.0 - th_p)

    P_LV = 6.0 * r * (pi / 3.0 - th_p)
    P_SL = 3.0 * a_prime - 6.0 * x

    values = l_arr * _triangle_energy(
        P_LV, P_SL, th, flag
    )
    S = np.where(cond, values, S)

    return _triangle_restore(S, scalar)


# ------------------------------------------------------------
# filament_w2c0
# ------------------------------------------------------------

def cond_filament_w2c0_triangle(l, phi, th, delta):
    l_arr, phi_arr, scalar = _triangle_broadcast(l, phi)
    a_prime = a_prime_triangle(delta)

    if a_prime <= 0 or th <= 2.0 * pi / 3.0:
        result = np.zeros(l_arr.shape, dtype=bool)
        return bool(result) if scalar else result

    th_p = pi - th
    denom = pi - 2.0 * th_p + sin(2.0 * th_p)

    if denom <= 0:
        result = np.zeros(l_arr.shape, dtype=bool)
        return bool(result) if scalar else result

    r = sqrt(
        phi_arr * (sqrt(3.0) / 4.0) / denom
    )

    L2 = 2.0 * r * sin(pi / 3.0 + th_p)
    cond1 = L2 <= a_prime

    max_dist = r * (2.0 * cos(th_p) + 1.0)
    cond2 = (
        max_dist
        <= (sqrt(3.0) / 2.0) * a_prime
    )

    cond_vol = (
        (phi_arr > 0.0)
        & (phi_arr < a_prime**2)
    )

    result = cond1 & cond2 & cond_vol
    return bool(result) if scalar else result


def area_filament_w2c0_triangle(l, phi, th, delta, flag=False):
    l_arr, phi_arr, S, scalar = _triangle_inf(l, phi)
    cond = np.asarray(
        cond_filament_w2c0_triangle(l_arr, phi_arr, th, delta),
        dtype=bool,
    )

    if th <= 2.0 * pi / 3.0:
        return _triangle_restore(S, scalar)

    th_p = pi - th
    denom = pi - 2.0 * th_p + sin(2.0 * th_p)
    if denom <= 0:
        return _triangle_restore(S, scalar)

    r = sqrt(
        phi_arr * (sqrt(3.0) / 4.0) / denom
    )

    P_LV = 2.0 * r * (pi - 2.0 * th_p)
    P_SL = 4.0 * r * sin(th_p)

    values = l_arr * _triangle_energy(
        P_LV, P_SL, th, flag
    )
    S = np.where(cond, values, S)

    return _triangle_restore(S, scalar)


# ------------------------------------------------------------
# filament_w3c2
# ------------------------------------------------------------

def cond_filament_w3c2_triangle(
    l, phi, th, delta
):
    l_arr, phi_arr, scalar = _triangle_broadcast(l, phi)
    a_prime = a_prime_triangle(delta)

    if a_prime <= 0 or th <= 2.0 * pi / 3.0:
        result = np.zeros(l_arr.shape, dtype=bool)
        return bool(result) if scalar else result

    th_p = pi - th
    denom = (
        th_p
        - pi / 3.0
        - sin(2.0 * th_p - pi / 3.0)
        + sqrt(3.0) / 2.0
    )

    if denom <= 0:
        result = np.zeros(l_arr.shape, dtype=bool)
        return bool(result) if scalar else result

    r = sqrt(
        np.maximum(0.0, a_prime**2 - phi_arr)
        * (sqrt(3.0) / 8.0)
        / denom
    )
    x = 2.0 * r * sin(pi / 3.0 - th_p)

    cond_vol = (
        (phi_arr > 0.0)
        & (phi_arr < a_prime**2)
        & (r > 0.0)
    )
    cond_x = x <= a_prime / 2.0

    result = cond_vol & cond_x
    return bool(result) if scalar else result


def area_filament_w3c2_triangle(
    l, phi, th, delta, flag=False
):
    l_arr, phi_arr, S, scalar = _triangle_inf(l, phi)
    cond = np.asarray(
        cond_filament_w3c2_triangle(
            l_arr, phi_arr, th, delta
        ),
        dtype=bool,
    )

    a_prime = a_prime_triangle(delta)
    if a_prime <= 0 or th <= 2.0 * pi / 3.0:
        return _triangle_restore(S, scalar)

    th_p = pi - th
    denom = (
        th_p
        - pi / 3.0
        - sin(2.0 * th_p - pi / 3.0)
        + sqrt(3.0) / 2.0
    )
    if denom <= 0:
        return _triangle_restore(S, scalar)

    r = sqrt(
        np.maximum(0.0, a_prime**2 - phi_arr)
        * (sqrt(3.0) / 8.0)
        / denom
    )
    x = 2.0 * r * sin(pi / 3.0 - th_p)

    P_LV = 4.0 * r * (pi / 3.0 - th_p)
    P_SL = 3.0 * a_prime - 4.0 * x

    values = l_arr * _triangle_energy(
        P_LV, P_SL, th, flag
    )
    S = np.where(cond, values, S)

    return _triangle_restore(S, scalar)


# ------------------------------------------------------------
# filament_w2c1
# ------------------------------------------------------------

def cond_filament_w2c1_triangle(l, phi, th, delta):
    l_arr, phi_arr, scalar = _triangle_broadcast(l, phi)
    a_prime = a_prime_triangle(delta)

    # Kept as in the notebook implementation:
    # only the lower theta bound is enforced in this function.
    if a_prime <= 0 or th <= pi / 2.0:
        result = np.zeros(l_arr.shape, dtype=bool)
        return bool(result) if scalar else result

    th_p = pi - th
    denom = (
        2.0 * pi / 3.0
        - th_p
        + sqrt(3.0) / 2.0
        + sin(pi / 3.0 + 2.0 * th_p)
    )

    if denom <= 0:
        result = np.zeros(l_arr.shape, dtype=bool)
        return bool(result) if scalar else result

    r = sqrt(
        phi_arr * (sqrt(3.0) / 4.0) / denom
    )
    x = 2.0 * r * sin(pi / 3.0 + th_p)

    cond1 = x <= a_prime

    max_dist = r * (2.0 * cos(th_p) + 1.0)
    cond2 = (
        max_dist
        <= (sqrt(3.0) / 2.0) * a_prime
    )

    cond_vol = (phi_arr > 0.0) & (phi_arr < 1.0)

    result = cond1 & cond2 & cond_vol
    return bool(result) if scalar else result


def area_filament_w2c1_triangle(
    l, phi, th, delta, flag=False
):
    l_arr, phi_arr, S, scalar = _triangle_inf(l, phi)
    cond = np.asarray(
        cond_filament_w2c1_triangle(
            l_arr, phi_arr, th, delta
        ),
        dtype=bool,
    )

    if th <= pi / 2.0:
        return _triangle_restore(S, scalar)

    th_p = pi - th
    denom = (
        2.0 * pi / 3.0
        - th_p
        + sqrt(3.0) / 2.0
        + sin(pi / 3.0 + 2.0 * th_p)
    )
    if denom <= 0:
        return _triangle_restore(S, scalar)

    r = sqrt(
        phi_arr * (sqrt(3.0) / 4.0) / denom
    )
    x = 2.0 * r * sin(pi / 3.0 + th_p)

    P_LV = r * (4.0 * pi / 3.0 - 2.0 * th_p)
    P_SL = 2.0 * x

    values = l_arr * _triangle_energy(
        P_LV, P_SL, th, flag
    )
    S = np.where(cond, values, S)

    return _triangle_restore(S, scalar)


# ------------------------------------------------------------
# filament_w3c1
# ------------------------------------------------------------

def cond_filament_w3c1_triangle(l, phi, th, delta):
    l_arr, phi_arr, scalar = _triangle_broadcast(l, phi)
    a_prime = a_prime_triangle(delta)

    if a_prime <= 0 or th <= 2.0 * pi / 3.0:
        result = np.zeros(l_arr.shape, dtype=bool)
        return bool(result) if scalar else result

    th_p = pi - th
    denom = (
        th_p
        - pi / 3.0
        - sin(2.0 * th_p - pi / 3.0)
        + sqrt(3.0) / 2.0
    )

    if denom <= 0:
        result = np.zeros(l_arr.shape, dtype=bool)
        return bool(result) if scalar else result

    r = sqrt(
        np.maximum(0.0, a_prime**2 - phi_arr)
        * (sqrt(3.0) / 4.0)
        / denom
    )
    x = 2.0 * r * sin(pi / 3.0 - th_p)

    cond_vol = (
        (phi_arr > 0.0)
        & (phi_arr < a_prime**2)
        & (r > 0.0)
    )
    cond_x = x <= a_prime

    result = cond_vol & cond_x
    return bool(result) if scalar else result


def area_filament_w3c1_triangle(
    l, phi, th, delta, flag=False
):
    l_arr, phi_arr, S, scalar = _triangle_inf(l, phi)
    cond = np.asarray(
        cond_filament_w3c1_triangle(
            l_arr, phi_arr, th, delta
        ),
        dtype=bool,
    )

    a_prime = a_prime_triangle(delta)
    if a_prime <= 0 or th <= 2.0 * pi / 3.0:
        return _triangle_restore(S, scalar)

    th_p = pi - th
    denom = (
        th_p
        - pi / 3.0
        - sin(2.0 * th_p - pi / 3.0)
        + sqrt(3.0) / 2.0
    )
    if denom <= 0:
        return _triangle_restore(S, scalar)

    r = sqrt(
        np.maximum(0.0, a_prime**2 - phi_arr)
        * (sqrt(3.0) / 4.0)
        / denom
    )
    x = 2.0 * r * sin(pi / 3.0 - th_p)

    P_LV = 2.0 * r * (pi / 3.0 - th_p)
    P_SL = 3.0 * a_prime - 2.0 * x

    values = l_arr * _triangle_energy(
        P_LV, P_SL, th, flag
    )
    S = np.where(cond, values, S)

    return _triangle_restore(S, scalar)


# ------------------------------------------------------------
# droplet_w1c0
# ------------------------------------------------------------

def cond_droplet_w1c0_triangle(l, phi, th, delta):
    l_arr, phi_arr, scalar = _triangle_broadcast(l, phi)
    a_prime = a_prime_triangle(delta)

    if a_prime <= 0 or th < pi / 2.0:
        result = np.zeros(l_arr.shape, dtype=bool)
        return bool(result) if scalar else result

    vol_factor = (
        pi / 3.0
        * (2.0 - 3.0 * cos(th) + cos(th)**3)
    )
    if vol_factor <= 0:
        result = np.zeros(l_arr.shape, dtype=bool)
        return bool(result) if scalar else result

    V = phi_arr * (sqrt(3.0) / 4.0) * l_arr
    R = (V / vol_factor) ** (1.0 / 3.0)

    cond_length = 2.0 * R <= l_arr
    cond_base = 2.0 * R * sin(th) <= a_prime

    h = R * (1.0 - cos(th))
    A = 4.0 / 3.0
    B = a_prime / sqrt(3.0) - 2.0 * R * cos(th)
    C = -0.25 * a_prime**2 + (R * sin(th))**2

    f_0 = C
    f_h = -A * h**2 + B * h + C
    z_v = B / (2.0 * A)
    f_v = -A * z_v**2 + B * z_v + C

    f_max = np.maximum(f_0, f_h)
    in_range = (z_v >= 0.0) & (z_v <= h)
    f_max = np.where(
        in_range,
        np.maximum(f_max, f_v),
        f_max,
    )

    cond_no_touch = f_max < 0.0
    cond_vol = (phi_arr > 0.0) & (phi_arr < 1.0)

    result = (
        cond_length
        & cond_base
        & cond_no_touch
        & cond_vol
    )
    return bool(result) if scalar else result


def area_droplet_w1c0_triangle(
    l, phi, th, delta, flag=False
):
    l_arr, phi_arr, S, scalar = _triangle_inf(l, phi)
    cond = np.asarray(
        cond_droplet_w1c0_triangle(
            l_arr, phi_arr, th, delta
        ),
        dtype=bool,
    )

    if th < pi / 2.0:
        return _triangle_restore(S, scalar)

    vol_factor = (
        pi / 3.0
        * (2.0 - 3.0 * cos(th) + cos(th)**3)
    )
    if vol_factor <= 0:
        return _triangle_restore(S, scalar)

    V = phi_arr * (sqrt(3.0) / 4.0) * l_arr
    R = (V / vol_factor) ** (1.0 / 3.0)

    S_LV = 2.0 * pi * R**2 * (1.0 - cos(th))
    S_SL = pi * R**2 * sin(th)**2

    values = _triangle_energy(
        S_LV, S_SL, th, flag
    )
    S = np.where(cond, values, S)

    return _triangle_restore(S, scalar)


# ------------------------------------------------------------
# droplet_w2c0
# ------------------------------------------------------------

def cond_droplet_w2c0_triangle(l, phi, th, delta):
    l_arr, phi_arr, scalar = _triangle_broadcast(l, phi)
    a_prime = a_prime_triangle(delta)
    H_prime = sqrt(3.0) * a_prime / 2.0

    if a_prime <= 0 or th <= 2.0 * pi / 3.0:
        result = np.zeros(l_arr.shape, dtype=bool)
        return bool(result) if scalar else result

    vol_factor = (
        2.0 * pi / 3.0
        * (-cos(th))
        * (3.0 - cos(th)**2)
    )
    if vol_factor <= 0:
        result = np.zeros(l_arr.shape, dtype=bool)
        return bool(result) if scalar else result

    V = phi_arr * (sqrt(3.0) / 4.0) * l_arr
    R = (V / vol_factor) ** (1.0 / 3.0)

    cond_length = 2.0 * R <= l_arr
    cond_no_touch_3rd = (
        R * (1.0 - 2.0 * cos(th)) <= H_prime
    )
    cond_width = (
        R * (sin(th) + sqrt(3.0) * cos(th))
        <= a_prime
    )
    cond_vol = (phi_arr > 0.0) & (phi_arr < 1.0)

    result = (
        cond_length
        & cond_no_touch_3rd
        & cond_width
        & cond_vol
    )
    return bool(result) if scalar else result


def area_droplet_w2c0_triangle(
    l, phi, th, delta, flag=False
):
    l_arr, phi_arr, S, scalar = _triangle_inf(l, phi)
    cond = np.asarray(
        cond_droplet_w2c0_triangle(
            l_arr, phi_arr, th, delta
        ),
        dtype=bool,
    )

    if th <= 2.0 * pi / 3.0:
        return _triangle_restore(S, scalar)

    vol_factor = (
        2.0 * pi / 3.0
        * (-cos(th))
        * (3.0 - cos(th)**2)
    )
    if vol_factor <= 0:
        return _triangle_restore(S, scalar)

    V = phi_arr * (sqrt(3.0) / 4.0) * l_arr
    R = (V / vol_factor) ** (1.0 / 3.0)

    S_LV = -4.0 * pi * R**2 * cos(th)
    S_SL = 2.0 * pi * R**2 * sin(th)**2

    values = _triangle_energy(
        S_LV, S_SL, th, flag
    )
    S = np.where(cond, values, S)

    return _triangle_restore(S, scalar)


# ------------------------------------------------------------
# droplet_w3c0
# ------------------------------------------------------------

def cond_droplet_w3c0_triangle(l, phi, th, delta):
    l_arr, phi_arr, scalar = _triangle_broadcast(l, phi)
    a_prime = a_prime_triangle(delta)

    if a_prime <= 0 or th <= 2.0 * pi / 3.0:
        result = np.zeros(l_arr.shape, dtype=bool)
        return bool(result) if scalar else result

    R = -a_prime / (2.0 * sqrt(3.0) * cos(th))
    if R <= 0:
        result = np.zeros(l_arr.shape, dtype=bool)
        return bool(result) if scalar else result

    V_3 = (
        pi * R**3 / 3.0
        * (
            3.0 * cos(th)**3
            - 9.0 * cos(th)
            - 2.0
        )
    )
    V_grid = phi_arr * (sqrt(3.0) / 4.0) * l_arr

    cond_length = 2.0 * R <= l_arr
    cond_vol = (
        np.abs(V_grid - V_3)
        / np.maximum(np.abs(V_3), 1e-12)
        < 0.015
    )
    cond_vol &= (phi_arr > 0.0) & (phi_arr < 1.0)

    result = cond_length & cond_vol
    return bool(result) if scalar else result


def area_droplet_w3c0_triangle(
    l, phi, th, delta, flag=False
):
    l_arr, phi_arr, S, scalar = _triangle_inf(l, phi)
    cond = np.asarray(
        cond_droplet_w3c0_triangle(
            l_arr, phi_arr, th, delta
        ),
        dtype=bool,
    )

    a_prime = a_prime_triangle(delta)
    if a_prime <= 0 or th <= 2.0 * pi / 3.0:
        return _triangle_restore(S, scalar)

    R = -a_prime / (2.0 * sqrt(3.0) * cos(th))
    if R <= 0:
        return _triangle_restore(S, scalar)

    S_LV = (
        -2.0 * pi * R**2
        * (1.0 + 3.0 * cos(th))
    )
    S_SL = 3.0 * pi * (R * sin(th))**2

    values = _triangle_energy(
        S_LV, S_SL, th, flag
    )
    S = np.where(cond, values, S)

    return _triangle_restore(S, scalar)


# ------------------------------------------------------------
# droplet_w2c1
# ------------------------------------------------------------

def cond_droplet_w2c1_triangle(
    l, phi, th, delta
):
    l_arr, phi_arr, scalar = _triangle_broadcast(l, phi)
    a_prime = a_prime_triangle(delta)
    H_prime = sqrt(3.0) * a_prime / 2.0

    if (
        a_prime <= 0
        or th <= pi / 2.0
        or th >= 2.0 * pi / 3.0
    ):
        result = np.zeros(l_arr.shape, dtype=bool)
        return bool(result) if scalar else result

    s_lv_factor = (
        2.0 * pi
        * (
            1.0
            + cos(th)
            - 2.0 * sin(th) / 3.0
        )
    )
    s_sl_factor = 2.0 * pi * sin(th)**2

    vol_factor = (
        s_lv_factor - s_sl_factor * cos(th)
    ) / 3.0
    if vol_factor <= 0:
        result = np.zeros(l_arr.shape, dtype=bool)
        return bool(result) if scalar else result

    V = phi_arr * (sqrt(3.0) / 4.0) * l_arr
    R = (V / vol_factor) ** (1.0 / 3.0)

    cond_length = 2.0 * R <= l_arr
    cond_no_touch_3rd = (
        R * (1.0 - 2.0 * cos(th)) <= H_prime
    )
    cond_width = (
        R * (sin(th) - sqrt(3.0) * cos(th))
        <= a_prime
    )
    cond_vol = (
        (phi_arr > 0.0)
        & (phi_arr < a_prime**2)
    )

    result = (
        cond_length
        & cond_no_touch_3rd
        & cond_width
        & cond_vol
    )
    return bool(result) if scalar else result


def area_droplet_w2c1_triangle(
    l, phi, th, delta, flag=False
):
    l_arr, phi_arr, S, scalar = _triangle_inf(l, phi)
    cond = np.asarray(
        cond_droplet_w2c1_triangle(
            l_arr, phi_arr, th, delta
        ),
        dtype=bool,
    )

    if (
        th <= pi / 2.0
        or th >= 2.0 * pi / 3.0
    ):
        return _triangle_restore(S, scalar)

    s_lv_factor = (
        2.0 * pi
        * (
            1.0
            + cos(th)
            - 2.0 * sin(th) / 3.0
        )
    )
    s_sl_factor = 2.0 * pi * sin(th)**2

    vol_factor = (
        s_lv_factor - s_sl_factor * cos(th)
    ) / 3.0
    if vol_factor <= 0:
        return _triangle_restore(S, scalar)

    V = phi_arr * (sqrt(3.0) / 4.0) * l_arr
    R = (V / vol_factor) ** (1.0 / 3.0)

    S_LV = s_lv_factor * R**2
    S_SL = s_sl_factor * R**2

    values = _triangle_energy(
        S_LV, S_SL, th, flag
    )
    S = np.where(cond, values, S)

    return _triangle_restore(S, scalar)


# ------------------------------------------------------------
# capsule_w3c3
# ------------------------------------------------------------

def cond_capsule_w3c3_triangle(l, phi, th, delta):
    l_arr, phi_arr, scalar = _triangle_broadcast(l, phi)
    a_prime = a_prime_triangle(delta)

    if (
        a_prime <= 0
        or th <= pi / 2.0
        or th >= 2.0 * pi / 3.0
    ):
        result = np.zeros(l_arr.shape, dtype=bool)
        return bool(result) if scalar else result

    c = cos(th)
    r = -a_prime / (2.0 * sqrt(3.0) * c)
    if r <= 0:
        result = np.zeros(l_arr.shape, dtype=bool)
        return bool(result) if scalar else result

    denom_sqrt = 1.0 - 4.0 * c**2
    if denom_sqrt <= 0:
        result = np.zeros(l_arr.shape, dtype=bool)
        return bool(result) if scalar else result

    arg = (
        (1.0 - 6.0 * c**2)
        / (2.0 - 6.0 * c**2)
    )
    arg = np.clip(arg, -1.0, 1.0)

    V_ends = (
        -(3.0 * np.arccos(arg) - pi)
        / (36.0 * sqrt(3.0) * c**3)
        + c / (3.0 * sqrt(denom_sqrt))
    )
    V_ends_scaled = V_ends * a_prime**3

    L_ins = (
        phi_arr * l_arr
        - 4.0 * V_ends_scaled / sqrt(3.0)
    )

    result = (
        (L_ins >= 0.0)
        & (L_ins + 2.0 * r <= l_arr)
        & (phi_arr > 0.0)
    )
    return bool(result) if scalar else result


def area_capsule_w3c3_triangle(
    l, phi, th, delta, flag=False
):
    l_arr, phi_arr, S, scalar = _triangle_inf(l, phi)
    cond = np.asarray(
        cond_capsule_w3c3_triangle(
            l_arr, phi_arr, th, delta
        ),
        dtype=bool,
    )

    a_prime = a_prime_triangle(delta)
    if (
        a_prime <= 0
        or th <= pi / 2.0
        or th >= 2.0 * pi / 3.0
    ):
        return _triangle_restore(S, scalar)

    c = cos(th)
    t = np.tan(th)
    denom_sqrt = 1.0 - 4.0 * c**2
    if denom_sqrt <= 0:
        return _triangle_restore(S, scalar)

    arg_acos = np.clip(
        (1.0 - 6.0 * c**2)
        / (2.0 - 6.0 * c**2),
        -1.0,
        1.0,
    )

    V_ends = (
        -(3.0 * np.arccos(arg_acos) - pi)
        / (36.0 * sqrt(3.0) * c**3)
        + c / (3.0 * sqrt(denom_sqrt))
    )
    V_ends_scaled = V_ends * a_prime**3

    L_ins = (
        phi_arr * l_arr
        - 4.0 * V_ends_scaled / sqrt(3.0)
    )

    arg_asin = np.clip(
        sqrt(3.0) / t,
        -1.0,
        1.0,
    )

    S_contact_ends = (
        -sqrt(3.0 * denom_sqrt) / c
        - t**2 * np.arcsin(arg_asin)
    ) * a_prime**2

    S_contact_body = 3.0 * L_ins

    S_free_ends = (
        (3.0 * np.arccos(arg_acos) - pi)
        / (6.0 * c**2)
    ) * a_prime**2

    if flag:
        values = (
            S_contact_ends
            + S_contact_body
            + S_free_ends
        )
    else:
        values = (
            -c
            * (
                S_contact_ends
                + S_contact_body
            )
            + S_free_ends
        )

    S = np.where(cond, values, S)
    return _triangle_restore(S, scalar)


# ------------------------------------------------------------
# capsule_w3c0
# ------------------------------------------------------------

def cond_capsule_w3c0_triangle(l, phi, th, delta=0.0, flag=False):
    """
    Условие существования формы capsule_w3c0
    в равносторонней треугольной поре
    с учетом пристеночного слоя delta.
    """
    del flag

    phi = np.asarray(phi, dtype=float)
    l = np.asarray(l, dtype=float)

    phi, l = np.broadcast_arrays(phi, l)

    a_prime = 1.0 - 2.0 * np.sqrt(3.0) * delta

    if a_prime <= 0:
        return np.zeros_like(phi, dtype=bool)

    if th <= 2.0 * np.pi / 3.0 or th > np.pi:
        return np.zeros_like(phi, dtype=bool)

    cos_t = np.cos(th)
    sin_t = np.sin(th)

    # Радиус мениска
    r = -a_prime / (
        2.0 * np.sqrt(3.0) * cos_t
    )

    if r <= 0:
        return np.zeros_like(phi, dtype=bool)

    # Объем двух торцевых частей
    V_ends = (
        np.pi
        * (
            2.0
            + 9.0 * cos_t
            - 3.0 * cos_t**3
        )
        / (
            72.0
            * np.sqrt(3.0)
            * cos_t**3
        )
    )

    V_ends *= a_prime**3

    # Площадь поперечного сечения
    A_cross = (
        3.0 * th
        - 2.0 * np.pi
        - 3.0 * sin_t * cos_t
    ) / (
        12.0 * cos_t**2
    )

    A_cross *= a_prime**2

    if A_cross <= 0:
        return np.zeros_like(phi, dtype=bool)

    # Длина цилиндрической части
    l_tilde = (
        (np.sqrt(3.0) / 4.0) * phi * l
        - V_ends
    ) / A_cross

    cond = (
        (l_tilde >= 0.0)
        & (l_tilde + 2.0 * r <= l)
        & (phi > 0.0)
    )

    return cond


def area_capsule_w3c0_triangle(
    l,
    phi,
    th,
    delta=0.0,
    flag=False,
):
    """
    Площадь поверхности / капиллярный функционал
    для формы capsule_w3c0.

    flag=False:
        S = S_lg - cos(th) * S_sl

    flag=True:
        S = S_lg + S_sl
    """
    phi = np.asarray(phi, dtype=float)
    l = np.asarray(l, dtype=float)

    phi, l = np.broadcast_arrays(phi, l)

    S = np.full_like(
        phi,
        np.inf,
        dtype=float,
    )

    cond = cond_capsule_w3c0_triangle(
        l,
        phi,
        th,
        delta=delta,
    )

    if not np.any(cond):
        return S

    # Расстояние от центра треугольника
    # до эффективной стенки:
    #
    # a' = 1 - 2*sqrt(3)*delta
    #
    # h = a'/(2*sqrt(3))
    #   = sqrt(3)/6 - delta
    h = np.sqrt(3.0) / 6.0 - delta

    # Радиус мениска
    r = -h / np.cos(th)

    # Объем жидкости
    V_water = (
        np.sqrt(3.0) / 4.0
        * phi[cond]
        * l[cond]
    )

    # Площадь поперечного сечения
    # центральной части capsule_w3c0
    S_cross = (
        np.pi * r**2
        - 3.0
        * (
            r**2 * np.arccos(h / r)
            - h * np.sqrt(r**2 - h**2)
        )
    )

    # Объем двух торцевых частей
    V_ends = (
        3.0 * np.pi * r**2 * h
        - (2.0 / 3.0) * np.pi * r**3
        - np.pi * h**3
    )

    # Длина цилиндрической части
    l_tilde = (
        V_water - V_ends
    ) / S_cross

    # ----------------------------------------
    # Liquid-gas / liquid-liquid interface
    # ----------------------------------------

    S_lg = (
        2.0
        * np.pi
        * r
        * (3.0 * h - r)
        + 2.0
        * r
        * (
            np.pi
            - 3.0 * np.arccos(h / r)
        )
        * l_tilde
    )

    # ----------------------------------------
    # Solid-liquid interface
    # ----------------------------------------

    S_sl = (
        6.0
        * l_tilde
        * np.sqrt(r**2 - h**2)
        + 3.0
        * np.pi
        * (r**2 - h**2)
    )

    if flag:
        # Геометрическая сумма площадей
        res = S_lg + S_sl
    else:
        # Капиллярный функционал
        res = (
            S_lg
            - S_sl * np.cos(th)
        )

    # Дополнительная защита от
    # нефизической отрицательной длины
    res = np.where(
        l_tilde >= 0.0,
        res,
        np.inf,
    )

    S[cond] = res

    return S


# ------------------------------------------------------------
# filament_w0c0_straight
# ------------------------------------------------------------

def cond_filament_w0c0_straight_triangle(l, phi, th, delta):
    del th
    l_arr, phi_arr, scalar = _triangle_broadcast(l, phi)
    a_prime = a_prime_triangle(delta)

    if a_prime <= 0:
        result = np.zeros(l_arr.shape, dtype=bool)
    else:
        result = (
            (phi_arr > 0.0)
            & (sqrt(phi_arr) <= a_prime)
        )

    return bool(result) if scalar else result


def area_filament_w0c0_straight_triangle(
    l, phi, th, delta, flag=False
):
    del th, flag
    l_arr, phi_arr, S, scalar = _triangle_inf(l, phi)
    cond = np.asarray(
        cond_filament_w0c0_straight_triangle(
            l_arr, phi_arr, 0.0, delta
        ),
        dtype=bool,
    )

    values = 3.0 * sqrt(phi_arr) * l_arr
    S = np.where(cond, values, S)

    return _triangle_restore(S, scalar)


# ------------------------------------------------------------
# PANDA-style S_* aliases
# ------------------------------------------------------------

S_filament_w0c0_triangle = area_filament_w0c0_triangle
S_filament_w1c0_triangle = area_filament_w1c0_triangle
S_filament_w2c0_triangle = area_filament_w2c0_triangle
S_filament_w3c0_triangle = area_filament_w3c0_triangle
S_filament_w2c1_triangle = area_filament_w2c1_triangle
S_filament_w3c2_triangle = area_filament_w3c2_triangle
S_filament_w3c1_triangle = area_filament_w3c1_triangle
S_droplet_w1c0_triangle = area_droplet_w1c0_triangle
S_droplet_w2c0_triangle = area_droplet_w2c0_triangle
S_droplet_w3c0_triangle = area_droplet_w3c0_triangle
S_droplet_w2c1_triangle = area_droplet_w2c1_triangle
S_capsule_w3c3_triangle = area_capsule_w3c3_triangle
S_capsule_w3c0_triangle = area_capsule_w3c0_triangle
S_filament_w0c0_straight_triangle = area_filament_w0c0_straight_triangle


# ------------------------------------------------------------
# Registries
# ------------------------------------------------------------

AREA_TRIANGLE = {
    "filament_w0c0": area_filament_w0c0_triangle,
    "filament_w1c0": area_filament_w1c0_triangle,
    "filament_w2c0": area_filament_w2c0_triangle,
    "filament_w3c0": area_filament_w3c0_triangle,
    "filament_w2c1": area_filament_w2c1_triangle,
    "filament_w3c2": area_filament_w3c2_triangle,
    "filament_w3c1": area_filament_w3c1_triangle,
    "droplet_w1c0": area_droplet_w1c0_triangle,
    "droplet_w2c0": area_droplet_w2c0_triangle,
    "droplet_w3c0": area_droplet_w3c0_triangle,
    "droplet_w2c1": area_droplet_w2c1_triangle,
    "capsule_w3c3": area_capsule_w3c3_triangle,
    "capsule_w3c0": area_capsule_w3c0_triangle,
    "filament_w0c0_straight": area_filament_w0c0_straight_triangle,
}

COND_TRIANGLE = {
    "filament_w0c0": cond_filament_w0c0_triangle,
    "filament_w1c0": cond_filament_w1c0_triangle,
    "filament_w2c0": cond_filament_w2c0_triangle,
    "filament_w3c0": cond_filament_w3c0_triangle,
    "filament_w2c1": cond_filament_w2c1_triangle,
    "filament_w3c2": cond_filament_w3c2_triangle,
    "filament_w3c1": cond_filament_w3c1_triangle,
    "droplet_w1c0": cond_droplet_w1c0_triangle,
    "droplet_w2c0": cond_droplet_w2c0_triangle,
    "droplet_w3c0": cond_droplet_w3c0_triangle,
    "droplet_w2c1": cond_droplet_w2c1_triangle,
    "capsule_w3c3": cond_capsule_w3c3_triangle,
    "capsule_w3c0": cond_capsule_w3c0_triangle,
    "filament_w0c0_straight": cond_filament_w0c0_straight_triangle,
}

RHO_TRIANGLE = {
    "filament_w0c0": rho_filament_w0c0_triangle,
    "droplet_w1c0": rho_droplet_w1c0_triangle,
    "filament_w3c0": rho_filament_w3c0_triangle,
    "capsule_w3c0": rho_capsule_w3c0_triangle,
}


def get_area_triangle(name):
    return AREA_TRIANGLE[name]


def get_cond_triangle(name):
    return COND_TRIANGLE[name]


def get_rho_triangle(name):
    return RHO_TRIANGLE[name]


__all__ = [
    # common geometry
    "BIG_NUM",
    "a_eff_triangle",
    "a_prime_triangle",
    "H_eff_triangle",
    "z_bottom_triangle",
    "z_top_triangle",

    # density profiles already implemented
    "rho_filament_w0c0_triangle",
    "rho_droplet_w1c0_triangle",
    "rho_filament_w3c0_triangle",
    "rho_capsule_w3c0_triangle",

    # area/condition registries
    "AREA_TRIANGLE",
    "COND_TRIANGLE",
    "RHO_TRIANGLE",
    "get_area_triangle",
    "get_cond_triangle",
    "get_rho_triangle",
]

__all__ += [
    name
    for name in globals()
    if name.startswith(("area_", "S_", "cond_", "r_", "d_"))
    and name.endswith("_triangle")
]
