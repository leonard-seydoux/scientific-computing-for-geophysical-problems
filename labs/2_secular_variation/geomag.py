"""Geomagnetic field module.

This module contains functions to compute the geomagnetic field and its
derivatives at a given point on the Earth's surface, given a set of spherical
harmonic coefficients.

The functions in this module are based on the formulas given in [1]. The Python
implementation only relies on the NumPy library. The module also contains
functions to plot the geomagnetic field on a map. The maps are created using the
Cartopy library.

References
----------
.. [1] Backus, G., Parker, R. L., & Constable, C. (1996). Foundations of
    geomagnetism. Cambridge University Press.

Notes
-----
This module was originally written in Fortran by Vincent Lesur, transformed into
Python by Vincent Lesur and Alexandre Fournier, and simplified by Leonard
Seydoux. It is part of the class Scientific Computing for Geophysical Problems,
for the master class at the institut de physique du globe de Paris.

Last update: 2023.
"""

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np

# Maximum spherical harmonic degree expension
MAX_DEGREE = 8


# Total number of spherical harmonic coefficients
N_COEFFICIENTS = MAX_DEGREE * (MAX_DEGREE + 2)

# Reference radius
EARTH_RADIUS_KM = 6371.2

__all__ = [
    "associate_legendre_functions_derivative",
    "associated_legendre_functions",
    "data_weights_matrix",
    "damping_matrix",
    "forward_matrix",
    "log_factorial",
    "weights_matrix",
]


def weights_matrix(dim, sigma=1):
    """Weights matrix.

    Compute the data weights matrix for a given number of data points and a
    given standard deviation of the data.

    Parameters
    ----------
    dim : int
        Number of data points.
    sigma : float
        Standard deviation of the data in nT/year. Default is 1 nT/year.

    Returns
    -------
    ndarray of shape (dim, dim)
        Data weights matrix.
    """
    return (1 / sigma) ** 2 * np.eye(dim)


def damping_matrix(
    n_coefficients,
    radius,
    reference_radius=EARTH_RADIUS_KM,
    max_degree=MAX_DEGREE,
):
    """Damping matrix.

    Compute the damping matrix for a given radius.

    Parameters
    ----------
    n_coefficients : int
        Number of spherical harmonic coefficients.
    radius : float
        Radius in kilometers.
    reference_radius : float, optional
        Reference radius in kilometers. Default is the Earth radius, as defined
        in the EARTH_RADIUS_KM constant in this module.
    max_degree : int, optional
        Maximum degree of the spherical harmonic expansion. Default is 8, as
        defined in the MAX_DEGREE constant in this module.

    Returns
    -------
    ndarray of shape (n_coefficients, n_coefficients)
        Damping matrix.
    """
    # Initialize the damping matrix
    matrix = np.zeros((n_coefficients), dtype=float)
    radius = reference_radius / radius

    # m = 0
    for l in range(1, max_degree + 1):
        k = l * l - 1
        matrix[k] = (l + 1) * np.power(radius, 2 * l + 4)

    #  m != 0
    for m in range(1, max_degree + 1):
        for l in range(m, max_degree + 1):
            k = l * l + 2 * m - 2
            matrix[k] = (l + 1) * np.power(radius, 2 * l + 4)
            matrix[k + 1] = (l + 1) * np.power(radius, 2 * l + 4)

    return np.diag(matrix)


def forward_matrix(
    radius,
    colatitude,
    longitude,
    component="x",
    reference_radius=EARTH_RADIUS_KM,
    max_degree=MAX_DEGREE,
):
    """Forward matrix.

    Computes the forward matrix for a given radius, colatitude and longitude.
    The forward matrix is the matrix that relates the spherical harmonic
    coefficients to the magnetic field components at a given point.

    Parameters
    ----------
    radius : float
        Radius in kilometers.
    colatitude : float
        Colatitude in degrees.
    longitude : float
        Longitude in degrees.
    component : str, optional
        Component of the magnetic field. Can be "x", "y" or "z". Default is
        "x".
    reference_radius : float, optional
        Reference radius in kilometers. Default is the Earth radius, as defined
        in the EARTH_RADIUS_KM constant in this module.
    max_degree : int, optional
        Maximum degree of the spherical harmonic expansion. Default is 8, as
        defined in the MAX_DEGREE constant in this module.

    Returns
    -------
    ndarray of shape (N_COEFFICIENTS,)
        Array containing the forward matrix for component of the magnetic field
        at the given point.
    """
    # Split the case where all components are computed
    if component.lower() == "all":
        ax = forward_matrix(radius, colatitude, longitude, component="x")
        ay = forward_matrix(radius, colatitude, longitude, component="y")
        az = forward_matrix(radius, colatitude, longitude, component="z")
        a = np.vstack((ax, ay, az))
        return a

    else:
        # Initialize the forward matrix
        a = np.zeros(N_COEFFICIENTS, dtype=float)

        # Normalize the radius
        radius = reference_radius / radius
        theta = np.radians(colatitude)
        phi = np.radians(longitude)

        if component.lower() == "x":
            # Compute the forward matrix
            m = 0
            f = associate_legendre_functions_derivative(m, theta)
            for l in range(1, max_degree + 1, 1):
                k = l * l - 1
                a[k] = f[l] * np.power(radius, l + 2)
            for m in range(1, max_degree + 1, 1):
                f = associate_legendre_functions_derivative(m, theta)
                cos_phi = np.cos(m * phi)
                sin_phi = np.sin(m * phi)
                for l in range(m, max_degree + 1, 1):
                    k = l * l + 2 * m - 2
                    w = f[l] * np.power(radius, l + 2)
                    a[k] = w * cos_phi
                    a[k + 1] = w * sin_phi

        elif component.lower() == "y":
            # Condition for colatitude = 0
            sin_theta = np.sin(theta)
            sin_theta = 1.0e-10 if sin_theta == 0 else sin_theta

            # Compute the forward matrix
            for m in range(1, max_degree + 1, 1):
                f = associated_legendre_functions(m, theta)
                cos_phi = -m * np.sin(m * phi)
                sin_phi = m * np.cos(m * phi)
                for il in range(m, max_degree + 1, 1):
                    k = il * il + 2 * m - 2
                    ww = f[il] * np.power(radius, il + 2) / sin_theta
                    a[k] = -ww * cos_phi
                    a[k + 1] = -ww * sin_phi

        elif component.lower() == "z":
            m = 0
            f = associated_legendre_functions(m, theta)
            for il in range(1, max_degree + 1, 1):
                k = il * il - 1
                a[k] = -float(il + 1) * f[il] * np.power(radius, float(il + 2))
            #  m != 0
            for m in range(1, max_degree + 1, 1):
                f = associated_legendre_functions(m, theta)
                cos_phi = np.cos(m * phi)
                sin_phi = np.sin(m * phi)
                for il in range(m, max_degree + 1, 1):
                    k = il * il + 2 * m - 2
                    ww = (
                        -float(il + 1) * f[il] * np.power(radius, float(il + 2))
                    )
                    a[k] = ww * cos_phi
                    a[k + 1] = ww * sin_phi

        else:
            raise ValueError(f"Invalid component: {component}")
        return a


def associate_legendre_functions_derivative(m, theta, max_degree=MAX_DEGREE):
    """Associated Legendre functions derivative.

    Computes the derivative of the associated Legendre functions for a given
    order up to a maximum degree for a given colatitude (in degree). The
    functions follows the Schmidt semi-normalization convention commonly used
    in geomagnetism. The functions are computed using the recursion formula
    given in [1].

    Parameters
    ----------
    m : int
        Order of the Legendre function.
    theta : float
        Colatitude in radians.
    max_degree : int, optional
        Maximum degree of the Legendre function. Default is 8, as defined in
        the MAX_DEGREE constant in this module.

    Returns
    -------
    ndarray of shape (max_degree + 1,)
        Array containing the derivative of the associated Legendre functions
        for the given order up to the maximum degree, evaluated at the given
        colatitude.

    References
    ----------
    [1] Backus, G., Parker, R. L., & Constable, C. (1996). Foundations of
    geomagnetism. Cambridge University Press.

    See Also
    --------
    associated_legendre_functions : Associated Legendre functions.
    """
    # Initialize the functions array
    dlf = np.zeros(max_degree + 1, dtype=float)

    # Pre-compute the sine and cosine of the colatitude
    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)

    # Compute the functions
    if sin_theta == 0:
        if m == 1:
            dlf[1] = -1
            dlf[2] = -np.sqrt(3)
            for l in range(3, max_degree + 1, 1):
                a = (2 * l - 1) / np.sqrt(l * l - 1)
                b = np.sqrt((l * (l - 2)) / (l * l - 1))
                dlf[l] = a * dlf[l - 1] - b * dlf[l - 2]
    else:
        dlf = associated_legendre_functions(m, theta)
        for l in range(max_degree, m, -1):
            a = np.sqrt((l - m) * (l + m))
            b = l
            dlf[l] = (b * cos_theta * dlf[l] - a * dlf[l - 1]) / sin_theta
        dlf[m] = m * cos_theta * dlf[m] / sin_theta

    return dlf


def associated_legendre_functions(m, theta, max_degree=MAX_DEGREE):
    """Associated Legendre functions.

    Computes the associated Legendre functions for a given order up to a maximum
    degree for a given colatitude (in degree). The functions follows the Schmidt
    semi-normalization convention commonly used in geomagnetism. The functions
    are computed using the recursion formula given in [1].

    Parameters
    ----------
    m : int
        Order of the Legendre function.
    theta : float
        Colatitude in radians.
    max_degree : int, optional
        Maximum degree of the Legendre function. Default is 8, as defined in
        the MAX_DEGREE constant in this module.

    Returns
    -------
    ndarray of shape (max_degree + 1,)
        Array containing the associated Legendre functions for the given order
        up to the maximum degree, evaluated at the given colatitude.

    References
    ----------
    [1] Backus, G., Parker, R. L., & Constable, C. (1996). Foundations of
    geomagnetism. Cambridge University Press.
    """
    # Initialize the functions array
    lf = np.zeros(max_degree + 1, dtype=float)

    # Pre-compute the sine and cosine of the colatitude
    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)

    # Compute the functions
    a = log_factorial(2 * m)
    b = log_factorial(m)
    b = 0.5 * a - b
    b = np.exp(b - m * np.log(2))
    if m != 0:
        b *= np.sqrt(2)
    a = sin_theta
    if a != 0:
        a = np.power(a, m)
    elif m == 0:
        a = 1

    # Compute P(m, m), P(m + 1, m)
    lf[m] = a * b
    if m != max_degree:
        lf[m + 1] = lf[m] * cos_theta * np.sqrt(2 * m + 1)

    # Compute P(m + 2, m), P(m + 3, m), etc.
    for l in range(2, max_degree - m + 1, 1):
        a = (l - 1) * (l + 2 * m - 1)
        b = l * (l + 2 * m)
        c = np.sqrt(a / b)
        a = 2 * (l + m) - 1
        d = a / np.sqrt(b)
        lf[m + l] = d * lf[m + l - 1] * cos_theta - c * lf[m + l - 2]

    return lf


def log_factorial(n):
    """Logarithm of factorial.

    Computes the natural logarithm of the factorial of a given integer,
    given by the following formula:

    .. math::
        \\ln(n!) = \\sum_{i=1}^{n} \\ln(i)


    Parameters
    ----------
    n : int
        Integer.

    Returns
    -------
    float
        Logarithm of the factorial of the given integer.
    """
    return np.sum(np.log(range(1, n + 1, 1)))


def create_map(central_longitude=180):
    """Create a Mollweide map with coastlines.

    Parameters
    ----------
    central_longitude : float, optional
        Central longitude of the map. The default is 180.

    Returns
    -------
    ax : cartopy GeoAxes
        The map axes.
    """
    # Create figure
    plt.figure()

    # Get projection
    projection = ccrs.Mollweide(central_longitude=central_longitude)
    ax = plt.axes(projection=projection)

    # Basic map properties
    ax.set_facecolor("whitesmoke")
    ax.set_global()
    ax.coastlines(lw=0.5)

    return ax


def plot_field(
    colatitude,
    longitude,
    field,
    contour_levels=11,
    central_longitude=180,
    ax=None,
):
    """Plot a field on a Mollweide map.

    Parameters
    ----------
    colatitude : array_like
        Colatitude in degrees.
    longitude : array_like
        Longitude in degrees.
    field : array_like
        Field to plot. Must have the same shape as `colatitude` and `longitude`.
    contour_levels : int, optional
        Number of contour levels. The default is 11.
    central_longitude : float, optional
        Central longitude of the map. The default is 180.
    ax : cartopy GeoAxes, optional
        The map axes. If not provided, a new map is created using `create_map`.

    Returns
    -------
    ax : cartopy GeoAxes
        The map axes.
    """
    # Transform colatitude to latitude
    latitude = 90 - colatitude

    # Recreeate the grid and reshape the field
    longitudes, latitudes = np.meshgrid(longitude, latitude)
    field = np.reshape(field, colatitude.shape + longitude.shape)

    # Define the contour levels
    fmax = np.abs(field).max() + 1e-1
    contour_levels = np.linspace(-fmax, fmax, contour_levels, endpoint=True)

    # Create figure
    ax = ax or create_map(central_longitude=central_longitude)

    # Plot the field
    mappable = ax.contourf(
        longitudes,
        latitudes,
        field,
        levels=contour_levels,
        transform=ccrs.PlateCarree(),
        cmap="RdBu_r",
        extend="both",
    )

    # Add colorbar
    cmap = ax.figure.colorbar(
        mappable,
        orientation="horizontal",
        shrink=0.7,
        pad=0.08,
    )
    cmap.set_label("nT/year")

    return ax


def plot(colatitude, longitude, ax=None, **kwargs):
    """Plot a Mollweide map.

    Parameters
    ----------
    colatitude : array_like
        Colatitude in degrees.
    longitude : array_like
        Longitude in degrees.
    ax : cartopy GeoAxes, optional
        The map axes. If not provided, a new map is created using `create_map`.
    **kwargs: dict, optional
        Keyword arguments passed to `ax.plot`.

    Returns
    -------
    ax : cartopy GeoAxes
        The map axes.
    """
    # Transform colatitude to latitude
    latitude = 90 - colatitude

    # Create figure
    ax = ax or create_map()

    # Plot the field
    ax.plot(
        longitude,
        latitude,
        transform=ccrs.PlateCarree(),
        **kwargs,
    )

    return ax


def scatter(colatitude, longitude, ax=None, **kwargs):
    """Scatter plot on a Mollweide map.

    Parameters
    ----------
    colatitude : array_like
        Colatitude in degrees.
    longitude : array_like
        Longitude in degrees.
    ax : cartopy GeoAxes, optional
        The map axes. If not provided, a new map is created using `create_map`.
    **kwargs: dict, optional
        Keyword arguments passed to `ax.scatter`.

    Returns
    -------
    ax : cartopy GeoAxes
        The map axes.
    """
    # Transform colatitude to latitude
    latitude = 90 - colatitude

    # Create figure
    ax = ax or create_map()

    # Plot the field
    mappable = ax.scatter(
        longitude,
        latitude,
        transform=ccrs.PlateCarree(),
        **kwargs,
    )

    # Add colorbar
    cmap = ax.figure.colorbar(
        mappable,
        orientation="horizontal",
        shrink=0.7,
        pad=0.08,
    )
    cmap.set_label("nT/year")

    return ax


def compute_spectrum(damping_matrix, model, max_degree=MAX_DEGREE):
    """Compute the spectrum of the model.

    Parameters
    ----------
    damping_matrix : ndarray of shape (N_COEFFICIENTS, N_COEFFICIENTS)
        Damping matrix.
    degree : int
        Degree of the spherical harmonic expansion.
    model : ndarray of shape (N_COEFFICIENTS,)
        Model spherical harmonic coefficients.

    Returns
    -------
    mynorm : ndarray of shape (degree + 1,)
        Spectrum of the model.
    """
    # Initialize the spectrum
    spectrum = np.zeros(max_degree + 1)

    # Compute the spectrum
    lm = -1
    for l in range(1, max_degree + 1):
        for m in range(0, l):
            if m == 0:
                lm = lm + 1
                spectrum[l] = (
                    spectrum[l] + damping_matrix[lm, lm] * model[lm] ** 2
                )
            elif m > 0:
                lm = lm + 1
                spectrum[l] = (
                    spectrum[l] + damping_matrix[lm, lm] * model[lm] ** 2
                )
                lm = lm + 1
                spectrum[l] = (
                    spectrum[l] + damping_matrix[lm, lm] * model[lm] ** 2
                )
    return spectrum
