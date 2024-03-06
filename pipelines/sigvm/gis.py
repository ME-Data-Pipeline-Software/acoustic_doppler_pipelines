import numpy as np

earth_radius = 6371000.  # m

"""
Many of the formulas below are from:
www.movable-type.co.uk/scripts/latlong.html

Their phi is latitude, and lambda (lmb) is longitude.
"""


# This was formerly 'poly_length' and had different inputs for
# lon, lat
def haversine(lonlat, units='degrees', radius=earth_radius):
    """Calculate the distance between lon/lat values.
    This function operates along the first dimension of lon/lat.

    Parameters
    ----------
    lonlat : array like (2, N, ...)
    units : str {'degrees' (default), 'radians'}
    radius : float (default: 6371000 [m])

    Returns
    -------
    dist : array like (N-1, ...)
           units match `radius`
    """
    lon, lat = _parse_inputs(lonlat, units)
    return radius * 2 * np.arcsin(np.sqrt(np.sin(np.diff(lat, 1, 0) / 2) ** 2 +
                                          np.cos(lat[:-1]) * np.cos(lat[1:]) *
                                          np.sin(np.diff(lon, 1, 0) / 2) ** 2))


def midpoint(lonlat, units='degrees'):
    """Calculate the midpoint of lon/lat values.
    This function operates along the first dimension of lon/lat.

    Parameters
    ----------
    lonlat : array like (2, N, ...)
    units : str {'degrees' (default), 'radians'}

    Returns
    -------
    lonlat_mid : array like (2, N - 1, ...)

    Notes
    -----
    Units of the returned arrays are consistent with `units`.
    """
    lon, lat = _parse_inputs(lonlat, units)
    dlon = np.diff(lon, 1, 0)
    Bx = np.cos(lat[1:])
    By = Bx * np.sin(dlon)
    Bx *= np.cos(dlon)
    tmp = np.cos(lat[:-1]) + Bx
    out = np.vstack((lon[:-1] + np.arctan2(By, tmp),
                     np.arctan2(np.sin(lat[:-1]) + np.sin(lat[1:]),
                                np.sqrt(tmp ** 2 + By ** 2)), ))
    if units == 'degrees':
        out *= 180 / np.pi
    return out


def bearing(lonlat, units='degrees'):
    """Calculate the bearing from each entry in lonlat to the next.

    This function operates along the first dimension of lon, lat.

    Parameters
    ----------
    lonlat : array like (2, N, ...)

    Returns
    -------
    bearing : array like (N - 1, ...)
              The bearing is in the same units as lon/lat, North is 0
              and increases clockwise.
    """
    lon, lat = _parse_inputs(lonlat, units)
    dlon = np.diff(lon, 1, 0)
    b = np.arctan2(np.sin(dlon) * np.cos(lat[1:]),
                   np.cos(lat[:-1]) * np.sin(lat[1:]) -
                   np.sin(lat[:-1]) * np.cos(lat[1:]) * np.cos(dlon))
    if units == 'degrees':
        b *= 180 / np.pi
    return b


def arg2bearing(val, units='degrees'):
    """
    Calculate the argu
    """
    if units == 'degrees':
        offset = 90
    elif units == 'radians':
        offset = np.pi / 2
    else:
        raise ValueError("Units must be 'degrees' or 'radians'.")
    return offset - val

# This operation is symmetric, so we can just do:
bearing2arg = arg2bearing


def bearing2point(lonlat_start, bearing, distance, units='degrees'):
    """Calculate the position of a new point if you head along a
    great-circle initiated by `bearing` for `distance` meters.

    Parameters
    ----------
    lonlat_start : array like (2, N, ...)
    bearing : array like (N, ...)
    dist : array like (N, ...)
           units: meters

    Returns
    -------
    lonlat_end : array like (2, N, ...)
    """
    lon, lat = _parse_inputs(lonlat_start, units)
    if units == 'degrees':
        bearing = bearing * np.pi / 180  # DO NOT USE *=
    delta = distance / earth_radius
    lat2 = np.arcsin(np.sin(lat) * np.cos(delta) +
                     np.cos(lat) * np.sin(delta) * np.cos(bearing))
    lon2 = lon + np.arctan2(np.sin(bearing) * np.sin(delta) * np.cos(lat),
                            np.cos(delta) - np.sin(lat) * np.sin(lat2))
    if units == 'degrees':
        lon2 *= 180 / np.pi
        lat2 *= 180 / np.pi
    return np.vstack((lon2, lat2))


def diffll(lonlat, units='degrees',
           radius=earth_radius,
           bearing_midpoint=False):
    """Calculate the length of segments, midpoints, and orientations
    between lon/lat values.

    This function operates along the first dimension of lon, lat.

    Parameters
    ----------

    lonlat : array like (2, N, ...)

    units : str {'degrees' (default), 'radians'}

    radius : float (default: 6371000 [m])

    bearing_midpoint : bool (default: True)
      Whether the bearing should be calculated at the midpoint (True),
      or at the origin point (False).

    Returns
    -------
    distC : array like (N - 1, ...)
            complex distance. real: east-component
                              imag: North-component
    lonlatmid : array like (2, N - 1, ...)
    """
    # This forces these variables into radians
    lonlat = np.array(_parse_inputs(lonlat, units))
    lonlatmid = midpoint(lonlat, units='radians')
    dist = haversine(lonlat, units='radians', radius=radius)
    if bearing_midpoint:
        b = bearing(np.hstack((lonlatmid[:, None],
                               lonlat[:, 1:][:, None])),
                    units='radians')[0]
    else:
        b = bearing(lonlat, units='radians')
    if units == 'degrees':
        lonlatmid *= 180 / np.pi
    return dist * np.exp(1j * arg2bearing(b, units='radians')), lonlatmid


def llspace2(lonlat_start, lonlat_stop, npt, units='degrees'):
    """Lat/Lon linspace function along great-circles.

    Parameters
    ----------
    lonlat_start : array like (2 )
    lonlat_stop : array like (2 )
    npt : int
          number of points to create between the two lat/lon points.
    units : str {'degrees' (default), 'radians'}

    Returns
    -------
    lonlat : array like (2, npt)
       lon/lat of points along a great circle between lonlat_start and
       lonlat_stop.

    """
    lonlat = np.hstack((np.array(_parse_inputs(lonlat_start, units)),
                        np.array(_parse_inputs(lonlat_stop, units)),))
    distC, midpoint = diffll(lonlat, units='radians', bearing_midpoint=False)
    lonlat = llspaceD(lonlat[:, :-1], distC,
                      npt, units='radians')
    if units == 'degrees':
        lonlat *= 180 / np.pi
    return lonlat


def llspace(lonlat, npt, units='degrees'):
    """Lat/Lon linspace function along great-circles.

    Parameters
    ----------
    lonlat : array like (2, N) (N>=2)

    npt : int
          number of points to create between successive lat/lon points
          (including endpoints)
    units : str {'degrees' (default), 'radians'}

    Returns
    -------
    lonlat : array like (2, (N - 1) * (npt - 1) + 1)
       lon/lat of points along a great circle between successive
       points in lonlat

    """
    lonlat = _parse_inputs(lonlat, units)
    distC, midpoint = diffll(lonlat, units='radians', bearing_midpoint=False)
    lonlat = llspaceD(lonlat[:, :-1], distC,
                      npt, units='radians')
    if units == 'degrees':
        lonlat *= 180 / np.pi
    return lonlat


class simple_proj(object):
    """
    Create a simple projection where distances from lonlat0 are equal.

    Parameters
    ----------

    lonlat0 : tuple (2)
        The (lon, lat) location of the center of the coordinate system.

    bearing : float
        The angle, in degrees clockwise from True North, that points
        along the x-axis of the coordinate system.

        Alternatively, if the bearing is a complex number, its angle
        specifies the direction of the x-axis counter-clockwise from East.

    radius : float
        The radius of the earth to use. The units of this determine
        the units of the data that is returned from the project. The
        default is 6371000 meters.

    Notes
    -----
    This projection uses a simple haversine formulation to map lat/lon
    points to x/y coordinates. It should not be used for distances
    larger than 100km or so.


    Example Usage
    -------

    # Create a coordinate system centered at 44 N 120 W, with an
    # x-axis that points 5.7 degrees south of East (and a y-axis that
    # is 5.7 degrees east of North)
    myproj = simple_proj((-120, 44), 95.7)

    # Get the x,y points of a set of lon/lat points, in meters
    x, y = myproj([[-121, -120.5, -120.3], [44.2, 44.22, 44.07]])

    """

    def __init__(self, lonlat0, bearing=90.0, radius=earth_radius):
        if np.iscomplex(bearing):
            self._rot = np.conj(bearing / np.abs(bearing))
        else:
            self._rot = np.exp(-1j * np.pi / 180 *
                               bearing2arg(bearing, units='degrees'))
        self.lonlat0 = lonlat0
        self.radius = radius

    def __call__(self, lonlat, units='degrees', inverse=False):
        if not inverse:
            tmp = diffll2(lonlat, self.lonlat0, units=units,
                          radius=self.radius)[..., 0] * self._rot
            return np.vstack((tmp.real, tmp.imag))
        else:
            npt = lonlat.shape[1]
            ll0 = np.tile(self.lonlat0, [npt, 1]).T
            cpt = (lonlat[0] + 1j * lonlat[1]) * np.conj(self._rot)
            b = np.angle(cpt)
            if units == 'degrees':
                b *= 180 / np.pi
            b = arg2bearing(b, units=units)
            dist = np.abs(cpt)
            ll = bearing2point(ll0, b, dist, units=units)
            return ll
            


def llspaceD(lonlat, distc, npt, units='degrees'):
    """Lat/Lon linspace function that uses a distance vector.

    Parameters
    ----------
    lonlat : array like (2, N, ...)
          The origin points.
    distC : array like (N, ...)
          The complex distance between the origin points, and the end
          point (as returned by diffll).
    npt : int or array like (N, ...)
          number of points to create between the lat/lon points
          (including them).
    units : str {'degrees' (default), 'radians'}

    Returns
    -------
    lonlat : array like (2, N * (npt - 1) + 1, ...)
       lon/lat of points along a great circle between lonlat and the
       final point.
    """
    lonlat = _parse_inputs(lonlat, units)
    try:
        len(distc)
    except TypeError:
        distc = np.array([distc])
    d = np.abs(distc)
    b = arg2bearing(np.angle(distc), units='radians')
    # Now calculate the shape of the output array (for arbitrary
    # values of npt).
    if isinstance(npt, int):
        npt = npt * np.ones(lonlat.shape[1], dtype=int)
    npt[npt < 2] = 2
    outshp = list(lonlat.shape)
    outshp[1] = int(np.sum(npt - 1) + 1)
    out = np.empty(outshp, dtype=lonlat.dtype)
    for i in range(len(npt)):
        inow = int(np.sum(npt[:i] - 1))
        iend = int(np.sum(npt[:i + 1] - 1) + 1)
        out[:, inow] = lonlat[:, i]
        out[:, inow + 1:iend] = bearing2point(lonlat[:, i], b[i],
                                              np.linspace(0, d[i], npt[i])[1:],
                                              units='radians')
    if units == 'degrees':
        out *= 180 / np.pi
    return out


# This was formerly 'great_circle_dist', and had separate inputs for
# lon0, lat0, lon1, lat1
def haversine2(lonlat0, lonlat1, units='degrees', radius=earth_radius):
    """Calculate the distance between all of the points lon0,lat0 and
    lon1,lat1.

    Parameters
    ----------
    lonlat0 : array like (2, N, )
    lonlat1 : array like (2, M, )
    units : str {'degrees' (default), 'radians'}
    radius : float (default: 6371000 [m])

    Returns
    -------
    dist : array like (N, M)
           units match `radius`
    """
    # First dim for lat/lon, second for diff, then separate arrays.
    ll1 = np.tile(lonlat1[:, None, None, :],
                  [1, 1, ] + list(lonlat0.shape[1:]) + [1])
    ll0 = np.tile(lonlat0[:, None, :, None],
                  [1, 1, 1, ] + list(lonlat1.shape[1:]))
    # Drop the lat/lon and diff dimensions
    return haversine(np.hstack((ll1, ll0)), units=units, radius=radius)[0]


def diffll2(lonlat0, lonlat1, units='degrees',
            radius=earth_radius, bearing_midpoint=False):
    """Calculate the length of segments, midpoints, and orientations
    between each point in lonlat0 and lonlat1.

    Parameters
    ----------
    lonlat0 : array like (2, N, )
    lonlat1 : array like (2, M, )
    units : str {'degrees' (default), 'radians'}
    radius : float (default: 6371000 [m])

    Returns
    -------
    distC : array like (N, M)
            complex distance. real: east-component
                              imag: North-component
    lonlatmid : array like (2, N, M)
    """
    lonlat0 = _parse_inputs(lonlat0, units=units)
    lonlat1 = _parse_inputs(lonlat1, units=units)
    ll0 = np.tile(lonlat0[:, None, :, None], (1, 1, 1, lonlat1.shape[1]))
    ll1 = np.tile(lonlat1[:, None, None, :], (1, 1, lonlat0.shape[1], 1))
    return diffll(np.hstack((ll1, ll0)),
                  units='radians', radius=radius,
                  bearing_midpoint=bearing_midpoint)[0]


def tovector(lonlat, units='degrees'):
    lon, lat = _parse_inputs(lonlat)
    return np.array((np.cos(lat) * np.cos(lon),
                     np.cos(lat) * np.sin(lon),
                     np.sin(lat)))


def tolonlat(vec, units='degrees'):
    lonlat = np.array((
        np.arctan2(vec[1], vec[0]),
        np.arctan2(vec[2], np.sqrt(vec[0] ** 2 + vec[1] ** 2)),
    ))
    if units == 'degrees':
        lonlat *= 180 / np.pi
    return lonlat


# def intersect_bearing(lonlat0_start, bearing0, lonlat1_start, bearing1):
#     dlat = np.
#     delta12 = 2 * np.arcsin(np.sqrt(np.sin(dlat / 2) ** 2 +
#                                     np.cos()
#     ))

def _within(vals, range):
    out = (range[:-1] <= vals) & (vals <= range[1:])
    out |= (range[1:] <= vals) & (vals <= range[:-1])
    out |= np.isclose(range[1:], vals)
    out |= np.isclose(range[:-1], vals)
    return out


def intersect2(lonlat0, lonlat1, units='degrees', outside=None, ):
    """Find the intersections of the segments lonlat0 and lonlat1.

    Parameters
    ----------
    lonlat0 : array_like (2, N)
              Longitude and latitude values of first line.
    lonlat1 : array_like (2, N)
              Longitude and latitude values of second line.
    units : {'degrees' (default), 'radians'}
    outside : {None, 'keep', scalar, NaN (default)}
        outside specifies the value to apply if the intersection
        is beyond the end of the segment:
          None : Drop all points that are not on both segments.
          'keep' : Keep all intersections of lines.
          NaN/scalar : insert this value where the segments don't
                       overlap.
    """
    vec0 = tovector(lonlat0, units=units)
    vec1 = tovector(lonlat1, units=units)
    c0 = np.cross(vec0[:, 1:], vec0[:, :-1], axis=0)
    c1 = np.cross(vec1[:, 1:], vec1[:, :-1], axis=0)
    vals = np.cross(c0[:, :, None], c1[:, None, :], axis=0)
    mid = (vec0[:, :-1][:, :, None] + vec0[:, 1:][:, :, None] +
           vec1[:, :-1][:, None, :] + vec1[:, 1:][:, None, :])
    flips = (vals * mid).sum(0) < 0
    vals[:, flips] *= -1
    vals = tolonlat(vals, units=units)
    if outside is 'keep':
        return vals
    # Check 0 dim
    good0 = _within(vals[0], lonlat0[0, :, None])
    # Check 1 dim
    good1 = _within(vals[0].swapaxes(1, 0), lonlat1[0, :, None]).swapaxes(1, 0)
    # We have to check Lat too, in case a line follows a meridian.
    good0 &= _within(vals[1],
                     lonlat0[1, :, None])
    good1 &= _within(vals[1].swapaxes(1, 0),
                     lonlat1[1, :, None]).swapaxes(1, 0)
    good = good0 & good1
    if outside is None:
        return vals[:, good]
    vals[:, ~good] = outside
    return vals


def _parse_inputs(lonlat, units='degrees'):
    if not isinstance(lonlat, np.ndarray):
        lonlat = np.array(lonlat)
    if lonlat.ndim == 1 and lonlat.shape[0] == 2:
        lonlat = lonlat.reshape((2, 1))
    elif lonlat.ndim >= 2 and lonlat.shape[0] == 2:
        pass
    else:
        raise ValueError("The shape of the lonlat input are unacceptable.")
    if units == 'degrees':
        lonlat = np.pi / 180 * lonlat
    elif units == 'radians':
        pass
    else:
        raise ValueError("Units must be 'degrees' or 'radians'.")
    return lonlat


class llbox(object):

    def __init__(self, lonlims, latlims):
        self.lonlims = np.sort(lonlims)
        self.latlims = np.sort(latlims)

    @property
    def limdict4basemap(self, ):
        return dict(llcrnrlon=self.lonlims[0],
                    llcrnrlat=self.latlims[0],
                    urcrnrlon=self.lonlims[1],
                    urcrnrlat=self.latlims[1])

    def slice_lon(self, lon, decim=1):
        loninds = np.nonzero((self.lonlims[0] < lon) &
                             (lon < self.lonlims[1]))[0]
        return slice(max(loninds[0] - decim // 2 - 1, 0),
                     max(loninds[-1] + decim, len(lon)),
                     decim)

    def slice_lat(self, lat, decim=1):
        latinds = np.nonzero((self.latlims[0] < lat) &
                             (lat < self.latlims[1]))[0]
        return slice(max(latinds[0] - decim // 2 - 1, 0),
                     min(latinds[-1] + decim, len(lat)),
                     decim)

    @property
    def bboxpath(self, ):
        """The lonlat array that can be used to draw the bounding box."""
        lo = self.lonlims
        la = self.latlims
        return np.array([[lo[0], lo[1], lo[1], lo[0], lo[0]],
                         [la[0], la[0], la[1], la[1], la[0]]])

    @property
    def dx(self, ):
        lat = np.abs(self.latlims).min() * np.ones(2)
        return haversine([self.lonlims, lat])

    @property
    def dy(self, ):
        return haversine([np.zeros(2), self.latlims])

    @property
    def aspect(self, ):
        return self.dx / self.dy

    def inside(self, lonlat):
        lon = lonlat[0]
        lat = lonlat[1]
        return ((self.lonlims[0] < lon) &
                (lon < self.lonlims[1]) &
                (self.latlims[0] < lat) &
                (lat < self.latlims[1]))


def argnearest(lls, ll, dll=None):
    """Find the nearest point in lls to ll.

    Parameters
    ----------
    lls : array_like(2, N)
        A 2xN array of lonlat points.
    ll : array_like(2)
        The point for which you want the nearest index.
    dll : float, array_like(2), or None (default)
        Crop the lls so that they are within +/- this many degrees of
        ll. This is used to accelerate the processing to ignore points
        outside of a certain box size. If `None` then no cropping is
        done (default).

    Notes
    -----
    The size of the crop box is doubled until points are found inside
    the box. Then the algorithm find the nearest one.
    """
    if dll is None:
        inds = range(len(lls[0]))
    else:
        # Check the form of dll
        if not hasattr(dll, '__iter__'):
            dll = [dll, dll]
        if len(dll) != 2:
            raise Exception('`dll` must either be a'
                            ' scalar, or a two-element list/tuple.')
        itmp = np.zeros(len(lls[0]), dtype='bool')
        dll = np.array(dll)
        while not itmp.any():
            bx = llbox([ll[0] - dll[0], ll[0] + dll[0]],
                       [ll[1] - dll[1], ll[1] + dll[1]])
            itmp |= bx.inside(lls)
            dll *= 2
        lls = lls[itmp]
        inds = np.nonzero(itmp)[0]
