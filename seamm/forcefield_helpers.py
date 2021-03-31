from enum import Enum
from seamm_util import Q_

class NonbondForms(Enum):
     SIGMA_EPS = 'sigma-eps'
     RMIN_EPS = 'rmin-eps'
     A_B = 'A-B'
     AR_BR = 'A/r-B/r'
 
def rmin_to_sigma(rmin):
    """Convert rmin to sigma for LJ potential."""
    return rmin / two_raised_to_one_sixth

 
def sigma_to_rmin(sigma):
    """Convert sigma to rmin for LJ potential."""
    return sigma * two_raised_to_one_sixth

 
def nonbond_transformation(
    in_form=NonbondForms.SIGMA_EPS,
    in1_units=None,
    in2_units=None,
    out_form=NonbondForms.SIGMA_EPS,
    out1_units='angstrom',
    out2_units='kcal/mol'
):
    """Return the transform method and unit conversions for the nonbonds.

    Parameters
    ----------
    in_form : enum (NonbondForms)
        The form of the parameters input ('sigma-eps', 'rmin-eps', 'A-B' or
        'A/r-B/r')
    in1_units : string
        The units for the first input parameter
    in2_units : string
        The units for the second input parameter
    out_form : enum (NonbondForms)
        The form of the parameters output.
    out1_units : string
        The units for the first output parameter
    out2_units : string
        The units for the second output parameter

    Returns
    -------
    function : function object
        Python function to transform the parameters
    factor1 : float
        Conversion factor to apply to first transformed parameter
    factor2 : float
        Conversion factor to apply to the second transformed parameter
    """
    if out_form == NonbondForms.SIGMA_EPS:
        if in_form == NonbondForms.SIGMA_EPS:
            transform = no_transform
            factor1 = Q_(1.0, in1_units).to(out1_units).magnitude
            factor2 = Q_(1.0, in2_units).to(out2_units).magnitude
        elif in_form == NonbondForms.RMIN_EPS:
            transform = rmin_eps_to_sigma_eps
            factor1 = Q_(1.0, in1_units).to(out1_units).magnitude
            factor2 = Q_(1.0, in2_units).to(out2_units).magnitude
        elif in_form == NonbondForms.A_B:
            transform = a_b_to_sigma_eps
            A = Q_(1.0, in1_units)
            B = Q_(1.0, in2_units)
            factor1 = (A / B)**(1 / 6).to(out1_units).magnitude
            factor2 = (B**2 / (4 * A)).to(out2_units).magnitude
        elif in_form == NonbondForms.AR_BR:
            transform = ar_br_to_sigma_eps
            A = Q_(1.0, in1_units)**12
            B = Q_(1.0, in2_units)**6
            sigma = (A / B)**(1 / 6)
            eps = B**2 / A
            factor1 = sigma.to(out1_units).magnitude
            factor2 = eps.to(out2_units).magnitude
        else:
            raise ValueError(
                "Cannot handle nonbond input form '" + str(in_form) + "'."
            )
    elif out_form == NonbondForms.RMIN_EPS:
        if in_form == NonbondForms.RMIN_EPS:
            transform = no_transform
            factor1 = Q_(1.0, in1_units).to(out1_units).magnitude
            factor2 = Q_(1.0, in2_units).to(out2_units).magnitude
        else:
            raise ValueError(
                "Cannot handle nonbond input form '" + str(in_form) + "'."
            )
    elif out_form == NonbondForms.A_B:
        raise NotImplementedError(
            "Nonbond output form '" + str(out_form) +
            "' not implemented yet."
        )
    elif out_form == NonbondForms.AR_BR:
        raise NotImplementedError(
            "Nonbond output form '" + str(out_form) +
            "' not implemented yet."
        )
    else:
        raise ValueError(
            "Cannot handle nonbond output form '" + str(out_form) + "'."
        )

    return lambda p1, p2: transform(p1, p2, factor1, factor2)

 
def no_transform(in1, in2, factor1, factor2):
    """No transformation of nonbond parameters, just units."""
    return in1 * factor1, in2 * factor2

 
def rmin_eps_to_sigma_eps(rmin, eps, factor1, factor2):
    """Transform nonbond parameters from rmin-eps to sigma-eps
    and apply the unit conversion factors
    """
    return Forcefield.rmin_to_sigma(rmin) * factor1, eps * factor2

 
def a_b_to_sigma_eps(A, B, factor1, factor2):
    """Transform nonbond parameters from A-B to sigma-eps
    and apply the unit conversion factors
    """
    if A == 0 and B == 0:
        return 0.0, 0.0
    else:
        sigma = (A / B)**(1 / 6)
        eps = B**2 / (4 * A)
        return sigma * factor1, eps * factor2

 
def ar_br_to_sigma_eps(A, B, factor1, factor2):
    """Transform nonbond parameters from A/r-B/r to sigma-eps
    and apply the unit conversion factors
    """
    if A == 0 and B == 0:
        return 0.0, 0.0
    else:
        A = A**12
        B = B**6
        sigma = (A / B)**(1 / 6)
        eps = B**2 / (4 * A)
        return sigma * factor1, eps * factor2
