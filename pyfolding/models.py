#!/usr/bin/env python

"""
Python implementation of common model fitting operations to
analyse protein folding data. Simply automates some fitting
and value calculation. Will be extended to include phi-value
analysis and other common calculations.

Allows for quick model evaluation and plotting.

Also tried to make this somewhat abstract and modular to
enable more interesting calculations, such as Ising models
and such.

Requirements (recommended python 2.7+):
    - numpy
    - scipy
    - matplotlib

Lowe, A.R. 2015

"""

import sys
import inspect
import numpy as np
import scipy as sp

import core
import constants

__author__ = "Alan R. Lowe"
__email__ = "a.lowe@ucl.ac.uk"


def list_models():
    """ List the kinetic of equilibrium models defined in this module.

    Returns a list of the names of the models, whose parent class is
    FitModel.
    """
    clsmembers = inspect.getmembers(sys.modules[__name__], inspect.isclass)
    verif = lambda cls: 'Verified: {0}'.format(cls[1]().verified)
    fit_models = [ (cls[0], verif(cls)) for cls in clsmembers if cls[1].__bases__[0] == core.FitModel ]
    return fit_models








class TemplateModel(core.FitModel):
    """ A template model for expansion
    """
    def __init__(self):
        core.FitModel.__init__(self)
        fit_args = self.fit_func_args
        self.params = tuple( [(fit_args[i],i) for i in xrange(len(fit_args))] )
        self.default_params = np.array([])


    def fit_func(self, x):
        raise NotImplementedError

    @property
    def equation(self):
        return r'F=f(x)'





# F = \frac{\exp( m(x-d_{50})) / RT} { 1+\exp(m(x-d_{50}))/RT}
"""
==========================================================
EQUILIBRIUM FOLDING models
==========================================================
"""

class TwoStateEquilibrium(core.FitModel):
    """ Two state equilibrium denaturation curve - No sloping baseline.

    Folding Scheme:
        N <-> D

    Params:
        F = Fraction unfolded
        m = m-value
        x = denaturant concentration (M)
        d50 = denaturant midpoint (M)
        R = Universal Gas Constant (kcal.mol-1.K-1)
        T = Temperature (Kelvin)

    Reference:
        Clarke and Fersht. Engineered disulfide bonds as probes of
        the folding pathway of barnase: Increasing the stability
        of proteins against the rate of denaturation.
        Biochemistry (1993) vol. 32 (16) pp. 4322-4329
    """
    def __init__(self):
        core.FitModel.__init__(self)
        fit_args = self.fit_func_args
        self.params = tuple( [(fit_args[i],i) for i in xrange(len(fit_args))] )
        self.default_params = np.array([1.5, 5.])
        self.verified = True


    def fit_func(self, x, m, d50):
        F = ( np.exp((m*(x-d50))/core.temperature.RT)) / (1.+np.exp((m*(x-d50))/core.temperature.RT))
        return F

    @property
    def equation(self):
        return r'F = \frac{\exp( m(x-d_{50})) / RT} { 1+\exp(m(x-d_{50}))/RT}'




class TwoStateEquilibriumSloping(core.FitModel):
    """ Two state equilibrium denaturation curve - Sloping baseline.

    Folding Scheme:
        N <-> D

    Params:
        F = Fraction unfolded
        alpha f = intercept of the native baseline at low denaturation concentrations
        beta f = slope/gradient of the native baseline at low denaturation concentrations
        alpha u = intercept of the denatured baseline at high denaturation concentrations
        beta u = slope/gradient of the denatured baseline at high denaturation concentrations
        m = m-value
        x = denaturant concentration (M)
        d50 = denaturant midpoint (M)
        R = Universal Gas Constant (kcal.mol-1.K-1)
        T = Temperature (Kelvin)

    Reference:
        Clarke and Fersht. Engineered disulfide bonds as probes of
        the folding pathway of barnase: Increasing the stability
        of proteins against the rate of denaturation.
        Biochemistry (1993) vol. 32 (16) pp. 4322-4329
    """
    def __init__(self):
        core.FitModel.__init__(self)
        fit_args = self.fit_func_args
        self.params = tuple( [(fit_args[i],i) for i in xrange(len(fit_args))] )
        self.default_params = np.array([1., 0.1, 0.0, 0.1, 1.5, 5.])
        self.verified = True


    def fit_func(self, x, alpha_f, beta_f, alpha_u, beta_u, m, d50):
        F = (alpha_f+beta_f*x) + (alpha_u+beta_u*x) * (\
        ( np.exp((m*(x-d50))/core.temperature.RT)) / (1.+np.exp((m*(x-d50))/core.temperature.RT)))
        return F

    @property
    def equation(self):
        return r'F = (\alpha_f+\beta_f x) + (\alpha_u+\beta_u x) \cdot \frac{\exp( m(x-d_{50})) / RT} { 1+\exp(m(x-d_{50}))/RT}'

# NOTE (ergm) added on 30/8/2017 and corrected incorrect asscii for running on PC 8/9/2017

class ThreeStateEquilibrium (core.FitModel):
    """ Three state equilbrium denaturation curve.

    Folding Scheme:
        N <-> I <-> D

    Params:
        Y_obs = The spectroscopic signal maximum as a function of denaturant concentration
        Y_N = spectroscopic signals of the  native state
        Y_D = spectroscopic signals of the denatured state
        F_D = fraction denatured
        F_N = fraction native
        F_I = fraction intermediate
        Kni = equilibrium contstant of unfolding native to intermediate state
        Kid = equilibrium contstant of unfolding intermediate to denatured state
        DGni = stability of native state relative to intermediate state
        m_ni = m-value of native to intermediate transition
        DGid = stability of intermediate state relative to denatured state
        m_id = m-value of intermediate to denatured transition
        x = denaturant concentration (M)
        R = Universal Gas Constant (kcal.mol-1.K-1)
        T = Temperature (Kelvin)

    Reference:
        Hecky J, Muller K.M. Structural Perturbation and Compensation by Directed
        Evolution at Physiological Temperature Leads to Thermostabilization of
        beta-Lactamase. (2005) Biochemistry 44. pp. 12640-12654
    """
    def __init__(self):
        core.FitModel.__init__(self)
        fit_args = self.fit_func_args
        self.params = tuple( [(fit_args[i],i) for i in xrange(len(fit_args))] )
        self.default_params = np.array([1., 0.5, 0.0, 5., 1.5, 5., 1])
        # NOTE (ergm) added on 3/11/2017
        self.verified = True


    def fit_func(self, x, Y_N, Y_I, Y_D, DGni, m_ni, DGid, m_id):
        F = (Y_N + Y_I*np.exp((-DGni + m_ni*x)/core.temperature.RT) + Y_D*np.exp((-DGni + m_ni*x)/core.temperature.RT) * np.exp((-DGid + m_id*x)/core.temperature.RT)) \
        / (1 + np.exp((-DGni + m_ni*x)/core.temperature.RT) + np.exp((-DGni + m_ni*x)/core.temperature.RT) * np.exp((-DGid + m_id*x)/core.temperature.RT))
        return F

    @property
    def equation(self):
        return r'\begin{equation} \
                \begin{aligned} \
                & \Upsilon_{obs} = \Upsilon_N F_N + \Upsilon_I F_I + \Upsilon_D F_D \ \\ \
                \text{where:}  \\ \
                & F_N = \frac{1} {1 + K_{NI} + K_{NI} K_{ID}}\\ \
                & F_I = \frac{K_{NI}} {1 + K_{NI} + K_{NI} K_{ID}}\\ \
                & F_D = \frac{K_{NI} K_{ID}} {1 + K_{NI} + K_{NI} K_{ID}}\\ \
                \text{and:}  \\ \
                & K_{NI} = \exp \frac{\Delta G_{NI}^{H_2O} + m_{NI} x} {RT}\\ \
                & K_{ID} = \exp \frac{\Delta G_{ID}^{H_2O} + m_{ID} x} {RT}\\ \
                \\ \
                \text{thus:} \\ \
                & \Upsilon_{obs} = \frac{ \Upsilon_N + \Upsilon_I \exp \frac {\Delta G_{NI}^{H_2O} + m_{NI} x} {RT} + \
                \Upsilon_D \exp \frac{\Delta G_{NI}^{H_2O} + m_{NI} x} {RT} \cdot \exp \frac{\Delta G_{ID}^{H_2O} + m_{ID} x} {RT}} {1 + \exp \
                \frac{\Delta G_{NI}^{H_2O} + m_{NI} x} {RT} + \exp \frac{\Delta G_{NI}^{H_2O} + m_{NI} x} {RT} \cdot  \
                 \exp \frac{\Delta G_{ID}^{H_2O} + m_{ID} x} {RT}}\
                \end{aligned}\
                \end{equation}'

# NOTE (ergm) added on 1/8/2017
class TwoStateDimerEquilibrium(core.FitModel):
    """ Two State model for a dimer denaturation Equilibrium - No Intermediate.

    Folding Scheme:
        N2 <-> 2D

    Params:
        Y_obs = spectroscopic signal at a given concentration of urea
        Y_N = spectroscopic signal for native monomeric subunits at a concentration of Pt
        Y_D = spectroscopic signal for denatured monomeric subunits at a concentration of Pt
        alpha_N = intercept of the native baseline at low denaturation concentrations
        beta_N = slope/gradient of the native baseline at low denaturation concentrations
        alpha_D = intercept of the denatured baseline at high denaturation concentrations
        beta_D = slope/gradient of the denatured baseline at high denaturation concentrations
        F_D = fraction of unfolded monomers
        K_U = Equilibrium Constant for Unfolding of dimer.
        Pt = total protein concentration. This variable needs to be set per denaturation curve.
        m = m-value
        x = denaturant concentration (M)
        d50 = denaturant midpoint (M)
        R = Universal Gas Constant (kcal.mol-1.K-1)
        T = Temperature (Kelvin)

    Reference:
        Mallam and Jackson. Folding studies on a knotted protein.
        Journal of Molecular Biology (2005) vol. 346 (5) pp. 1409-1421
    """

    def __init__(self):
        core.FitModel.__init__(self)
        fit_args = self.fit_func_args
        self.params = tuple( [(fit_args[i],i) for i in xrange(len(fit_args))] )
        self.default_params = np.array([1., 0.1, 0.0, 0.1, 1.5, 5., 1e-6])
        self.constants = (('Pt',1e-6),)
        # NOTE (ergm) added on 3/11/2017
        self.verified = True


     # NOTE (ergm) added on 25/8/2017
    def fit_func(self, x, alpha_N, beta_N, alpha_D, beta_D, m, d50, Pt):
        K_U = np.exp(((core.temperature.RT * np.log(Pt))-m*(d50-x)) / core.temperature.RT)
        F_D = (np.sqrt((np.square(K_U) + (8 * K_U * Pt))) - K_U) / (4*Pt)
        Y_0 = ((alpha_N + beta_N*x)*(1-F_D)) +  ((alpha_D + beta_D*x)*(F_D))
        return Y_0


    @property
    def equation(self):
        return r'\begin{equation} \
                \begin{aligned} \
                & \Upsilon_{obs} = \Upsilon_N \cdot (1-F_D)  + \Upsilon_D \cdot F_D \\ \
                \text{where} \\ \
                & \Upsilon_N = \alpha_N+\beta_N x \\ \
                & \Upsilon_D = \alpha_D+\beta_D x \\ \
                & F_D = \frac{\sqrt{((K_U^2 + (8 K_U Pt)) - K_U}}  {4 Pt} \\ \
                & K_U = \exp \frac{(RT \ln(Pt - m(d_{50} - x))} {RT}\
                \end{aligned}\
                \end{equation}'


# NOTE (ergm) added on 1/8/2017
# NOTE (ergm) updated Folding Scheme - was wrong 7/9/2017
class ThreeStateMonoIEquilibrium(core.FitModel):
    """ Three State model for a dimer denaturation Equilibrium - Monomeric intermediate.

    Folding Scheme:
        N2 <-> 2I <-> 2D

    Params:
        Y_rel = spectroscopic signal at a given concentration of urea
        Y_N = spectroscopic signal for native state
        Y_D = spectroscopic signal for denatured state
        Y_I = spectroscopic signal for intermediate state
        F_D = fraction denatured monomers
        F_N = fraction native dimers
        F_I = fraction intermediate dimers
        Pt = total protein concentration. This variable needs to be set per denaturation curve.
        K1 = equilibrium constant of unfolding for native to intermediate state
        K2 = equilibrium constant of unfolding for intermediate to denatured state
        DG1 = stability of native state relative to intermediate state
        m1 = m-value of native to intermediate transition
        DG2 = stability of intermediate state relative to denatured state
        m2 = m-value of intermediate to denatured transition
        x = denaturant concentration (M)
        R = Universal Gas Constant (kcal.mol-1.K-1)
        T = Temperature (Kelvin)

    Reference:
        Mallam and Jackson. Folding studies on a knotted protein.
        Journal of Molecular Biology (2005) vol. 346 (5) pp. 1409-1421
    """
    def __init__(self):
        core.FitModel.__init__(self)
        fit_args = self.fit_func_args
        self.params = tuple( [(fit_args[i],i) for i in xrange(len(fit_args))] )
        self.default_params = np.array([1., 0.1, 1.0, 0.1, 1.5, 5., 3., 1e-6])
        self.constants = (('Pt',1e-6),)
        # NOTE (ergm) added on 3/11/2017
        self.verified = True

    def fit_func(self, x, DG1, m1, DG2, m2, Y_N, Y_I, Y_D, Pt):
        K1 = np.exp((-DG1 + (m1*x)) / core.temperature.RT)
        K2 = np.exp((-DG2 + (m2*x)) / core.temperature.RT)
        F_I = -(K1*(1+K2) + (np.sqrt(np.square(K1) * np.square(1+K2) +(8*Pt*K1)))) / (4*Pt)
        Y_rel = (Y_N * ((2 * Pt * np.square(F_I))/K1)) + (Y_I * F_I) + (Y_D * (K2*F_I))
        return Y_rel


    @property
    def equation(self):
        return r'\begin{equation} \
                \begin{aligned} \
                & \Upsilon_{rel} = \Upsilon_N F_N + \Upsilon_I F_I + \Upsilon_D F_D \\ \
                \text{expanded:} \\ \
                & \Upsilon_{rel} = \Upsilon_N \cdot \frac{2PtF_I^2} {K_1} + \Upsilon_I F_I + \Upsilon_D * K_2F_I \\ \
                \\ \
                \text{where:} \\ \
                & F_I = \frac {- K_1 (1+K_2) + \sqrt{(K_1^2 (1+K_2)^2 + (8 Pt K_1))}} {4Pt} \\ \
                & K_1 = \exp \frac{-\Delta G_{H_20}^1 + m_1 x} {RT} \\ \
                & K_2 = \exp \frac{-\Delta G_{H_20}^2 + m_2 x} {RT}\
                \end{aligned}\
                \end{equation}'

# NOTE (ergm) added on 1/8/2017
# NOTE (ergm) updated Folding Scheme - was wrong 7/9/2017
class ThreeStateDimericIEquilibrium(core.FitModel):
    """ Three State model for a dimer denaturation Equilibrium - Dimeric Intermediate.

    Folding Scheme:
        N2 <-> I2 <-> 2D

    Params:
        Y_rel = spectroscopic signal at a given concentration of urea
        Y_N = spectroscopic signal for native state
        Y_D = spectroscopic signal for denatured state
        Y_I = spectroscopic signal for intermediate state
        F_D = fraction denatured monomers
        F_N = fraction native dimers
        F_I = fraction intermediate dimers
        Pt = total protein concentration. This variable needs to be set per denaturation curve.
        K1 = equilibrium contstant of unfolding native to intermediate state
        K2 = equilibrium contstant of unfolding intermediate to denatured state
        DG1 = stability of native state relative to intermediate state
        m1 = m-value of native to intermediate transition
        DG2 = stability of intermediate state relative to denatured state
        m2 = m-value of intermediate to denatured transition
        x = denaturant concentration (M)
        R = Universal Gas Constant (kcal.mol-1.K-1)
        T = Temperature (Kelvin)

    Reference:
        Mallam and Jackson. Folding studies on a knotted protein.
        Journal of Molecular Biology (2005) vol. 346 (5) pp. 1409-1421
    """

    def __init__(self):
        core.FitModel.__init__(self)
        fit_args = self.fit_func_args
        self.params = tuple( [(fit_args[i],i) for i in xrange(len(fit_args))] )
        self.default_params = np.array([1., 0.1, 0.0, 0.1, 1.5, 5., 2., 1e-6])
        self.constants = (('Pt',1e-6),)
        # NOTE (ergm) added on 3/11/2017
        self.verified = True

    def fit_func(self, x, DG1, m1, DG2, m2, Y_N, Y_I, Y_D, Pt):
        K1 = np.exp((-DG1 + (m1*x)) / core.temperature.RT)
        K2 = np.exp((-DG2 + (m2*x)) / core.temperature.RT)
        F_D = (-(K1*K2) + np.sqrt(np.square(K1*K2) + 8*(1+K1)*(K1*K2)*Pt)) / (4*Pt*(1+K1))
        Y_rel = (Y_N * ((2 * Pt * np.square(F_D))/(K1*K2))) + (Y_I * ((2 * Pt * np.square(F_D))/K2)) + (Y_D * F_D)
        return Y_rel



    @property
    def equation(self):
        return r'\begin{equation} \
                \begin{aligned} \
                & \Upsilon_{rel} = \Upsilon_N F_N + \Upsilon_I F_I + \Upsilon_D F_D \\ \
                \text{expanded:} \\ \
                & \Upsilon_{rel} = \Upsilon_N \cdot \frac{2PtF_D^2} {K_1 K_2} + \Upsilon_I \frac{2PtF_D^2} {K_2} + \Upsilon_D * (F_D) \\ \
                \\ \
                \text{where:} \\ \
                & F_D = \frac {- K_1 K_2 + \sqrt{((K_1 K_2)^2 + 8(1+K_1)(K_1 K_2)Pt)}} {4Pt (1 + K_1)} \\ \
                & K_1 = \exp \frac{-\Delta G_{H_20}^1 + m_1 x} {RT} \\ \
                & K_2 = \exp \frac{-\Delta G_{H_20}^2 + m_2 x} {RT}\
                \end{aligned}\
                \end{equation}'


class HomozipperIsingEquilibrium(core.FitModel):
    """ Homopolymer Zipper Ising model

    Params:
        q = partition function
        f = fraction of folded protein
        Kappa = equilibrium constant of folding for a given repeating unit
        Tau = equilibrium constant of association between 2 repeating units
        n = number of repeating units
        x = denaturant concentration (M)
        Gi = intrinsic stability (folding energy) of a repeating unit i
        mi = denaturant sensitivity of the intrinsic stability of a repeating unit i
        Gi,i+1 = interface interaction energy between 2 repeating units
        R = Universal Gas Constant (kcal.mol-1.K-1)
        T = Temperature (Kelvin)

    Reference:
        Aksel and Barrick. Analysis of repeat-protein folding using
        nearest-neighbor statistical mechanical models.
        Methods in enzymology (2009) vol. 455 pp. 95-125
    """
    def __init__(self):
        core.FitModel.__init__(self)
        fit_args = self.fit_func_args
        self.params = tuple( [(fit_args[i],i) for i in xrange(len(fit_args))] )
        self.default_params = np.array([7, 0.1, -.53, -4.6])
        self.constants = (('n',7),)
        self.verified = True

    def fit_func(self, x, n, DG_intrinsic, m_intrinsic, DG_interface):

        # # clamp to prevent instability
        # if DG_intrinsic<0. or DG_interface>0.:
        #     return core.FIT_ERROR(x)

        k = np.exp(-(DG_intrinsic - m_intrinsic*x) / core.temperature.RT )
        #t = np.exp(-(DG_interface - m_interface*x) / core.temperature.RT )
        t = np.exp(-(DG_interface) / core.temperature.RT )
        pre_factor = (k/(n*(k*t-1)))
        numerator = n*(k*t)**(n+2) - (n+2)*(k*t)**(n+1) + (n+2)*k*t-n
        denominator = (k*t-1)**2 + k*((k*t)**(n+1) - (n+1)*k*t+n )
        theta = pre_factor * (numerator / denominator)
        return 1.-theta

    # NOTE (ergm) changed on 4/9/2017
    @property
    def equation(self):
        return r'\text{the partition function } (q) \text{ and thus fraction of folded protein } (f) \text{ of n arrayed repeats are given by:}\\ \
                \begin{equation} \\ \
                \begin{aligned} \
                & q = 1 + \frac{\kappa([\kappa \tau]^{n+1} - [n+1]\kappa \tau - n)} {(\kappa \tau + 1)^2}  \\ \
                \\ \
                & f = \frac{1} {n} \sum^{n}_{i=0}i\frac{(n-i+1)\kappa^i\tau^{i-1}} {q} \\ \
                \\ \
                \text{where:} \\ \
                & \kappa (x) = \exp\frac{-G_i} {RT} = \exp\frac{-G_{i,H_20} + m_i x} {RT}  \\ \
                \\ \
                & \tau (x) = \exp\frac{-G_{i,i+1}} {RT} \
                \end{aligned}\
                \end{equation}'


class HeteropolymerIsingEquilibrium(core.FitModel):
    """ Heteropolymer Ising model

    Params:
        q = partition function
        f = fraction of folded protein
        Kappa = equilibrium constant of folding for a given repeating unit
        Tau = equilibrium constant of association between 2 repeating units
        n = number of repeating units
        x = denaturant concentration (M)
        DG_intrinsic = intrinsic stability (folding energy) of a repeating unit i
        m_intrinsic = denaturant sensitivity of the intrinsic stability of a repeating unit i
        DG_interface = interface interaction energy between 2 repeating units
        R = Universal Gas Constant (kcal.mol-1.K-1)
        T = Temperature (Kelvin)

    Reference:
        Aksel and Barrick. Analysis of repeat-protein folding using
        nearest-neighbor statistical mechanical models.
        Methods in enzymology (2009) vol. 455 pp. 95-125
    """
    def __init__(self):
        core.FitModel.__init__(self)

    def fit_func(self, x):
        raise NotImplementedError('This is a dummy model.')

    # NOTE (ergm) changed on 4/9/2017
    @property
    def equation(self):
        return r'\text{the partition function } (q) \text{ and thus fraction of folded protein } (f) \text{ of n arrayed repeats are given by:}  \\ \
            \begin{equation} \\ \
            \begin{aligned} \\ \
            \kappa(x) &= \exp(-(\Delta G_{intrinsic} - m_{intrinsic}x) / RT) \\ \
            \tau(x) &= \exp(-\Delta G_{interface}) / RT) \\ \
            q(i) &=  \
            \begin{bmatrix} 0 & 1\end{bmatrix}  \
            \begin{bmatrix} \kappa_1\tau_{-1} & 1\\ \kappa & 1 \end{bmatrix} \
            \ldots  \
            \begin{bmatrix} \kappa_n\tau_{n-1} & 1\\ \kappa & 1 \end{bmatrix} \
            \begin{bmatrix} 1 \\ 1 \end{bmatrix} \\ \
            \theta &=  \frac{1}{nq(n)} \sum_{i=0}^{n}{q(i)} \
            \end{aligned} \
            \end{equation}'








"""
==========================================================
KINETIC FOLDING models
==========================================================
"""

class TwoStateChevron(core.FitModel):
    """ Two state chevron plot.

    Folding Scheme:
        N <-> D

    Params:
        k obs = rate constant of unfolding or refolding at a particular denaturant concentration
        kf = rate constant of refolding at a particular denaturant concentration
        mf = the gradient of refolding arm of the chevron
        ku = rate constant of unfolding at a a particular denaturant concentration
        mu = the gradient of unfolding arm of the chevron
        x = denaturant concentration (M)

    Reference:
        Jackson SE and Fersht AR.  Folding of chymotrypsin inhibitor 2.
        1. Evidence for a two-state transition.
        Biochemistry (1991) 30(43):10428-10435.
    """
    def __init__(self):
        core.FitModel.__init__(self)
        fit_args = self.fit_func_args
        self.params = tuple( [(fit_args[i],i) for i in xrange(len(fit_args))] )
        self.default_params = np.array([100., 1.3480, 5e-4, 1.])
        #self.constants = (('mf',1.76408),('mu',1.13725))
        self.verified = True


    def fit_func(self, x, kf, mf, ku, mu):
        k_obs = kf*np.exp(-mf*x) + ku*np.exp(mu*x)
        return k_obs

    def error_func(self, y):
        return np.log(y)

    # NOTE (ergm) added on 24/8/2017
    # def components(self, x, kf, mf, ku, mu):
    #     k_f = kf*np.exp(-mf*x)
    #     k_u = ku*np.exp(mu*x)
    #     k_obs = k_f + k_u
    #     return {'k_f':k_f, 'k_u':k_u}

    @property
    def equation(self):
        return r'\begin{equation} \
                \begin{aligned} \
                & k_{obs} = k_f + k_u \\ \
                \\ \
                \text{where:} \\ \
                & k_f = k_f^{H_2O}\exp(-m_{kf}x)\\ \
                & k_u = k_u^{H_2O}\exp(m_{ku}x) \\ \
                \text{thus:} \\ \
                & k_{obs} = k_f^{H_2O}\exp(-m_{kf}x) + k_u^{H_2O}\exp(m_{ku}x)\\ \
                \end{aligned} \
                \end{equation}'


class ThreeStateChevron(core.FitModel):
    """ Three state chevron with single intermediate.

    Folding Scheme:
        N <-> I <-> D

    Params:
        k obs = rate constant of unfolding or refolding at a particular denaturant concentration
        kfi = microscopic rate constant for the conversion of folded to intermediate
        kif = microscopic rate constant for the conversion of intermediate to folded
        i.e. k_if = kif(H20) * exp((mi - mif)*x)
        Kiu = equilibrium constant for the rapid equilibration between intermediate & unfolded
        i.e. Kiu = Kiu(H2O) * exp((mu-mi)*x)
        mif = m-value associated with the kinetic transition between intermediate & folded
        mi = m-value associated with the equilibrium transition between intermediate & folded
        mu = m-value associated with the equilibrium transition between unfolded & folded
        x = denaturant concentration (M)

    Reference:
        Parker et al. An integrated kinetic analysis of
        intermediates and transition states in protein folding reactions.
        Journal of molecular biology (1995) vol. 253 (5) pp. 771-86
    """
    def __init__(self):
        core.FitModel.__init__(self)
        fit_args = self.fit_func_args
        self.params = tuple( [(fit_args[i],i) for i in xrange(len(fit_args))] )
        self.default_params = np.array([4.5e-4, -9.5e-1, 1.3e9, -6.9,  1.4e-8, -1.6])
        #self.constants = (('mif',-0.97996),('mi',-6.00355),('mu',-1.66154))
        self.verified = True

    def fit_func(self, x, kfi, mif, kif, mi, Kiu, mu):
        k_fi = kfi*np.exp(-mif*x)
        k_if = kif*np.exp((mi - mif)*x)
        K_iu = Kiu*np.exp((mu - mi)*x)
        k_obs = k_fi + k_if / (1.+1./K_iu)
        return k_obs

    def error_func(self, y):
        return np.log(y)

    def components(self, x, kfi, mif, kif, mi, Kiu, mu):
        k_fi = kfi*np.exp(-mif*x)
        k_if = kif*np.exp((mi - mif)*x)
        k_obs_I = k_fi + k_if
        return {'kobs_I':k_obs_I}

    @property
    def equation(self):
        return r'\begin{equation} \
                \begin{aligned} \
                & k_{obs} = \frac{k_{fi} + k_{if}} {(1+1/K_{iu})} \\ \
                \\ \
                \text{where:} \\ \
                & k_{fi} = k_{fi}^{H_2O}\exp(-m_{fi}x)\\ \
                & k_{if} = k_{if}^{H_2O}\exp((m_i - m_{if})x)\\ \
                & K_{iu} = K_{iu}^{H_2O}\exp((m_u - m_i)x)\\ \
                \text{thus:} \\ \
                & k_{obs} = k_{fi}^{H_2O}\exp(-m_{if}x) + k_{if}^{H_2O}\exp((m_i - m_{if})x) /(1 + 1 / (K_{iu}^{H_2O}\exp((m_u-m_i)x)))\\ \
                \end{aligned} \
                \end{equation}'



class ThreeStateFastPhaseChevron(core.FitModel):
    """ Three state chevron with single intermediate.

    Folding Scheme: N <-> I <-> D

    Params:
        k obs = rate constant of unfolding or refolding at a particular denaturant concentration
        kfi = microscopic rate constant for the conversion of folded to intermediate
        kif = microscopic rate constant for the conversion of intermediate to folded
        kiu = microscopic rate constant for the conversion of intermediate to unfolded
        kui = microscopic rate constant for the conversion of unfolded to intermediate
        Kiu = equilibrium constant for the rapid equilibration between intermediate & unfolded
        mfi = m-value associated with the kinetic transition between folded & intermediate
        mif = m-value associated with the kinetic transition between intermediate & folded
        miu = m-value associated with the kinetic transition between intermediate & unfolded
        mui = m-value associated with the kinetic transition between unfolded & intermediate
        x = denaturant concentration (M)

    Reference:
        Parker et al. An integrated kinetic analysis of
        intermediates and transition states in protein folding reactions.
        Journal of molecular biology (1995) vol. 253 (5) pp. 771-86
    """

    def __init__(self):
        core.FitModel.__init__(self)
        fit_args = self.fit_func_args
        self.params = tuple( [(fit_args[i],i) for i in xrange(len(fit_args))] )
        self.default_params = np.array([172., 1.42, .445, .641, 1e4, 2.71313, 1.83e-3, 1.06])
        #self.constants = (('kui',172.), ('mui',1.42), ('kiu',.445), ('miu',.641), ('mif',-2.71313),('mfi',1.06534))
        self.verified = True

    def fit_func(self, x, kui, mui, kiu, miu, kif, mif, kfi, mfi):
        k_iu = kiu*np.exp(miu*x)
        k_ui = kui*np.exp(-mui*x)
        k_if = kif*np.exp(-mif*x)
        k_fi = kfi*np.exp(mfi*x)
        K_iu = k_iu / (k_iu+k_ui)
        k_obs = k_fi + k_if / (1.+1./K_iu)
        return k_obs

    def error_func(self, y):
        return np.log(y)

    def components(self, x, kui, mui, kiu, miu, kif, mif, kfi, mfi):
        k_iu = kiu*np.exp(miu*x)
        k_ui = kui*np.exp(-mui*x)
        k_if = kif*np.exp(-mif*x)
        k_fi = kfi*np.exp(mfi*x)
        k_obs_I = k_iu + k_ui
        k_obs_N = k_fi + k_if
        return {'kobs_I':k_obs_I} #, 'kobs_N':k_obs_N}

    # NOTE (ergm) added on 23/8/2017
    @property
    def equation(self):
        return r'\begin{equation} \
                \begin{aligned} \
                & k_{obs} = \frac{k_{fi} + k_{if}} {(1+1/K_{iu})} \\ \
                \\ \
                \text{where:} \\ \
                & k_{fi} = k_{fi}^{H_2O}\exp(m_{fi}x)\\ \
                & k_{if} = k_{if}^{H_2O}\exp(-m_{if}x)\\ \
                & k_{iu} = k_{iu}^{H_2O}\exp(m_{iu}x)\\ \
                & k_{ui} = k_{ui}^{H_2O}\exp(-m_{ui}x)\\ \
                & K_{iu} = \frac{k_{iu}} {k_{iu} + k_{ui}}\\ \
                \end{aligned} \
                \end{equation}'


class ThreeStateSequentialChevron(core.FitModel):
    """ Three state metastable intermediate chevron plot.

    Folding Scheme: N <-> I <-> D

    Params:
        k obs = rate constant of unfolding or refolding at a particular denaturant concentration
        kfi = microscopic rate constant for the conversion of folded to intermediate
        kif = microscopic rate constant for the conversion of intermediate to folded
        kiu = microscopic rate constant for the conversion of intermediate to unfolded
        kui = microscopic rate constant for the conversion of unfolded to intermediate
        mfi = m-value associated with the kinetic transition between folded & intermediate
        mif = m-value associated with the kinetic transition between intermediate & folded
        miu = m-value associated with the kinetic transition between intermediate & unfolded
        mui = m-value associated with the kinetic transition between unfolded & intermediate
        x = denaturant concentration (M)

    Reference:
        Bachmann and Kiefhaber. Apparent two-state tendamistat
        folding is a sequential process along a defined route.
        J Mol Biol (2001) vol. 306 (2) pp. 375-386
    """

    def __init__(self):
        core.FitModel.__init__(self)
        fit_args = self.fit_func_args
        self.params = tuple( [(fit_args[i],i) for i in xrange(len(fit_args))] )
        self.default_params = np.array([2e4, 0.3480, 1e4, 0, 20.163, 1.327, 0.3033, 0.2431])
        # NOTE (ergm) changed constants on 3/10/2017
        self.constants = (('kiu', 1.e4),('miu',0.))
        self.verified = True

    def fit_func(self, x, kui, mui, kiu, miu, kif, mif, kfi, mfi):
        k_ui = kui*np.exp(-mui*x)
        k_iu = kiu*np.exp(miu*x)
        k_if = kif*np.exp(-mif*x)
        k_fi = kfi*np.exp(mfi*x)
        lam_1 = -(k_ui + k_iu + k_if + k_fi)
        lam_2 = k_ui * (k_if+k_fi) + k_iu*k_fi
        k_obs = 0.5 * (-lam_1 - np.sqrt(lam_1**2 - 4*lam_2))
        return k_obs

    def error_func(self, y):
        return np.log(y)

    def components(self, x, kui, mui, kiu, miu, kif, mif, kfi, mfi):
        k_ui = kui*np.exp(-mui*x)
        k_iu = kiu*np.exp(miu*x)
        k_if = kif*np.exp(-mif*x)
        k_fi = kfi*np.exp(mfi*x)
        k_TS1 = k_ui + (k_fi/kif)*k_iu
        k_TS2 = (k_ui/k_iu)*k_if + k_fi
        return {'kTS1':k_TS1, 'kTS2':k_TS2}

    @property
    def equation(self):
        return r'\begin{equation} \
                \begin{aligned} \
                & k_{obs} = 0.5(-A_2 \pm \sqrt{A_2^2 - 4A_1}) \\ \
                \\ \
                \text{where:}\\ \
                & A_1 = -(k_{ui} + k_{iu} + k_{if} + k_{fi}) \\ \
                & A_2 = k_{ui}(k_{if} + k_{fi}) + k_{iu}k_{if} \\ \
                \text{and:} \\ \
                & k_{fi} = k_{fi}^{H_2O}\exp(m_{fi}x)\\ \
                & k_{if} = k_{if}^{H_2O}\exp(-m_{if}x)\\ \
                & k_{iu} = k_{iu}^{H_2O}\exp(m_{iu}x)\\ \
                & k_{ui} = k_{ui}^{H_2O}\exp(-m_{ui}x)\\ \
                \end{aligned} \
                \end{equation}'


class ParallelTwoStateChevron(core.FitModel):
    """ Parallel Two state chevron plot.

    Folding Scheme:
        N <-> D
        ^     ^
        |_____|

    Params:
        k obs = rate constant of unfolding or refolding at a particular denaturant concentration
        k_obs_A = rate constant of unfolding or refolding of pathway A at a particular denaturant concentration
        k_obs_B = rate constant of unfolding or refolding  of pathway B at a particular denaturant concentration
        mf_A = the gradient of refolding arm of pathway A
        mf_B = the gradient of refolding arm of pathway B
        mu_A = the gradient of unfolding arm of pathway A
        mu_B = the gradient of unfolding arm of pathway B
        x = denaturant concentration (M)

    Reference:
        Lowe & Itzhaki. Rational redesign of the folding pathway of a modular protein.
        PNAS (2007) vol. 104 (8) pp. 2679-2684
    """
    def __init__(self):
        core.FitModel.__init__(self)
        fit_args = self.fit_func_args
        self.params = tuple( [(fit_args[i],i) for i in xrange(len(fit_args))] )
        self.default_params = np.array([50., 1.3480, 5e-4, 1., 150., 3.5])


    def fit_func(self, x, kf_A, mf_A, ku_A, mu_A, kf_B, mf_B):

        if mf_A < 0. or mf_B < 0. or mu_A < 0.:
            return core.FIT_ERROR(x)

        if kf_A <0. or ku_A <0. or kf_B < 0.:
            return core.FIT_ERROR(x)

        deltaG_A = kf_A / ku_A
        ku_B = kf_B / deltaG_A
        mu_B = np.abs(mf_A + mu_A) - np.abs(mf_B)
        k_obs_A = kf_A*np.exp(-mf_A*x) + ku_A*np.exp(mu_A*x)
        k_obs_B = kf_B*np.exp(-mf_B*x) + ku_B*np.exp(mu_B*x)
        k_obs = k_obs_A + k_obs_B
        return k_obs

    def error_func(self, y):
        return np.log(y)

    def components(self, x, kf_A, mf_A, ku_A, mu_A, kf_B, mf_B):
        deltaG_A = kf_A / ku_A
        ku_B = kf_B / deltaG_A
        mu_B = np.abs(mf_A + mu_A) - np.abs(mf_B)
        k_obs_A = kf_A*np.exp(-mf_A*x) + ku_A*np.exp(mu_A*x)
        k_obs_B = kf_B*np.exp(-mf_B*x) + ku_B*np.exp(mu_B*x)
        k_obs = k_obs_A + k_obs_B
        return {'kobs_A':k_obs_A, 'kobs_B':k_obs_B}

    # NOTE (ergm) added on 23/8/2017
    @property
    def equation(self):
        return r'\begin{equation} \
                \begin{aligned} \
                & k_{obs} = k_{obs}^A + k_{obs}^B \\ \
                \\ \
                \text{where:}\\ \
                & \Delta G^A = k_f^A / k_u^A \\ \
                & k_u^B = k_f^B / \Delta G^A \\ \
                & m_u^B = (m_f^A + m_u^A) - (m_f^B) \\ \
                & k_{obs}^A = k_f^A exp(-m_f^A x) + k_u^A exp(m_u^A x) \\ \
                & k_{obs}^B = k_f^B exp(-m_f^B x) + k_u^B exp(m_u^B x) \\ \
                \end{aligned} \
                \end{equation}'




class ParallelTwoStateUnfoldingChevron(core.FitModel):
    """ Parallel Two state unfolding chevron plot.

    Folding Scheme:
        N -> D
        |    ^
        |____|

    Params:
        k obs = rate constant of unfolding at a particular denaturant concentration
        k_obs_A = rate constant of unfolding of pathway A at a particular denaturant concentration
        k_obs_B = rate constant of unfolding of pathway B at a particular denaturant concentration
        mu_A = the gradient of unfolding arm of pathway A
        mu_B = the gradient of unfolding arm of pathway B
        x = denaturant concentration (M)

    Reference:
        Hutton et al. Mapping the Topography of a Protein Energy Landscape.
        JACS (2015) vol. 137 (46) pp. 14610-14625
    """
    def __init__(self):
        core.FitModel.__init__(self)
        fit_args = self.fit_func_args
        self.params = tuple( [(fit_args[i],i) for i in xrange(len(fit_args))] )
        self.default_params = np.array([5e-4, 1., 1e-5, 1.5])


    def fit_func(self, x, ku_A, mu_A, ku_B, mu_B):

        if mu_A < 0. or mu_B < 0.:
            return core.FIT_ERROR(x)

        k_obs_A = ku_A*np.exp(mu_A*x)
        k_obs_B = ku_B*np.exp(mu_B*x)
        k_obs = k_obs_A + k_obs_B
        return k_obs

    def error_func(self, y):
        return np.log(y)

    def components(self, x, ku_A, mu_A, ku_B, mu_B):
        k_obs_A = ku_A*np.exp(mu_A*x)
        k_obs_B = ku_B*np.exp(mu_B*x)
        k_obs = k_obs_A + k_obs_B
        return {'kobs_A':k_obs_A, 'kobs_B':k_obs_B}

    # NOTE (ergm) added on 23/8/2017
    @property
    def equation(self):
        return r'\begin{equation} \
                \begin{aligned} \
                & k_{obs} = k_obs^A + k_obs^B \\ \
                \\ \
                \text{where:}\\ \
                & k_obs^A = k_u^A exp(m_u^A x) \\ \
                & k_obs^B = k_u^B exp(m_u^B x) \\ \
                \end{aligned} \
                \end{equation}'



class TwoStateChevronMovingTransition(core.FitModel):
    """ Two state chevron with moving transition state.

    Folding Scheme:
        N <-> D

    Params:
        k obs = rate of unfolding or refolding at a particular denaturant concentration
        kf = rate constant of refolding at a particular denaturant concentration
        mf = refolding coefficient for the first order [D] term.
        ku = rate constant of unfolding at a particular denaturant concentration
        mu = unfolding coefficient for the first order [D] term.
        m' = coefficient for the second-order [D] term (both unfolding and refolding).
        x = denaturant concentration (M)

        Reference:
        Ternstrom et al. From snapshot to movie: phi analysis
        of protein folding transition states taken one step
        further. PNAS (1999) vol. 96 (26) pp. 14854-9

"""
    def __init__(self):
        core.FitModel.__init__(self)
        fit_args = self.fit_func_args
        self.params = tuple( [(fit_args[i],i) for i in xrange(len(fit_args))] )
        # NOTE (ergm) changed on 23/8/2017
        self.default_params = np.array([5e-5, 0.2, 10., 0.2, -1.])
        # NOTE (ergm) added on 3/11/2017
        self.verified = True

    # NOTE (ergm) changed on 23/8/2017
    def fit_func(self, x, ku, mu, kf, mf, m_prime):
        k_obs = ku*(np.exp(mu*x))*(np.exp(m_prime*x*x)) + kf*(np.exp(mf*x))*(np.exp(m_prime*x*x))
        return k_obs

    def error_func(self, y):
        return np.log(y)

    # NOTE (ergm) added on 23/8/2017
    @property
    def equation(self):
        return r'\begin{equation} \
                \begin{aligned} \
                & k_{obs} = k_u + k_f \\ \
                \\ \
                \text{where:}\\ \
                & k_u = k_u^{H_2O} \cdot \exp(m_{u} x) \cdot \exp(m^{*} x^2) \\ \
                & k_f = k_f^{H_2O} \cdot \exp(m_{f} x) \cdot \exp(m^{*} x^2) \\ \
                \end{aligned} \
                \end{equation}'



# NOTE (ergm) added on 24/8/2017 & modified on 7/11/2017
class ChevronPolynomialFit(core.FitModel):
    """ Chevron fit with 2 different second order polynomials for kf & ku.

    Folding Scheme:
        N <-> D

    Params:
        k obs = rate of unfolding or refolding at a particular denaturant concentration
        kf = rate constant of refolding at a particular denaturant concentration
        mf & mf* = are the refolding coefficients for the first and second-order [D] terms, respectively.
        ku = rate constant of unfolding at a particular denaturant concentration
        mu & mu* = are the unfolding coefficients for the first and second-order [D] terms, respectively.
        x = denaturant concentration (M)

    Reference:
    Modified version of equation found in:
    Ternstrom et al. From snapshot to movie: phi analysis
    of protein folding transition states taken one step
    further. PNAS (1999) vol. 96 (26) pp. 14854-9

    """
    def __init__(self):
        core.FitModel.__init__(self)
        fit_args = self.fit_func_args
        self.params = tuple( [(fit_args[i],i) for i in xrange(len(fit_args))] )
        self.default_params = np.array([5e-5, 1., -0.5, 100., 1., -0.5])
        # NOTE (ergm) changed on 3/11/2017
        self.verified = True

    def fit_func(self, x, ku, mu, mu_prime, kf, mf, mf_prime):
        k_obs = ku*(np.exp(mu*x))*(np.exp(mu_prime*x*x)) + kf*(np.exp(mf*x))*(np.exp(mf_prime*x*x))
        return k_obs

    def error_func(self, y):
        return np.log(y)

    @property
    def equation(self):
        return r'\begin{equation} \
                \begin{aligned} \
                & k_{obs} = k_u + k_f\\ \
                \\ \
                \text{where:}\\ \
                & k_u = k_u^{H_2O} \cdot \exp(m_{u} x) \cdot \exp(m_{u}^{*} x^2) \\ \
                & k_f = k_f^{H_2O} \cdot \exp(m_{f} x) \cdot \exp(m_{f}^{*} x^2) \\ \
                \end{aligned} \
                \end{equation}'


if __name__ == "__main__":
    get_models()
