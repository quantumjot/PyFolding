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
	fit_models = [ cls[0] for cls in clsmembers if cls[1].__bases__[0] == core.FitModel ]
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






"""
==========================================================
EQUILIBRIUM FOLDING models
==========================================================
"""

class TwoStateEquilibrium(core.FitModel):
	""" Two state equilbrium denaturation curve.

	F = \frac{\exp( m(x-d_{50})) / RT} { 1+\exp(m(x-d_{50}))/RT}

	Notes:
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


	def fit_func(self, x, m, d50):
		F = ( np.exp((m*(x-d50))/constants.RT)) / (1.+np.exp np.exp((m*(x-d50))/constants.RT))
		return F

	@property
	def equation(self):
		return r'F = \frac{\exp( m(x-d_{50})) / RT} { 1+\exp(m(x-d_{50}))/RT}'




class TwoStateEquilibriumSloping(core.FitModel):
	""" Two state equilbrium denaturation curve.

	F = (\alpha_f+\beta_f x) + (\alpha_u+\beta_u x) \cdot \frac{\exp( m(x-d_{50})) / RT} { 1+\exp(m(x-d_{50}))/RT}

	Notes:
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


	def fit_func(self, x, alpha_f, beta_f, alpha_u, beta_u, m, d50):
		F = (alpha_f+beta_f*x) + (alpha_u+beta_u*x) * (\
		( np.exp((m*(x-d50))/constants.RT)) / (1.+np.exp((m*(x-d50))/constants.RT)))
		return F

	@property
	def equation(self):
		return r'F = (\alpha_f+\beta_f x) + (\alpha_u+\beta_u x) \cdot \frac{\exp( m(x-d_{50})) / RT} { 1+\exp(m(x-d_{50}))/RT}'




class TwoStateDimerEquilibrium(core.FitModel):
	""" Two State model for a dimer denaturation Equilibrium.
	i.e.  N2 = 2D

	Y_0 = Y_N * (1 - F_D) + Y_D * F_D
	Y_N = alpha_N + beta_N * x
	Y_D = alpha_D + beta_D * x
	F_D = (((K_U)^2 + (8 * K_U * P_t - K_U)^0.5)/(4 * P_t)
	K_U = exp((RT * ln(P_t)-m(d_{50}-x)))/RT)

	Notes:
		Mallam and Jackson. Folding studies on a knotted protein. 
		Journal of Molecular Biology (2005) vol. 346 (5) pp. 1409-1421

	Comments:
		Pt is a variable that needs to be set, So it like Ising model when you need 
		to define a specific value to a curve.	
		Could this be the self-constants code line ?

	"""
	
	def __init__(self):
		core.FitModel.__init__(self)
		fit_args = self.fit_func_args
		self.params = tuple( [(fit_args[i],i) for i in xrange(len(fit_args))] )
		self.default_params = np.array([1., 0.1, 0.0, 0.1, 1.5, 5.])
		self.constants = (('Pt',7),)

	def fit_func(self, x, alpha_N, beta_N, alpha_D, beta_D, m, d50):
		K_U = np.exp(((constants.RT * np.ln(Pt))-(m*d50)-x) / constants.RT)
		F_D = ((K_U^2 + (8 * K_U * Pt))^0.5 - K_U) / (4*Pt)
		Y_0 = ((alpha_N + beta_N*x)*(1-F_D)) +  ((alpha_D + beta_D*x)*(F_D))
		return Y_0

	"""
	@property
	def equation(self):
		return r'Y_0 = (\alpha_N+\beta_N x) \cdot (1-F_D) + Y_D \cdot F_D \\ 
				\text{where}\\ 
				F_D = \frac{(K_U^2 + (8 * K_U * Pt))^0.5 - K_U / 4*Pt} (\alpha_u+\beta_u x) \\ 
				K_U = \frac{\exp(((constants.RT * np.ln(Pt))-(m*d50)-x) / constants.RT}'
	"""



class ThreeStateMonoIEquilibrium(core.FitModel):
	""" Three State model for a dimer denaturation Equilibrium.
	i.e.  N2 = I2 = 2D

	Y_rel = (Y_N * ((2*Pt*F_D^2)/(K1*K2))) + (Y_I * ((2*Pt*F_D^2)/K2)) + (Y_D * F_D)
	F_D = -((K1*K2) + ((K1*K2)^2 + (8*(1+K1)*(K1*K2)*Pt))^0.5) / (4*Pt*(1+K1))
	K1 = exp((DG1 + m1*x)/RT)
	K2 = exp((DG2 + m2*x)/RT)

	Notes:
		Mallam and Jackson. Folding studies on a knotted protein. 
		Journal of Molecular Biology (2005) vol. 346 (5) pp. 1409-1421

	Comments:
		Pt is a variable that needs to be set, So it like Ising model when you need 
		to define a specific value to a curve.	
		also needs to fit to multiple datasets 
	
	"""
	
	def __init__(self):
		core.FitModel.__init__(self)
		fit_args = self.fit_func_args
		self.params = tuple( [(fit_args[i],i) for i in xrange(len(fit_args))] )
		self.default_params = np.array([1., 0.1, 0.0, 0.1, 1.5, 5., 2.])
		self.constants = (('Pt',7),)

	def fit_func(self, x, DG1, m1, DG2, m2, Y_N, Y_I, Y_D, Pt):
		K1 = np.exp((-DG1 + (m1*x)) / constants.RT)
		K2 = np.exp((-DG2 + (m2*x)) / constants.RT)
		F_D = -((K1*K2) + ((K1*K2)^2 + (8*(1+K1)*(K1*K2)*Pt))^0.5) / (4*Pt*(1+K1))
		Y_rel = (Y_N * ((2*Pt*F_D^2)/(K1*K2))) + (Y_I * ((2*Pt*F_D^2)/K2)) + (Y_D * F_D)
		return Y_rel

	
	"""
	@property
	def equation(self):
		return r'Y_rel = (Y_N \cdot \frac{(2PtF_D^2/K1K2)} + (Y_I \cdot \frac{(2PtF_D^2)/K2)} + (Y_D * F_D) \\ 
			\text{where}\\ F_D = \frac {- K1K2 + ((K1K2)^2 + (8(1+K1)(K1K2)Pt))^0.5) / 4Pt(1+K1)} \\ 
			K1 = \frac {\exp((-DG1 + (m1 x)) / RT)} \\ K2 = \frac {\exp((-DG2 + (m2 x)) / RT)} '
	"""

class ThreeStateDimericIEquilibrium(core.FitModel):
	""" Three State model for a dimer denaturation Equilibrium.
	i.e.  N2 = 2I = 2D

	Y_rel = (Y_N * ((2*Pt*F_D^2)/(K1*K2))) + (Y_I * ((2*Pt*F_D^2)/K2)) + (Y_D * F_D)
	F_D = -((K1*K2) + ((K1*K2)^2 + (8*(1+K1)*(K1*K2)*Pt))^0.5) / (4*Pt*(1+K1))
	K1 = exp((DG1 + m1*x)/RT)
	K2 = exp((DG2 + m2*x)/RT)

	Notes:
		Mallam and Jackson. Folding studies on a knotted protein. 
		Journal of Molecular Biology (2005) vol. 346 (5) pp. 1409-1421

	Comments:
		Pt is a variable that needs to be set, So it like Ising model when you need 
		to define a specific value to a curve.	
		also needs to fit to multiple datasets.
		
	"""
	
	def __init__(self):
		core.FitModel.__init__(self)
		fit_args = self.fit_func_args
		self.params = tuple( [(fit_args[i],i) for i in xrange(len(fit_args))] )
		self.default_params = np.array([1., 0.1, 0.0, 0.1, 1.5, 5., 3.])
		self.constants = (('Pt',7),)

	def fit_func(self, x, DG1, m1, DG2, m2, Y_N, Y_I, Y_D, Pt):
		K1 = np.exp((-DG1 + (m1*x)) / constants.RT)
		K2 = np.exp((-DG2 + (m2*x)) / constants.RT)
		F_I = -(K1*(1+K2) + ((K1^2*(1+K2)^2 +(8*Pt*K1))^0.5)) / (4*Pt)
		Y_rel = (Y_N * ((2*Pt*F_I^2)/K1)) + (Y_I * F_I) + (Y_D * (K2*F_I))
		return Y_rel

	
	"""
	@property
	def equation(self):
		return r'Y_rel = (Y_N \cdot \frac{(2PtF_I^2/K1)} + (Y_I F_I} + (Y_D * (K2F_I)) \\ 
			\text{where} \\ 
			F_I = \frac {- K1(1+K2) + (K1^2 \cdot(1+K2)^2 + (8 Pt K1))^0.5) / 4Pt} \\ 
			K1 = \frac {\exp((-DG1 + (m1 x)) / RT)} \\ 
			K2 = \frac {\exp((-DG2 + (m2 x)) / RT)}'
	"""



class HomozipperIsingEquilibrium(core.FitModel):
	""" Homopolymer Zipper Ising model

	Notes:
		Aksel and Barrick. Analysis of repeat-protein folding using 
		nearest-neighbor statistical mechanical models. 
		Methods in enzymology (2009) vol. 455 pp. 95-125
	"""
	def __init__(self):
		core.FitModel.__init__(self)
		fit_args = self.fit_func_args
		self.params = tuple( [(fit_args[i],i) for i in xrange(len(fit_args))] )
		self.default_params = np.array([7, 0.1, -.53, -4.6, -0.6])
		self.constants = (('n',7),)

	def fit_func(self, x, n, DG_intrinsic, m_intrinsic, DG_interface, m_interface):
		
		# clamp to prevent instability
		if DG_intrinsic<0. or DG_interface>0.:
			return core.FIT_ERROR(x)

		k = np.exp(-(DG_intrinsic - m_intrinsic*x) / core.temperature.RT)
		t = np.exp(-(DG_interface - m_interface*x) / core.temperature.RT)
		pre_factor = (k/(n*(k*t-1))) 
		numerator = n*(k*t)**(n+2) - (n+2)*(k*t)**(n+1) + (n+2)*k*t-n
		denominator = (k*t-1)**2 + k*((k*t)**(n+1) - (n+1)*k*t+n )
		theta = pre_factor * numerator / denominator 
		return theta


"""
==========================================================
KINETIC FOLDING models
==========================================================
"""

class TwoStateChevron(core.FitModel):
	""" Two state chevron plot. 

	k_{obs} = k_u^{H_2O}\exp(m_{ku}x) + k_f^{H_2O}\exp(m_{kf}x)

	Notes:
		Jackson SE and Fersht AR.  Folding of chymotrypsin inhibitor 2. 
		1. Evidence for a two-state transition.
		Biochemistry (1991) 30(43):10428-10435.
	"""
	def __init__(self):
		core.FitModel.__init__(self)
		fit_args = self.fit_func_args
		self.params = tuple( [(fit_args[i],i) for i in xrange(len(fit_args))] )
		self.default_params = np.array([50., 1.3480, 5e-4, 1.])
		#self.constants = (('mf',1.76408),('mu',1.13725))


	def fit_func(self, x, kf, mf, ku, mu):
		k_obs = kf*np.exp(-mf*x) + ku*np.exp(mu*x)
		return k_obs

	def error_func(self, y): 
		return np.log(y)

	@property 
	def equation(self):
		return r'k_{obs} = k_f^{H_2O}\exp(-m_{kf}x) + k_u^{H_2O}\exp(m_{ku}x)'


class ThreeStateChevron(core.FitModel):
	""" Three state chevron with single intermediate.

	k_{obs} = k_{fi}^{H_2O} * exp(-m_{if}*x) +
				k_{if}^{H_2O} * exp((m_i - m_{if})*x) /
				(1 + 1 / K_{iu})

	where:

	K_{iu} = K_{iu}^{H_2O} * exp((m_u-m_i)*x)

	Notes:
		Parker et al. An integrated kinetic analysis of 
		intermediates and transition states in protein folding 
		reactions. 
		Journal of molecular biology (1995) vol. 253 (5) pp. 771-86
	"""
	def __init__(self):
		core.FitModel.__init__(self)
		fit_args = self.fit_func_args
		self.params = tuple( [(fit_args[i],i) for i in xrange(len(fit_args))] )
		self.default_params = np.array([4.5e-4, -9.5e-1, 1.3e9, -6.9,  1.4e-8, -1.6])
		#self.constants = (('mif',-0.97996),('mi',-6.00355),('mu',-1.66154))
		
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
		return r'k_{obs} = k_{fi}^{H_2O}\exp(-m_{if}x) + k_{if}^{H_2O}\exp((m_i - m_{if})x) /(1 + 1 / (K_{iu}^{H_2O}\exp((m_u-m_i)x)))'



class ThreeStateFastPhaseChevron(core.FitModel):
	""" Three state chevron with single intermediate.

	k_{obs} = k_{fi}^{H_2O} * exp(-m_{if}*x) +
				k_{if}^{H_2O} * exp((m_i - m_{if})*x) /
				(1 + 1 / K_{iu})

	where:

	K_{iu} = K_{iu}^{H_2O} * exp((m_u-m_i)*x)

	Notes:
		Parker et al. An integrated kinetic analysis of 
		intermediates and transition states in protein folding 
		reactions. 
		Journal of molecular biology (1995) vol. 253 (5) pp. 771-86
	"""
	def __init__(self):
		core.FitModel.__init__(self)
		fit_args = self.fit_func_args
		self.params = tuple( [(fit_args[i],i) for i in xrange(len(fit_args))] )
		self.default_params = np.array([172., 1.42, .445, .641, 9.41e3, 2.71313, 1.83e-4, 1.06])
		#self.constants = (('kui',172.), ('mui',1.42), ('kiu',.445), ('miu',.641), ('mif',-2.71313),('mfi',1.06534))
		
	def fit_func(self, x, kui, mui, kiu, miu, kif, mif, kfi, mfi):
		k_iu = kiu*np.exp(miu*x)
		k_ui = kui*np.exp(-mui*x)
		k_if = kif*np.exp(-mif*x)
		k_fi = kfi*np.exp(mfi*x)
		K_iu = k_ui / (k_ui+k_iu)
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

	@property 
	def equation(self):
		return r'k_{obs} = k_{fi}^{H_2O}\exp(-m_{if}x) + k_{if}^{H_2O}\exp((m_i - m_{if})x) / (1 + 1 /K_{iu}^{H_2O}\exp((m_u-m_i)x))'



class ThreeStateSequentialChevron(core.FitModel):
	""" Three state metastable intermediate chevron plot.

	A_1 = -(k_{ui}+k_{iu}+k_{if}+k_{fi})
	A_2 = k_{ui}*(k_{if}+k_{fi}) + k_{iu}*k_{uf}

	k_{obs} = 0.5 * (-A_2 - sqrt(A_2^2 - 4*A_1))

	Notes:
		Bachmann and Kiefhaber. Apparent two-state tendamistat 
		folding is a sequential process along a defined route. 
		J Mol Biol (2001) vol. 306 (2) pp. 375-386
	"""
	def __init__(self):
		core.FitModel.__init__(self)
		fit_args = self.fit_func_args
		self.params = tuple( [(fit_args[i],i) for i in xrange(len(fit_args))] )
		self.default_params = np.array([2e4, 0.3480, 20.163, 1.327, 0.3033, 0.2431])
		# self.constants = (('mui',4.34965),('mif',0.68348),('mfi',0.97966))
		
	def fit_func(self, x, kui, mui, kif, mif, kfi, mfi):
		kiu = 1.e4
		miu = 0.
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

	def components(self, x, kui, mui, kif, mif, kfi, mfi):
		kiu = 1.e4
		miu = 0.
		k_ui = kui*np.exp(-mui*x)
		k_iu = kiu*np.exp(miu*x)
		k_if = kif*np.exp(-mif*x)
		k_fi = kfi*np.exp(mfi*x)
		k_TS1 = k_ui + (k_fi/kif)*k_iu
		k_TS2 = (k_ui/k_iu)*k_if + k_fi
		return {'kTS1':k_TS1, 'kTS2':k_TS2}

	@property
	def equation(self):
		return r'k_{obs} = 0.5(-A_2 \pm \sqrt{A_2^2 - 4A_1}) \\ \text{where}\\ A_1 = -(k_{ui} + k_{iu} + k_{if} + k_{fi}) \\A_2 = k_{ui}(k_{if} + k_{fi}) + k_{iu}k_{if} \\ \text{and} \\k_{ui} = k_{ui}^{H_2O}\exp(-m_{ui}x) \\k_{iu} = k_{iu}^{H_2O}\exp(-m_{iu}x) \\ etc...' 



class ParallelTwoStateChevron(core.FitModel):
	""" Two state chevron plot. 

	k_{obs} = k_u^{H_2O}\exp(m_ku*x) + k_u^{H_2O}\exp(m_kf*x)

	Notes:
		[Reference]
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

	@property 
	def equation(self):
		return r'k_{obs} = k_u^{H_2O}\exp(m_{ku}x) + k_u^{H_2O}\exp(m_{kf}x)' 




class ParallelTwoStateUnfoldingChevron(core.FitModel):
	""" Two state chevron plot. 

	k_{obs} = k_u^{H_2O}\exp(m_ku*x) + k_u^{H_2O}\exp(m_kf*x)

	Notes:
		[Reference]
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




class TwoStateChevronMovingTransition(core.FitModel):
	""" Two state chevron with moving transition state.
	Second order polynomial.

	k_u = k_u^{H_2O} * \exp(m_{ku}*x) * \exp(m_{ku}^' * x^2)
	k_f = k_f^{H_2O} * \exp(m_{kf}*x) * \exp(m_{kf}^' * x^2)

	k_{obs} = k_u + k_f

	Notes:
		Ternstrom et al. From snapshot to movie: phi analysis
		of protein folding transition states taken one step 
		further. 
		PNAS (1999) vol. 96 (26) pp. 14854-9
	"""
	def __init__(self):
		core.FitModel.__init__(self)
		fit_args = self.fit_func_args
		self.params = tuple( [(fit_args[i],i) for i in xrange(len(fit_args))] )
		self.default_params = np.array([50., 1.3480, 5e-4, 1., 1.])

	def fit_func(self, x, ku, mu, kf, mf, m_prime):
		k_obs = ku * np.exp(mu*x)*np.exp(m_prime*x*x) + kf*np.exp(-mf*x)*np.exp(-m_prime*x*x)

	def error_func(self, y): 
		return np.log(y)

	@property
	def equation(self):
		return r'k_u = k_u^{H_2O} \cdot \exp(m_{ku}*x) \cdot \exp(m^{\'} x^2) \\ \
				k_f = k_f^{H_2O} \cdot \exp(m_{kf}*x) \cdot \exp(m^{\'} x^2) \\ \
				k_{obs} = k_u + k_f'

	



if __name__ == "__main__":
	get_models()	




