import inspect
import numpy as np 
import scipy as sp

import constants

__author__ = "Alan R. Lowe"
__email__ = "a.lowe@ucl.ac.uk"


class FitModel(object):
	""" 
	FitModel class

	A generic fit model to enable a common core set of functions
	but specific new member functions to be enabled in derived
	classes.

	Can define parameters in this manner:
	('kf',0), ('mf',1), ... in order to enable paramater sharing
	in global fitting. By default the model just gets the params
	in the order they are defined in the function defininition

	Args:

	Methods:

	Notes:

	"""
	def __init__(self):
		self.__params = None
		self.__param_names = None
		self.__default_params = None

		self.fit_params = None
		self.fit_covar = None
		self.constants = None

	@property 
	def params(self): return self.__params
	@params.setter 
	def params(self, params=None):
		if isinstance(params, tuple):
			self.__params, self.__param_names = [], []
			for key,value in params:
				self.__param_names.append(key)
				self.__params.append(value)
		else:
			raise Warning("Fit parameters must be a tuple")


	@property 
	def name(self): return self.__class__.__name__

	def __call__(self, x, *args):
		""" Parse the fit arguments and pass onto the 
		fitting function
		"""
		fit_args = self.get_fit_params(x, *args)
		return self.error_func( self.fit_func(x, *fit_args) )


	def fit_func(self, x, *args):
		""" The fit function should be defined here """
		raise Exception("Fit function must be defined")

	def error_func(self, y):
		""" The error function should be defined here """
		return y

	def get_fit_params(self, x, *args):
		fit_args = [args[v] for v in self.__params]

		# if we have constants replace the arguments
		# with the constants
		if self.constants: 
			for arg, constant in self.constants:
				if arg in self.__param_names:
					idx = self.__params[ self.__param_names.index(arg) ]
					fit_args[idx] = constant
		return fit_args

	@property 
	def default_params(self):
		""" Give back either the set starting parameters,
		or set all to 1.0 
		"""
		if isinstance(self.__default_params, np.ndarray):
			return self.__default_params
		else:
			return np.ones((len(self.params),1))
	@default_params.setter
	def default_params(self, params):
		if isinstance(params, np.ndarray):
			self.__default_params = params

	@property 
	def fit_func_args(self):
		return inspect.getargspec(self.fit_func).args[2:]

	def test(self):
		print self.fit_func_args
		print self.params






class TwoStateChevron(FitModel):
	""" Two state chevron plot. 

	k_{obs} = k_u^{H_2O} + exp(m_ku*x) + k_u^{H_2O} + exp(m_kf*x)

	Notes:
		[Reference]
	"""
	def __init__(self):
		FitModel.__init__(self)
		fit_args = self.fit_func_args
		self.params = tuple( [(fit_args[i],i) for i in xrange(len(fit_args))] )
		self.default_params = np.array([50., 1.3480, 5e-4, 1.])
		self.constants = (('mf',1.76408),('mu',1.13725))


	def fit_func(self, x, kf, mf, ku, mu):
		k_obs = kf*np.exp(-mf*x) + ku*np.exp(mu*x)
		return k_obs

	def error_func(self, y): 
		return np.log(y)


class ThreeStateChevron(FitModel):
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
		FitModel.__init__(self)
		fit_args = self.fit_func_args
		self.params = tuple( [(fit_args[i],i) for i in xrange(len(fit_args))] )
		self.default_params = np.array([4.5e-4, -9.5e-1, 1.3e9, -6.9,  1.4e-8, -1.6])
		self.constants = (('mif',-0.97996),('mi',-6.00355),('mu',-1.66154))
		
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



class ThreeStateFastPhaseChevron(FitModel):
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
		FitModel.__init__(self)
		fit_args = self.fit_func_args
		self.params = tuple( [(fit_args[i],i) for i in xrange(len(fit_args))] )
		self.default_params = np.array([75.47454, 0.99163, 0.34078, 0.67511, 100., -1.5, 4.5e-4, 1.])
		self.constants = (('kui',75.47454), ('mui',0.99163), ('kiu',0.34078), ('miu',0.67511), \
			('mif',-2.71313),('mfi',1.07534))
		
	def fit_func(self, x, kui, mui, kiu, miu, kif, mif, kfi, mfi):
		k_iu = kiu*np.exp(miu*x)
		k_ui = kui*np.exp(-mui*x)
		k_if = kif*np.exp(mif*x)
		k_fi = kfi*np.exp(mfi*x)
		K_iu = k_iu / (k_ui)
		k_obs = k_fi + k_if / (1.+1./K_iu)
		return k_obs

	def error_func(self, y): 
		return np.log(y)

	def components(self, x, kui, mui, kiu, miu, kif, mif, kfi, mfi):
		k_iu = kiu*np.exp(miu*x)
		k_ui = kui*np.exp(-mui*x)
		k_if = kif*np.exp(mif*x)
		k_fi = kfi*np.exp(mfi*x)
		k_obs_I = k_iu + k_ui
		k_obs_N = k_fi + k_if
		return {'kobs_I':k_obs_I, 'kobs_N':k_obs_N}



class ThreeStateSequentialChevron(FitModel):
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
		FitModel.__init__(self)
		fit_args = self.fit_func_args
		self.params = tuple( [(fit_args[i],i) for i in xrange(len(fit_args))] )
		self.default_params = np.array([2e4, 0.3480, 20.163, 1.327, 0.3033, 0.2431])
		self.constants = (('mui',4.34965),('mif',0.68348),('mfi',0.97966))
		
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



class ParallelTwoStateChevron(FitModel):
	""" Two state chevron plot. 

	k_{obs} = k_u^{H_2O} + exp(m_ku*x) + k_u^{H_2O} + exp(m_kf*x)

	Notes:
		[Reference]
	"""
	def __init__(self):
		FitModel.__init__(self)
		fit_args = self.fit_func_args
		self.params = tuple( [(fit_args[i],i) for i in xrange(len(fit_args))] )
		self.default_params = np.array([50., 1.3480, 5e-4, 1., 150., 3.5])


	def fit_func(self, x, kf_A, mf_A, ku_A, mu_A, kf_B, mf_B):

		if mf_A < 0. or mf_B < 0. or mu_A < 0.:
			return x * 1e10

		if kf_A <0. or ku_A <0. or kf_B < 0.:
			return x * 1e10

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



class TwoStateEquilibrium(FitModel):
	""" Two state equilbrium denaturation curve.

	Notes:
		Clarke and Fersht. Engineered disulfide bonds as probes of
		the folding pathway of barnase: Increasing the stability 
		of proteins against the rate of denaturation. 
		Biochemistry (1993) vol. 32 (16) pp. 4322-4329
	"""
	def __init__(self):
		FitModel.__init__(self)
		fit_args = self.fit_func_args
		self.params = tuple( [(fit_args[i],i) for i in xrange(len(fit_args))] )
		self.default_params = np.array([1., 0.1, 0.0, 0.1, 1.5, 5.])


	def fit_func(self, x, alpha_f, beta_f, alpha_u, beta_u, m, d50):
		F = (alpha_f+beta_f*x) + (alpha_u+beta_u*x) * \
		( np.exp((m*(x-d50)))/constants.RT) / (1.+np.exp((m*(x-d50)))/constants.RT)
		return F





class ParallelTwoStateUnfoldingChevron(FitModel):
	""" Two state chevron plot. 

	k_{obs} = k_u^{H_2O} + exp(m_ku*x) + k_u^{H_2O} + exp(m_kf*x)

	Notes:
		[Reference]
	"""
	def __init__(self):
		FitModel.__init__(self)
		fit_args = self.fit_func_args
		self.params = tuple( [(fit_args[i],i) for i in xrange(len(fit_args))] )
		self.default_params = np.array([5e-4, 1., 1e-5, 1.5])


	def fit_func(self, x, ku_A, mu_A, ku_B, mu_B):

		if mu_A < 0. or mu_B < 0.:
			return x * 1e10

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





def two_state_moving_transition_chevron(x, p0, p=None):
	""" Two state chevron with moving transition state.
	Second order polynomial.

	k_u = k_u^{H_2O} + exp(m_{ku}*x) + exp(m_{ku}^' * x^2)
	k_f = k_f^{H_2O} + exp(m_{kf}*x) + exp(m_{kf}^' * x^2)

	k_{obs} = k_u + k_f

	Notes:
		Ternstrom et al. From snapshot to movie: phi analysis
		of protein folding transition states taken one step 
		further. 
		PNAS (1999) vol. 96 (26) pp. 14854-9
	"""

	if not p:
		p = {'ku': p0[0], 'mu': p0[1], \
			 'kf': p0[2], 'mf': p0[3], 'm_prime': p0[4]}

	k_u = p['ku'] + np.exp(p['mu']*x) + np.exp(p['m_prime']*x**2)
	k_f = p['kf'] + np.exp(p['mf']*x) + np.exp(p['m_prime']*x**2)
	k_obs = k_u + k_f
	return k_obs
	




def three_state_chevron_test(x, p0, p=None):
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
	if not p:
		p = {'kif':p0[0], 'mif':p0[1], 'kfi':p0[2], 'mfi':p0[3],}


	p['kui'] = 100.58949362
	p['mui'] = -1.1097177
	p['kiu'] = 0.37152382
	p['miu'] = 0.66305123
	k_iu = p['kiu']*np.exp(p['miu']*x)
	k_ui = p['kui']*np.exp(p['mui']*x)
	k_if = p['kif']*np.exp(p['mif']*x)
	k_fi = p['kfi']*np.exp(p['mfi']*x)
	K_iu = k_iu / (k_iu + k_ui)
	k_obs = k_fi + k_if / (1.+1./K_iu)
	return k_obs



if __name__ == "__main__":
	pass

