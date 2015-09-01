import numpy as np 
import scipy as sp

import constants

__author__ = "Alan R. Lowe"
__email__ = "a.lowe@ucl.ac.uk"


def two_state_chevron(x, p0, p=None):
	""" Two state chevron plot. 

	k_{obs} = k_u^{H_2O} + exp(m_ku*x) + k_u^{H_2O} + exp(m_kf*x)

	Notes:
		[Reference]
	"""

	# these options are to allow global fitting
	if not p:
		p = {'kf': p0[0], 'mf': p0[1], 'ku': p0[2], 'mu': p0[3]}

	k_u = p['ku']*np.exp(p['mu']*x)
	k_f = p['kf']*np.exp(p['mf']*x)
	k_obs = k_f + k_u
	return k_obs

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
	

def three_state_chevron(x, p0, p=None):
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
		p = {'kfi':p0[0], 'mif':p0[1], 'kif':p0[2], 'mi':p0[3], \
			'mif':p0[4], 'Kiu':p0[5], 'mu':p0[6]}

	k_fi = p['kfi']*np.exp(-p['mif']*x)
	k_if = p['kif']*np.exp((p['mi'] - p['mif'])*x)
	K_iu = p['Kiu']*np.exp((p['mu'] - p['mi'])*x)
	k_obs = k_fi + k_if / (1.+1./K_iu)
	return k_obs

def three_state_sequential_chevron(x, p0, p=None):
	""" Three state metastable intermediate chevron plot.

	A_1 = -(k_{ui}+k_{iu}+k_{if}+k_{fi})
	A_2 = k_{ui}*(k_{if}+k_{fi}) + k_{iu}*k_{uf}

	k_{obs} = 0.5 * (-A_2 - sqrt(A_2^2 - 4*A_1))

	Notes:
		Bachmann and Kiefhaber. Apparent two-state tendamistat 
		folding is a sequential process along a defined route. 
		J Mol Biol (2001) vol. 306 (2) pp. 375-386
	"""

	if not p:
		p = {'kui': p0[0], 'mui': p0[1], 'kiu': p0[2], 'miu': p0[3], \
			 'kif': p0[4], 'mif': p0[5], 'kfi': p0[6], 'mfi': p0[7]}

	k_ui = p['kui']*np.exp(p['mui']*x)
	k_iu = p['kiu']*np.exp(p['miu']*x)
	k_if = p['kif']*np.exp(p['mif']*x)
	k_fi = p['kfi']*np.exp(p['mfi']*x)
	lam_1 = -(k_ui + k_iu + k_if + kfi)
	lam_2 = k_ui * (k_if+k_fi) + k_iu*kuf
	k_obs = 0.5 * (-lam_2 - np.sqrt(-lam2**2 - 4*lam_1))
	return k_obs

def two_state_equilibrium(x, p0, p=None):
	""" Two state equilbrium denaturation curve.

	Notes:
		Clarke and Fersht. Engineered disulfide bonds as probes of
		the folding pathway of barnase: Increasing the stability 
		of proteins against the rate of denaturation. 
		Biochemistry (1993) vol. 32 (16) pp. 4322-4329
	"""
	
	if not p:
		p = {'alpha_f':p0[0], 'alpha_u':p0[1], 'beta_f':p0[2], \
			 'beta_u':p0[3], 'm':p0[4], 'd50':p0[5]}

	F = (p['alpha_f']+p['beta_f']*x) + (p['alpha_u']+p['beta_u']*x) * \
		( np.exp((p['m']*(x-p['d50'])))/constants.RT) / (1.+np.exp((p['m']*(x-p['d50'])))/constants.RT)

	return F


if __name__ == "__main__":
	pass