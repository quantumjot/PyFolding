from numpy import linspace

# version
VERSION = '0.3.0'

# constants
IDEAL_GAS_CONSTANT_KCAL = 1.987204118E-3
TEMPERATURE_CELSIUS = 25.

# calculated constants
ZERO_KELVIN = 273.15
TEMPERATURE_KELVIN = ZERO_KELVIN + TEMPERATURE_CELSIUS
RT = IDEAL_GAS_CONSTANT_KCAL * TEMPERATURE_KELVIN

# error
FITTING_PENALTY = 1e10

# display properties
FONT_SIZE = 18
MARKER_SIZE = 8
LABEL_SIZE = 18
LINE_WIDTH = 3

# error options
CONFIDENCE_INTERVAL = 95

# default range to simulate fit function
XSIM = linspace(0.,10.,100)

if __name__ == "__main__":
	pass
