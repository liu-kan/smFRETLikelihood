from lmfit import minimizer
def mylorentzian(x, amp, cen, wid):
    "lorentzian function: wid = half-width at half-max"
    return (amp / (1 + ((x-cen) / wid)**2))

fitter = minimizer.Minimizer()
fitter.asteval.symtable['lorentzian'] = mylorentzian
print(lorentzian(2,4,6,7))