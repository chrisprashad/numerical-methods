# Newton's method for optimization

newton_optimize <- function(g, x0, tol=1e-09, max.iter=100) {
  x = x0  # initial estimate
  g.x = g(x)  # g is the vector with components equal to 0th, 1st and 2nd derivatives
  iter = 0
  while ((abs(g.x[2]) > tol) & (iter < max.iter)) {
    x = x - g.x[2]/g.x[3]
    g.x = g(x)
    iter = iter + 1
  }
  if (iter == max.iter) {
    cat('Algorithm failed to converge. \n')
  } else {
    return(x)
  }
}

# Optimize the Maxwell-Boltzmann distribution for diatomic nitrogen

kB = 1.38064852e-23  # Boltzmann constant
L = 6.022140857e23  # Avagadro constant
R = kB*L  # Gas constant (kJ/(K*mol))
M.N2 = 2*14.007e-03  # molar mass of diatomic nitrogen (kg/mol)
t = 300  # room temperature (K)
a = sqrt(R*t/M.N2)

maxwell_boltzmann = function(x) {
  if (x < 0) return(c(0, 0, 0))
  if (x == 0) return(c(0, 0, NaN))
  prefactor = sqrt(2/pi)*(1/a^3)*exp(-x^2/(2*a^2))
  return( c(prefactor*x^2,
            prefactor*2*x*(1-x^2/(2*a^2)),
            prefactor*(1/a^4)*(2*a^4 - 5*a^2*x^2 + x^4)) )
}

# The speed probability density function of the speed of diatomic nitrogen N2
# at a temperature of 300K (27C)

v = seq(0, 2000, length=100)
plot(v, sqrt(2/pi)*(1/a^3)*v^2*exp(-v^2/(2*a^2)), type="l", lty=1,
     xlab="Speed (m/s)", ylab="Probability density (s/m)")

# Find the maximum using Newton's method starting from x0=500

x.opt.newton = newton_optimize(maxwell_boltzmann, 500)

# Analytical solution

x.opt = sqrt(2)*a

# Relative error of the approximation

abs(1 - x.opt.newton / x.opt)