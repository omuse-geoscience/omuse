"""
This file contain implimentations of the continuation and timesteppers developed for IEMIC,
but applicable to any code which implements the necessary methods.

"""


from math import isnan, sqrt

def newton(interface, x0, tol=1.e-7, maxit=1000):
    x = x0
    for k in range(maxit):
        fval = interface.rhs(x)
        interface.jacobian(x)
        dx = -interface.solve(fval)

        x = x + dx

        dxnorm = dx.norm()
        if dxnorm < tol:
            print('Newton converged in %d steps with norm %e' % (k, dxnorm))
            break

    return x

def newtoncorrector(interface, par, ds, x, x0, l, l0, tol):
    # Set some parameters
    maxit = 20
    zeta = 1 / x.length()
    delta = 1e-6

    print("ds, l, l0, tol:", ds, l, l0, tol)

    # Do the main iteration
    for k in range(maxit):
        # Set the parameter value and compute F (RHS of 2.2.9)
        interface.set_parameter(par, l)
        fval = interface.rhs(x)

        # Compute F_mu (bottom part of the RHS of 2.2.9)
        interface.set_parameter(par, l + delta)
        dflval = (interface.rhs(x) - fval) / delta
        interface.set_parameter(par, l)

        # Compute the jacobian at x
        interface.jacobian(x)

        # Solve twice with F_x (2.2.9)
        z1 = -interface.solve(fval)
        z2 = interface.solve(dflval)

        # Compute r (2.2.8)
        diff = x - x0
        rnp1 = zeta*diff.dot(diff) + (1-zeta)*(l-l0)**2 - ds**2

        # Compute dl (2.2.13)
        dl = (-rnp1 - 2*zeta*diff.dot(z1)) / (2*(1-zeta)*(l-l0) - 2*zeta*diff.dot(z2))

        # Compute dx (2.2.12)
        dx = z1 - dl*z2

        # Compute a new x and l (2.2.10 - 2.2.11)
        x = x + dx
        l = l + dl

        dxnorm = dx.norm()
        print(k,dxnorm)
        if dxnorm < tol:
            print('Newton corrector converged in %d steps with norm %e' % (k, dxnorm))
            return (x, l)

    raise Exception('No convergence achieved by Newton corrector')

def continuation(interface, x0, par_name, target, ds, maxit, tol=1.e-4):
    x = x0

    # Get the initial tangent (2.2.5 - 2.2.7). 'l' is called mu in Erik's thesis.
    delta = 1e-6
    l = interface.get_parameter(par_name)
    fval = interface.rhs(x)
    interface.set_parameter(par_name, l + delta)
    dl = (interface.rhs(x) - fval) / delta
    interface.set_parameter(par_name, l)

    # Compute the jacobian at x and solve with it (2.2.5)
    interface.jacobian(x)
    dx = -interface.solve(dl)

    # Scaling of the initial tangent (2.2.7)
    dl = 1
    zeta = 1 / x.length()
    nrm = sqrt( zeta*dx.dot(dx) + dl**2)
    print("nrm, zeta:", nrm, zeta, dx.norm())
    dl = dl / nrm
    dx = dx / nrm

    dl0 = dl
    dx0 = dx

    # Perform the continuation
    for j in range(maxit):
        l0 = l
        x0 = x

        # Predictor (2.2.3)
        l = l0 + ds * dl0
        x = x0 + ds * dx0

        # Corrector (2.2.9 and onward)
        x2, l2 = newtoncorrector(interface, par_name, ds, x, x0, l, l0, tol)

        print("%s:" % par_name, l2)

        if (l2 >= target and l0 < target) or (l2 <= target and l0 > target):
            # Converge onto the end point (we usually go past it, so we
            # use Newton to converge)
            l = target;
            interface.set_parameter(par_name, l);
            x = newton(interface, x, tol);

            return x

        # Set the new values computed by the corrector
        dl = l2 - l0
        l = l2
        dx = x2 - x0
        x = x2

        if abs(dl) < 1e-16:
            return x

        # Compute the tangent (2.2.4)
        dx0 = dx / ds
        dl0 = dl / ds

    return x

def newton_time_stepper(interface, x0, theta, dt, tol=1.e-7, maxit=1000):
    x = x0
    b0 = interface.rhs(x)
    sigma = -1 / (theta * dt)
    for k in range(maxit):
        # J - 1 / (theta * dt) * M
        interface.jacobian_with_mass_matrix(x, sigma)
        # M * u_n + dt * theta * F(u_(n+1)) + dt * (1 - theta) * F(u_n) - M * u_(n+1) = 0
        fval = interface.apply_mass_matrix(x0 - x) + dt * theta * interface.rhs(x) + dt * (1 - theta) * b0
        dx = -interface.solve(fval / (theta * dt))

        x = x + dx

        dxnorm = dx.norm()
        if dxnorm < tol:
            print('Newton converged in %d steps with norm %e' % (k, dxnorm))
            break

    return x

def time_stepper(interface, x0, theta, dt, tmax):
    x = x0
    t = 0
    while t < tmax - dt / 2:
        x = newton_time_stepper(interface, x, theta, dt)
        t += dt

        print("t=%f:" % t, interface.rhs(x).norm())

    return x
