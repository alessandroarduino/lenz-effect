import matplotlib.pyplot as plt
import numpy as np
import os
import scipy.integrate
import seaborn as sns

def ComputeDynamics(Fext, flenz, t_max, dt, q0, q_min, q_max):
    '''Solve the dynamics in a given time interval and spatial domain.

    Parameters
    ----------
    Fext : callable ``Fext(t, q)``
        Forcing term representing the external forces acting on the body
        divided by the body inertia. Both `t` and `q` are scalars. `Fext`
        should return a scalar.
    flenz : callable ``flenz(q)``
        Function representing the unitary Lenz effect divided by the body
        inertia. `q` is a scalar. `flenz` should return a scalar.
    t_max : float
        Endpoint of the integration time interval.
    dt : float
        Time steps at which the output is stored.
    q0 : float
        Initial condition of `q`. The initial condition of its derivative is
        always zero.
    q_min : float
        Minimum value that `q` can assume. The computation is interrupted if
        a lower value is reached.
    q_max : float
        Maximum value that `q` can assume. The computation is interrupted if
        a larger value is reached.

    Returns
    -------
    tt : ndarray
        Time points where the solution is stored. ``tt.shape == (n,)``, being
        `n` the number of stored time points.
    qq : ndarray
        Solution values in each time point. ``qq.shape == (n,2)``, the first
        column contains the actual solution, the second column its derivative.
    '''
    n_max = int(t_max/dt) + 2
    tt = np.zeros((n_max,))
    qq = np.zeros((n_max,2))

    f = lambda t,q: [q[1], Fext(t,q[0]) + flenz(q[0])*q[1]]

    ode = scipy.integrate.ode(f)
    ode.set_integrator('vode', method='adams')
    ode.set_initial_value([q0,0.0], 0)

    n = 0
    tt[n] = ode.t
    qq[n,:] = ode.y
    while ode.successful() and ode.t<t_max and ode.y[0]>q_min and ode.y[0]<q_max:
        ode.integrate(ode.t + dt)
        n += 1
        tt[n] = ode.t
        qq[n,:] = ode.y
    tt = tt[:n+1]
    qq = qq[:n+1,:]

    return tt, qq

def PlotResults(tt_w,qq_w, tt_wo,qq_wo, Fext, flenz, is_angle = False):
    '''Plot the dynamics results.

    Parameters
    ----------
    tt_w, tt_wo : ndarray
        Time points where the solutions with (w) or without (wo) magnetic field
        are stored. ``tt_w.shape == (n_w,)`` and ``tt_wo.shape == (n_wo,)``.
    qq_w, qq_wo : ndarray
        Solution values in each time point with (w) or without (wo) the magnetic
        field. ``qq_w.shape == (n_w,2)`` and ``qq_wo.shape == (n_wo,2)``.
    Fext : callable ``Fext(t, q)``
        Forcing term representing the external forces acting on the body. Both
        `t` and `q` are scalars. `Fext` should return a scalar.
    flenz : callable ``flenz(q)``
        Function representing the unitary Lenz effect. `q` is a scalar. `flenz`
        should return a scalar.
    is_angle : bool, optional
        It is ``True`` if the solution is an angle, in radians. By default it is
        ``False``.
    '''
    identity = lambda x: x
    transform = np.rad2deg if is_angle else identity

    fig,axes = plt.subplots(3,1, figsize=(3.5,7), sharex=True)
    colors = sns.color_palette("Set2")

    axes[0].semilogx(tt_w, transform(qq_w[:,0]), linewidth=2, color=colors[0])
    axes[0].semilogx(tt_wo, transform(qq_wo[:,0]), linewidth=2, color=colors[1])
    axes[0].set_ylabel("Angle (deg)" if is_angle else "Position (m)")
    axes[0].grid()
    axes[0].legend(("With magnet", "Without magnet"))

    axes[1].loglog(tt_w, transform(qq_w[:,1]), linewidth=2, color=colors[0])
    axes[1].loglog(tt_wo, transform(qq_wo[:,1]), linewidth=2, color=colors[1])
    axes[1].set_ylabel("Angular velocity (deg/s)" if is_angle else "Velocity (m/s)")
    axes[1].grid()
    
    FFext = np.array([Fext(t,q) for t,q in zip(tt_w,qq_w[:,0])])
    FFlenz = np.array([flenz(q[0])*q[1] for q in qq_w])
    axes[2].semilogx(tt_w, FFext, linewidth=2, color=colors[2])
    axes[2].semilogx(tt_w, -FFlenz, linewidth=2, color=colors[3])
    axes[2].set_xlabel("Time (s)")
    axes[2].set_ylabel("Moment (Nm)" if is_angle else "Force (N)")
    axes[2].grid()
    axes[2].legend(("Forcing term", "Lenz effect"))

    return fig

def SaveResults(filename, tt,qq, Fext,flenz, is_angle):
    identity = lambda x: x
    transform = np.rad2deg if is_angle else identity

    with open(filename, "w") as ofile:
        ofile.write("time,dof,velocity,ext_force,lenz_force\n")
        for idx in range(len(tt)):
            ofile.write("{:12e},{:.12e},{:.12e},{:.12e},{:.12e}\n".format(
                tt[idx],transform(qq[idx,0]),transform(qq[idx,1]),
                Fext(tt[idx],qq[idx,0]),-flenz(qq[idx,0])*qq[idx,1]))

    return

def ReadLenzEffectFromFile(address, col=1):
    '''Read the Lenz effect with unitary velocity from a text file.

    Parameters
    -----------
    address : string
        Address of the file to read.
    col : int, optional
        Index of the column containing the relevant component of the force
        (or torque). By default it is 1.
    
    Notes
    -----
    The text file is structured in four columns. The first column contains
    the values of the degree of freedom where the forces (or torques) are
    evaluated. The remaining columns contain the components of the relevant
    force or torque.
    '''
    data = np.genfromtxt(address)
    xp = data[:,0]
    fp = data[:,col]
    flenz = lambda x: np.interp(x, xp, fp)
    return flenz
