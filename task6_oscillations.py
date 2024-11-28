## TASK 6: OSCILLATIONS
def oscillations(a, b, c, d, e, pars):
    # Make a copy of the pars dictionary to avoid modifying the original
    local_pars = pars.copy()
    
    # Update the copy with the new gain values
    local_pars['wEE'] = a   
    local_pars['wEI'] = b   
    local_pars['wIE'] = c   
    local_pars['wII'] = d
    local_pars['I_ext_E'] = e
    

    # Generate nullclines using the updated local parameters
    Exc_null_rE = np.linspace(-0.01, 1, 100)
    Exc_null_rI = nullcline_E(Exc_null_rE, a_E=local_pars['a_E'], theta_E=local_pars['theta_E'], wEE=local_pars['wEE'], wEI=local_pars['wEI'], I_ext_E=local_pars['I_ext_E'])
    Inh_null_rI = np.linspace(-0.01, 1, 100)
    Inh_null_rE = nullcline_I(Inh_null_rI, a_I=local_pars['a_I'], theta_I=local_pars['theta_I'], wIE=local_pars['wIE'], wII=local_pars['wII'], I_ext_I=local_pars['I_ext_I'])


    def system_derivatives(rE, rI, tau_E, a_E, theta_E, wEE, wEI, I_ext_E,
                       tau_I, a_I, theta_I, wIE, wII, I_ext_I, **other_pars):
        drEdt = (-rE + F(wEE * rE - wEI * rI + I_ext_E, a_E, theta_E)) / tau_E
        drIdt = (-rI + F(wIE * rE - wII * rI + I_ext_I, a_I, theta_I)) / tau_I
        return np.array([drEdt, drIdt])

    plt.figure(figsize=(12, 10))
    # Generate a grid
    x, y = np.meshgrid(np.linspace(-2, 2, 40), np.linspace(-2, 2, 40))

    # Compute derivatives for each point on the grid
    u, v = np.zeros(x.shape), np.zeros(y.shape)
    for i in range(x.shape[0]):
        for j in range(x.shape[1]):
            derivatives = system_derivatives(x[i, j], y[i, j], **local_pars)
            u[i, j], v[i, j] = derivatives

    initial_guess = [0.55, 0.2]
    initial_guess2 = [0,0]

    args = (local_pars['a_E'], local_pars['theta_E'], local_pars['wEE'], local_pars['wEI'], local_pars['I_ext_E'],
            local_pars['a_I'], local_pars['theta_I'], local_pars['wIE'], local_pars['wII'], local_pars['I_ext_I'])

    # Find the intersection point using fsolve
    intersection_point = fsolve(nullcline_difference, initial_guess, args=args)
    intersection_point2 = fsolve(nullcline_difference, initial_guess2, args=args)


    # Simulate trajectories using parameters from local_pars
    rE1, rI1 = simulate_wc(local_pars['tau_E'], local_pars['a_E'], local_pars['theta_E'], local_pars['tau_I'], local_pars['a_I'], local_pars['theta_I'],
                           local_pars['wEE'], local_pars['wEI'], local_pars['wIE'], local_pars['wII'], local_pars['I_ext_E'], local_pars['I_ext_I'],
                           rE_init=0.32, rI_init=0.15, dt=local_pars['dt'], range_t=local_pars['range_t'])

    rE2, rI2 = simulate_wc(local_pars['tau_E'], local_pars['a_E'], local_pars['theta_E'], local_pars['tau_I'], local_pars['a_I'], local_pars['theta_I'],
                           local_pars['wEE'], local_pars['wEI'], local_pars['wIE'], local_pars['wII'], local_pars['I_ext_E'], local_pars['I_ext_I'],
                           rE_init=0.32, rI_init=0.15, dt=local_pars['dt'], range_t=local_pars['range_t'])

    plt.figure(figsize=(10, 8))
    plt.xlim(-0.5, 1)
    plt.ylim(-0.5, 0.5)
     # Plot the vector field
    plt.quiver(x, y, u, v, color='black', scale=4, width=0.008)
    plt.plot(Exc_null_rE, Exc_null_rI, 'b', label='E nullcline')
    plt.plot(Inh_null_rE, Inh_null_rI, 'r', label='I nullcline')
    plt.xlabel('E')
    plt.ylabel('I')
    plt.legend(loc='best')
    plt.scatter(intersection_point[0], intersection_point[1], color='green', zorder=5, label='Intersection point')
    plt.scatter(intersection_point2[0], intersection_point2[1], color='green', zorder=5, label='Intersection point')
    print(f"Intersection point 2 (rE, rI): {intersection_point}")
    print(f"Intersection point 1 (rE, rI): {intersection_point2}")
    plt.show()
    my_test_plot(local_pars['range_t'], rE1, rI1, rE2, rI2)

oscillations(6.4, 4.8, 6.0, 1.2, 0.7, pars)
