'''
MuRoSim (Multi-Rotor Simulation) Coupling Script

Authors: Rob Rau and Daniel Weitsman
8/3/21

This module utilizes the wake model in MuRoSim to calculate and return the induced velocities about the rotor disk.
'''

def MuRoSim_wrap(UserIn,geomParams,XsecPolarExp,alphaShaft,omega,U,CT,th,phi,iterations,wake_history_length):
    '''
    This function computes and returns the wake velocities about the rotor disk.
    :param UserIn: UserIn dictionary
    :param geomParams: geomParam dictionary
    :param XsecPolarExp: XsecPolarExp dictionary
    :param omega: rotational rate [rad/s]
    :param U: freestream velocity [m/s]
    :param CT: thrust coefficient
    :param th: computed cyclic pitch controls (th0,th1c,th1s) [rad]
    :param phi: azimuthal vector [rad]
    :return:
    '''

    print('Running MuRoSim')
    import murosim.libmurosim as libmurosim
    import matplotlib.pyplot as plt
    import numpy as np
    from time import monotonic


    def py_step(ac_state, aircraft, ac_input_state, inflows, wake_history, iteration, dt):
        time = iteration * dt

        for rotor_idx in range(aircraft.rotors.length()):
            ac_input_state.rotor_inputs[rotor_idx].cos_aoa = np.cos(
                ac_input_state.rotor_inputs[rotor_idx].angle_of_attack)
            ac_input_state.rotor_inputs[rotor_idx].sin_aoa = np.sin(
                ac_input_state.rotor_inputs[rotor_idx].angle_of_attack)

        for rotor_idx in range(aircraft.rotors.length()):
            libmurosim.compute_rotor_properties(
                aircraft.rotors[rotor_idx],
                ac_state.rotor_states[rotor_idx],
                ac_input_state.rotor_inputs[rotor_idx],
                ac_state,
                inflows[rotor_idx],
                wake_history.history[0],
                ac_state.rotor_states[rotor_idx].C_T,
                time,
                dt
            )

        libmurosim.update_wake(aircraft, ac_state, ac_input_state, wake_history, inflows, iteration, dt)

    dt = (omega * 2 * np.pi / 60) ** -1 / 120
    num_rotors = 1
    num_blades = UserIn['Nb']
    angle_of_attack = alphaShaft
    elements = geomParams['nXsecs']
    d_azimuth = (2.0 * np.pi) / num_rotors
    c = geomParams['chordDist'] / geomParams['R']
    R = geomParams['R']
    r = geomParams['r']
    twist = geomParams['twistDist'] + th[0]
    alpha_0 = XsecPolarExp['Alpha0'][0]
    C_l_alpha = XsecPolarExp['Lift Slope'][0]
    phiRes = len(phi)

    # Create our outer aircraft geometry container
    aircraft = libmurosim.Aircraft(num_rotors)

    # The origins of our 2 rotors
    origins = [libmurosim.Vec3([0, 0, 0])]

    # aircraft.rotors is basically read only iterable so
    # to create all the members we iterate the old fashioned
    # way
    for r_idx in range(aircraft.rotors.length()):
        # Construct the rotor geometry container. We iterate over
        # the blades later to actual set the geometry
        aircraft.rotors[r_idx] = libmurosim.RotorGeometry(
            num_blades,
            origins[r_idx],
            R,
            0  # solidity, we set this later
        )

        # same here.
        for b_idx in range(aircraft.rotors[r_idx].blades.length()):
            # Build the geom of the blades
            aircraft.rotors[r_idx].blades[b_idx] = libmurosim.BladeGeometry(
                elements,
                b_idx * d_azimuth,
                np.sum(c) / len(c)
            )
            # todo set twist on each iteration
            # Convert from linear array format to murosim chunk format
            libmurosim.set_r(aircraft.rotors[r_idx].blades[b_idx], r)
            libmurosim.set_twist(aircraft.rotors[r_idx].blades[b_idx], twist)
            libmurosim.set_alpha_0(aircraft.rotors[r_idx].blades[b_idx], alpha_0)
            libmurosim.set_C_l_alpha(aircraft.rotors[r_idx].blades[b_idx], C_l_alpha)
            libmurosim.set_chord(aircraft.rotors[r_idx].blades[b_idx], c)

        aircraft.rotors[r_idx].solidity = num_blades * aircraft.rotors[r_idx].blades[0].average_chord / (
                np.pi * aircraft.rotors[r_idx].radius)

    # AircraftState is the top level container for holding the current
    # aerodynamic state of the aircraft. It breaks down into rotors
    # the the individual blades. There is a series of functions provided
    # to turn internal state data into a linear array.
    ac_state = libmurosim.AircraftState(num_rotors, num_blades, elements)

    # Create and setup the input state. This would be the
    # sort of input a dynamics simulator might feed into
    # the aero model.
    ac_input_state = libmurosim.AircraftInputState(num_rotors)

    ac_input_state.rotor_inputs[0].angle_of_attack = -angle_of_attack
    ac_input_state.rotor_inputs[0].angular_velocity = omega
    ac_input_state.rotor_inputs[0].angular_accel = 0
    ac_input_state.rotor_inputs[0].freestream_velocity = U
    ac_state.rotor_states[0].C_T = CT

    # Create our inflow models. One for each rotor
    inflows = [libmurosim.BeddosInflow()]

    # Setup the wake history. We need at minimum 2 timesteps worth of history for the update.
    # Increasing the history increases computation time with the current implementation
    wake_history = libmurosim.WakeHistory(num_rotors, num_blades, wake_history_length, 2)

    start_time = monotonic()

    for iteration in range(iterations):
        if iteration % 100 == 0:
            now = monotonic()
            print(now - start_time, ": iteration: ", iteration)
            start_time = now

        py_step(ac_state, aircraft, ac_input_state, inflows, wake_history, iteration, dt)

    print("rotor 0 C_T: ", ac_state.rotor_states[0].C_T)

    chunk_size = 8
    coord = np.zeros((phiRes, elements, 3))
    coord[:, :, 0] = np.expand_dims(np.cos(phi), axis=1) * r
    coord[:, :, 1] = np.expand_dims(np.sin(phi), axis=1) * r
    coord = coord.reshape(int(phiRes * elements / chunk_size), chunk_size, 3)

    # This how the induced velocities produced by the wake are computed. Both the
    # input coords and the output velocities are in the global coordinate space
    # both the xyz inputs and velocity outputs are non-dimensionalized
    induced_velocities = list(map(
        lambda x, y, z: libmurosim.compute_wake_induced_velocities(wake_history.history[0], list(x), list(y), list(z),
                                                                   ac_state, omega), coord[:, :, 0], coord[:, :, 1],coord[:, :, 2]))
    induced_velocities = np.transpose(np.array([[x.v_x, x.v_y, x.v_z] for x in induced_velocities]),
                                      axes=(0, 2, 1)).reshape((phiRes, elements, 3))

    if UserIn['plotWake']:
        max_dim_1 = -np.inf
        min_dim_1 = np.inf

        max_dim_2 = -np.inf
        min_dim_2 = np.inf

        # Iterating directly on the collections like this is read only
        for r_idx, rotor in enumerate(ac_state.rotor_states):
            for b_idx, blade in enumerate(rotor.blade_states):
                x = libmurosim.get_wake_x_component(
                    wake_history.history[0].rotor_wakes[r_idx].vortex_filaments[b_idx])
                z = libmurosim.get_wake_z_component(
                    wake_history.history[0].rotor_wakes[r_idx].vortex_filaments[b_idx])

                max_dim_2 = max(max_dim_2, max(z))
                min_dim_2 = min(min_dim_2, min(z))

                max_dim_1 = max(max_dim_1, max(x))
                min_dim_1 = min(min_dim_1, min(x))
                plt.plot(x, z, linewidth=0.5)

                x_r = np.cos(blade.azimuth)
                x_b = aircraft.rotors[r_idx].origin[0] + x_r * np.cos(
                    ac_input_state.rotor_inputs[r_idx].angle_of_attack)
                z_b = aircraft.rotors[r_idx].origin[2] + x_r * np.sin(
                    ac_input_state.rotor_inputs[r_idx].angle_of_attack)

                plt.plot([aircraft.rotors[r_idx].origin[0], x_b], [aircraft.rotors[r_idx].origin[2], z_b], "k-",
                         linewidth=1.5)
        plt.axis("square")

        span_1 = max_dim_1 - min_dim_1
        span_2 = max_dim_2 - min_dim_2
        plt.xlim(left=min_dim_1 - .1 * span_1, right=max_dim_1 + .1 * span_1)
        plt.ylim(bottom=min_dim_2 - .9 * span_2, top=max_dim_2 + .9 * span_2)
        plt.show()

    return induced_velocities