# Class Parameter Docstrings (in alphabetical order)

## Airfoil

### main (pterasoftware\geometry.py)
```python
        """
        :param name: str, optional
            This is the name of the airfoil. It should correspond to the name in the
            airfoils directory unless you are passing in your own coordinates. The
            default is "Untitled Airfoil".
        :param coordinates: array, optional
            This is an N x 2 array of the airfoil's coordinates, where N is the
            number of coordinates. Treat this as an immutable, don't edit directly
            after initialization. If you wish to load coordinates from the airfoil
            directory, leave this as None. The default is None. Make sure that any
            airfoil coordinates used range in x from 0 to 1.
        :param repanel: bool, optional
            This is the variable that determines whether you would like to repanel
            the airfoil coordinates. This applies to coordinates passed in by the
            user or to the directory coordinates. I highly recommended setting this
            to True. The default is True.
        :param n_points_per_side: int, optional
            This is number of points to use when repaneling the airfoil. It is
            ignored if the repanel is False. The default is 400.
        """
```

### feature/improved_geometry_definitions (pterasoftware\geometry\airfoil.py)
```python
        """
        :param name: str, optional

            This is the name of the Airfoil. It should correspond to the name of a
            file the airfoils directory (once converted to lower-case and stripped of
            leading and trailing whitespace) unless you are passing in your own
            array of points using outline_A_lp. The default is "NACA0012".

        :param outline_A_lp: array-like of shape (N,2), optional

            This is an array of the 2D points making up the Airfoil's outline (in
            airfoil axes, relative to the leading point). If you wish to load
            coordinates from the airfoils directory, leave this as None, which is the
            default. If not, it must be an array-like object of numbers (int or
            float) with shape (N,2). It can be a tuple, list, or numpy array. Values
            are converted to floats internally.it must be a Nx2 numpy array of
            numbers (int or float). Make sure all x-component values are in the range
            [ 0.0, 1.0]. The default value is None.

        :param resample: boolLike, optional

            This is the variable that determines whether you would like to resample
            the points defining the Airfoil's outline. This applies to points passed
            in by the user or to those from the airfoils directory. I highly
            recommended setting this to True. It can be a boolean or a NumPy boolean
            and will be converted internally to a boolean. The default is True.

        :param n_points_per_side: int or None, optional

            This is number of points to use when creating the Airfoil's MCL and when
            resampling the upper and lower parts of the Airfoil's outline. It must be
            a positive int greater than or equal to 3. The resampled outline will
            have a total number of points equal to (2 * n_points_per_side) - 1. I
            highly recommend setting this to at least 100. The default value is 400.
        """
```

## Airplane

### main (pterasoftware\geometry.py)
```python
        """
        :param wings: list of Wing objects
            This is a list of the airplane's wings defined as Wing objects. It must
            contain at least one Wing object.
        :param name: str, optional
            A sensible name for your airplane. The default is "Untitled Airplane".
        :param x_ref: float, optional
            This is the x coordinate of the moment reference point. It should be the
            x coordinate of the center of gravity. The default is 0.0.
        :param y_ref: float, optional
            This is the y coordinate of the moment reference point. It should be the
            y coordinate of the center of gravity. The default is 0.0.
        :param z_ref: float, optional
            This is the z coordinate of the moment reference point. It should be the
            z coordinate of the center of gravity. The default is 0.0.
        :param weight: float, optional
            This parameter holds the weight of the aircraft in Newtons. This is used
            by the trim functions. The default value is 0.0.
        :param s_ref: float, optional if more than one wing is in the wings list.
            This is the reference wetted area. If not set, it populates from first
            wing object.
        :param c_ref: float, optional if more than one wing is in the wings list.
            This is the reference chord length. If not set, it populates from first
            wing object.
        :param b_ref: float, optional if more than one wing is in the wings list.
            This is the reference calculate_span. If not set, it populates from first
            wing object.
        """
```

### feature/improved_geometry_definitions (pterasoftware\geometry\airplane.py)
```python
        """
        :param wings: list of Wings

            This is a list of the airplane's wings defined as Wings. It must contain
            at least one Wing. Wings with symmetric=True and non-coincident symmetry
            planes will be automatically processed into separate Wings during
            initialization (type 5 symmetry).

        :param name: str, optional

            A sensible name for your airplane. The default is "Untitled Airplane".

        :param Cgi_E_I: array-like of 3 numbers, optional

            Position [x, y, z] of this Airplane's starting point (in Earth axes,
            relative to the simulation's starting point). Can be a list, tuple,
            or numpy array of numbers (int or float). Values are converted to floats
            internally. For the first Airplane in a simulation, this must be [0.0,
            0.0, 0.0] since the simulation's starting point is defined as the first
            Airplane's starting point (the location of its CG at t=0). The default is
            (0.0, 0.0, 0.0).

        :param angles_E_to_B_izyx: array-like of 3 numbers, optional

            Angles [angle1, angle2, angle3] from Earth axes to body axes using an
            intrinsic 3-2'-1" sequence. Can be a tuple, list, or numpy array of
            numbers (int or float). Values are converted to floats internally. This
            defines the orientation of the airplane's body axes relative to Earth
            axes. Note that body axes differ from geometry axes: body axes point
            forward/right/down while geometry axes point aft/right/up. The units are
            degrees. All angles must lie in the range (-180.0, 180.0] degrees. The
            default is (0.0, 0.0, 0.0).

        :param weight: number, optional

            This parameter is a number (int or float) that represents the weight of
            the aircraft in Newtons. This is used by the trim functions. It must be
            greater than or equal to zero. The default value is 0.0.

        :param s_ref: number, optional

           This parameter is a number (int or float) that represents the reference
           wetted area. If not set or set to None (the default value), it populates
           from first Wing. If set, it must be greater than zero. The units are
           square meters.

        :param c_ref: float, optional

            This parameter is a number (int or float) that represents the reference
            chord length. If not set or set to None (the default value), it populates
            from first Wing. If set, it must be greater than zero. The units are meters.

        :param b_ref: float, optional

            This parameter is a number (int or float) that represents the reference
            span. If not set or set to None (the default value), it populates from
            first Wing. If set, it must be greater than zero. The units are meters.
        """
```

## AirplaneMovement

### main (pterasoftware\movement.py)
```python
        """
        :param base_airplane: Airplane
            This is the first airplane object, from which the others will be created.
        :param wing_movements: list of WingMovement objects
            This is a list of the WingMovement objects associated with each of the
            base airplane's wings.
        :param x_ref_amplitude: float, optional
            This is the amplitude of the airplane's change in its x reference point.
            Its units are meters and its default value is 0 meters.
        :param x_ref_period: float, optional
            This is the period of the airplane's change in its x reference point. Its
            units are seconds and its default value is 0 seconds.
        :param x_ref_spacing: string, optional
            This value determines the spacing of the airplane's change in its x
            reference point. The options are "sine", and "uniform". The default value
            is "sine".
        :param y_ref_amplitude: float, optional
            This is the amplitude of the airplane's change in its y reference point.
            Its units are meters and its default value is 0 meters.
        :param y_ref_period: float, optional
            This is the period of the airplane's change in its y reference point. Its
            units are seconds and its default value is 0 seconds.
        :param y_ref_spacing: string, optional
            This value determines the spacing of the airplane's change in its y
            reference point. The options are "sine", and "uniform". The default value
            is "sine".
        :param z_ref_amplitude: float, optional
            This is the amplitude of the airplane's change in its z reference point.
            Its units are meters and its default value is 0 meters.
        :param z_ref_period: float, optional
            This is the period of the airplane's change in its z reference point. Its
            units are seconds and its default value is 0 seconds.
        :param z_ref_spacing: string, optional
            This value determines the spacing of the airplane's change in its z
            reference point. The options are "sine", and "uniform". The default value
            is "sine".
        """
```

### feature/improved_geometry_definitions (pterasoftware\movements\airplane_movement.py)
```python
        """
        :param base_airplane: Airplane

            This is the base Airplane, from which the Airplane at each time step will
            be created.

        :param wing_movements: list of WingMovements

            This is a list of the WingMovement associated with each of the base
            Airplane's Wings. It must have the same length as the base Airplane's
            list of Wings.

        :param ampCgi_E_I: array-like of 3 numbers, optional

            The amplitudes of the AirplaneMovement's changes in its Airplanes'
            Cgi_E_I parameters. Can be a tuple, list, or numpy array of non-negative
            numbers (int or float). Values are converted to floats internally. The
            default value is ( 0.0, 0.0, 0.0). The units are in meters.

        :param periodCgi_E_I: array-like of 3 numbers, optional

            The periods of the AirplaneMovement's changes in its Airplanes' Cgi_E_I
            parameters. Can be a tuple, list, or numpy array of non-negative numbers
            (int or float). Values are converted to floats internally. The default
            value is (0.0, 0.0, 0.0). Each element must be 0.0 if the corresponding
            element in ampCgi_E_I is 0.0 and non-zero if not. The units are in seconds.

        :param spacingCgi_E_I: array-like of 3 strs or callables, optional

            The value determines the spacing of the AirplaneMovement's change in its
            Airplanes' Cgi_E_I parameters. Can be a tuple, list, or numpy array. Each
            element can be the string "sine", the string "uniform", or a callable
            custom spacing function. Custom spacing functions are for advanced users
            and must start at 0, return to 0 after one period of 2*pi radians, have
            zero mean, have amplitude of 1, be periodic, return finite values only,
            and accept a ndarray as input and return a ndarray of the same shape. The
            custom function is scaled by ampCgi_E_I, shifted by phaseCgi_E_I, and
            centered around the base value, with the period controlled by
            periodCgi_E_I. The default value is ("sine", "sine", "sine").

        :param phaseCgi_E_I: array-like of 3 numbers, optional

            The phase offsets of the elements in the first time step's Airplane's
            Cgi_E_I parameter relative to the base Airplane's Cgi_E_I parameter. Can
            be a tuple, list, or numpy array of non-negative numbers (int or float)
            in the range [0.0, 360.0). Values are converted to floats internally. The
            default value is (0.0, 0.0, 0.0). Each element must be 0.0 if the
            corresponding element in ampCgi_E_I is 0.0 and non-zero if not. The units
            are in degrees.

        :param ampAngles_E_to_B_izyx: array-like of 3 numbers, optional

            The amplitudes of the AirplaneMovement's changes in its Airplanes'
            angles_E_to_B_izyx parameters. Can be a tuple, list, or numpy array of
            numbers (int or float) in the range [0.0, 360.0). Values are converted to
            floats internally. The default value is (0.0, 0.0, 0.0). The units are in
            degrees.

        :param periodAngles_E_to_B_izyx: array-like of 3 numbers, optional

            The periods of the AirplaneMovement's changes in its Airplanes'
            angles_E_to_B_izyx parameters. Can be a tuple, list, or numpy array of
            non-negative numbers (int or float). Values are converted to floats
            internally. The default value is (0.0, 0.0, 0.0). Each element must be
            0.0 if the corresponding element in ampAngles_E_to_B_izyx is 0.0 and
            non-zero if not. The units are in seconds.

        :param spacingAngles_E_to_B_izyx: array-like of 3 strs or callables, optional

            The value determines the spacing of the AirplaneMovement's change in its
            Airplanes' angles_E_to_B_izyx parameters. Can be a tuple, list, or numpy
            array. Each element can be the string "sine", the string "uniform", or a
            callable custom spacing function. Custom spacing functions are for advanced
            users and must start at 0, return to 0 after one period of 2*pi radians,
            have zero mean, have amplitude of 1, be periodic, return finite values
            only, and accept a ndarray as input and return a ndarray of the same shape.
            The custom function is scaled by ampAngles_E_to_B_izyx, shifted by
            phaseAngles_E_to_B_izyx, and centered around the base value, with the
            period controlled by periodAngles_E_to_B_izyx. The default value is
            ("sine", "sine", "sine").

        :param phaseAngles_E_to_B_izyx: array-like of 3 numbers, optional

            The phase offsets of the elements in the first time step's Airplane's
            angles_E_to_B_izyx parameter relative to the base Airplane's
            angles_E_to_B_izyx parameter. Can be a tuple, list, or numpy array of
            numbers (int or float) in the range [0.0, 360.0). Values are converted to
            floats internally. The default value is (0.0, 0.0, 0.0). Each element
            must be 0.0 if the corresponding element in ampAngles_E_to_B_izyx is 0.0
            and non-zero if not. The units are in degrees.
        """
```

## HorseshoeVortex

### main (pterasoftware\aerodynamics.py)
```python
        """
        :param front_right_vortex_vertex: 1D array
            This is a vector containing the x, y, and z coordinates of the vortex's
            front right point. It's a (,3) array with units of meters.
        :param front_left_vortex_vertex: 1D array
            This is a vector containing the x, y, and z coordinates of the vortex's
            front left point. It's a (,3) array with units of meters.
        :param strength: float
            This is the strength of the vortex in meters squared per second.
        :param left_vortex_leg_direction: 1D array
            This is a vector that points in the direction of the left leg of the
            horseshoe vortex. It should be a (,3) array.
        :param infinite_leg_length: float
            This is the length of the infinite legs in meters.
        """
```

### feature/improved_geometry_definitions (pterasoftware\aerodynamics.py)
```python
    """
    :param Frhvp_G_Cg: array-like of 3 numbers

        Position [x, y, z] of the HorseshoeVortex's front-right point (in geometry
        axes, relative to the CG). The front-right point is defined as the start
        point of the HorseshoeVortex's front leg, which is also its one finite leg.
        Can be a list, tuple, or numpy array of numbers (int or float). Values are
        converted to floats internally. The units are in meters.

    :param Flhvp_G_Cg: array-like of 3 numbers

        Position [x, y, z] of the HorseshoeVortex's front-left point (in geometry
        axes, relative to the CG). The front-left point is defined as the end point
        of the HorseshoeVortex's front leg, which is also its one finite leg. Can be
        a list, tuple, or numpy array of numbers (int or float). Values are converted
        to floats internally. The units are in meters.

    :param leftLegVector_G: array-like of 3 numbers

        Direction vector of the HorseshoeVortex's left leg (in geometry axes). The
        left leg starts from the front-left point and ends at the back-left point. It
        is one of the HorseshoeVortex's two quasi-infinite legs, the other being the
        right leg. It can be a list, tuple, or numpy array of numbers (int or float).
        Values are converted to floats internally. If this isn't already a unit
        vector, it will be converted to one during initialization. The right leg's
        vector (in geometry axes) is defined as -1 times this vector. The units are
        in meters.

    :param left_right_leg_lengths: number

        This is the length of the HorseshoeVortex's left and right quasi-infinite
        legs. It must be a positive number and will be converted internally to a
        float. I recommend setting it to at least 20 times the length of the finite
        leg. The units are in meters.

    :param strength: number

        This is the strength of the HorseshoeVortex. It must be a number and will be
        converted internally to a float. Its units are in meters squared per second.
    """
```

## LineVortex/_LineVortex

### main (pterasoftware\aerodynamics.py)
```python
        """
        :param origin: 1D array
            This is a vector containing the x, y, and z coordinates of the line
            vortex's origin. It's a (,3) array with units of meters.
        :param termination: 1D array
            This is a vector containing the x, y, and z coordinates of the line
            vortex's terminus. It's a (,3) array with units of meters.
        :param strength: float
            This is the strength of the vortex in meters squared per second.
        """
```

### feature/improved_geometry_definitions (pterasoftware\aerodynamics.py)
```python
        """
        :param Slvp_G_Cg: array-like of 3 numbers

            Position [x, y, z] of the LineVortex's start point (in geometry axes,
            relative to the CG). Can be a list, tuple, or numpy array of numbers (int
            or float). Values are converted to floats internally. The units are in
            meters.

        :param Elvp_G_Cg: array-like of 3 numbers

            Position [x, y, z] of the LineVortex's end point (in geometry axes,
            relative to the CG). Can be a list, tuple, or numpy array of numbers (int
            or float). Values are converted to floats internally. The units are in
            meters.

        :param strength: number

            This is the strength of the LineVortex. It must be a number and will be
            converted internally to a float. Its units are in meters squared per second.
        """
```

## Movement

### main (pterasoftware\movement.py)
```python
        """
        :param airplane_movements: list of AirplaneMovement objects
            This is a list of objects which characterize the movement of
            each airplane in the problem.
        :param operating_point_movement: OperatingPointMovement
            This object characterizes the movement of the operating point.
        :param num_steps: int, optional
            This integer is the number of time steps of the unsteady simulation. If
            not given a value, and the movement is dynamic, this method will
            calculate one such that the simulation will cover some number of cycles
            of the maximum period movement. The number of cycles defaults to three,
            but can be changed with the num_cycles parameter. If not given a value,
            and the movement is static, the number of steps will default to the
            number of time steps such that the wake extends back by some number of
            reference chord lengths. The number of chord lengths defaults to ten,
            but can be changed with the num_chords parameter.
        :param num_cycles: int, optional
            This integer is the number of cycles of the maximum period movement used
            to calculate a non-populated num_steps parameter. This parameter is only
            used if the num_steps parameter is None, and the movement isn't static.
            The default value is None.
        :param num_chords: int, optional
            This integer is the number of reference chord lengths used to calculate a
            non-populated num_steps parameter. This parameter is only used if the
            num_steps parameter is None, and the movement is static. The default
            value is None.
        :param delta_time: float, optional
            This float is the time, in seconds, between each time current_step. If
            not given a value, this method will calculate one such the ring vortices
            shed off the main wing will have roughly the same chord length as the
            panels on the main wing. This is based on the base airplane's reference
            chord length, its main wing's number of chordwise panels, and its base
            operating point's velocity.
        """
```

### feature/improved_geometry_definitions (pterasoftware\movements\movement.py)
```python
        """
        :param airplane_movements: list of AirplaneMovements

            This is a list of objects which characterize the movement of each
            of the airplanes in the UnsteadyProblem.

        :param operating_point_movement: OperatingPointMovement

            This object characterizes changes to the UnsteadyProblem's the operating
            point.

        :param delta_time: number or None, optional

            delta_time is the time, in seconds, between each time step. If left as
            None, which is the default value, Movement will calculate a value such
            that RingVortices shed from the first Wing will have roughly the same
            chord length as the RingVortices on the first Wing. This is based on
            first base Airplane's reference chord length, its first Wing's number of
            chordwise panels, and its base OperatingPoint's velocity. If set,
            delta_time must be a positive number (int or float). It will be converted
            internally to a float.

        :param num_cycles: int or None, optional

            num_cycles is the number of cycles of the maximum period motion used to
            calculate a non-populated num_steps parameter if Movement isn't static.
            If num_steps is set or Movement is static, this must be left as None,
            which is the default value. If num_steps isn't set and Movement isn't
            static, num_cycles must be a positive int. In that case, I recommend
            setting num_cycles to 3.

        :param num_chords: int or None, optional

            num_chords is the number of chord lengths used to calculate a
            non-populated num_steps parameter if Movement is static. If num_steps is
            set or Movement isn't static, this must be left as None, which is the
            default value. If num_steps isn't set and Movement is static, num_chords
            must be a positive int. In that case, I recommend setting num_chords to
            10. For cases with multiple Airplanes, the num_chords will reference the
            largest reference chord length.

        :param num_steps: int or None, optional

            num_steps is the number of time steps of the unsteady simulation. It must
            be a positive int. The default value is None. If left as None,
            and Movement isn't static, Movement will calculate a value such that the
            simulation will cover some number of cycles of the maximum period of all
            the motion described in Movement's sub-movement objects, sub-sub-movement
            objects, etc. If num_steps is left as None, and Movement is static,
            it will default to the number of time steps such that the wake extends
            back by some number of reference chord lengths.
        """
```

## OperatingPoint

### main (pterasoftware\operating_point.py)
```python
        """
        :param density: float, optional
            This parameter is the density. The units are kilograms per meters cubed.
            The default value is 1.225.
        :param velocity: float, optional
            This parameter is the freestream speed in the positive x direction. The
            units are meters per second. The
            default value is 10.0.
        :param alpha: float, optional
            This parameter is the angle of attack. The units are degrees. The default
            value is 5.0.
        :param beta: float, optional
            This parameter is the sideslip angle. The units are degrees. The default
            value is 0.0.
        :param nu: float, optional
            This parameter is the air's kinematic viscosity. The units are meters
            squared per second. This parameter is only used in the unsteady ring
            vortex lattice method's vortex core growth model. The default value is
            15.06e-6 meters squared per second, which corresponds to air's kinematic
            viscosity at 20 degrees Celsius [source:
            https://www.engineeringtoolbox.com].
        """
```

### feature/improved_geometry_definitions (pterasoftware\operating_point.py)
```python
        """
        :param rho: number, optional

            This parameter is the fluid's density. It must be a positive number and
            will be converted internally to a float. The units are kilograms per
            meters cubed. The default value is 1.225.

        :param vCg__E: number, optional

            This parameter is the speed of the Airplane's CG (observed from the Earth
            frame). Given that (1) this is the magnitude of a vector, and (2) we
            always assume a still fluid in our simulations, this value is equivalent
            to the freestream speed (the speed of the apparent wind, infinitely far
            away from the Airplane, observed while moving at the same speed as the
            Airplane's non-accelerating CG). It must be a positive number and will be
            converted internally to a float. Its units are in meters per second. The
            default value is 10.0.

        :param alpha: number, optional

            This parameter is the angle of attack. For more details on the exact
            interpretation of this value, see the description of wind axes in
            docs/AXES_POINTS_AND_FRAMES.md. It must be a number in the range (-180.0,
            180.0] and will be converted internally to a float. The units are
            degrees. The default value is 5.0.

        :param beta: number, optional

            This parameter is the sideslip angle. For more details on the exact
            interpretation of this value, see the description of wind axes in
            docs/AXES_POINTS_AND_FRAMES.md. It must be a number in the range (-180.0,
            180.0] and will be converted internally to a float. The units are
            degrees. The default value is 0.0.

        :param externalFX_W: number, optional

            This parameter is for any additional thrust or drag on the Airplane's
            body (in wind axes) not due to the Airplane's Wings. It is useful for
            trim analyses. It must be a number and will be converted internally to a
            float. The units are Newtons. The default value is 0.0.

        :param nu: number, optional

            This parameter is the fluid's kinematic viscosity. The units are meters
            squared per second. This parameter is only used in the unsteady ring
            vortex lattice method's vortex core growth model. It must be a positive
            number and will be converted internally to a float. Its units are in
            meters squared per second. The default value is 15.06e-6,
            which corresponds to air's kinematic viscosity at 20 degrees Celsius
            [source: https://www.engineeringtoolbox.com].
        """
```

## OperatingPointMovement

### main (pterasoftware\movement.py)
```python
        """
        :param base_operating_point: OperatingPoint
            This is the operating point object, from which the others will be created.
        :param velocity_amplitude: float, optional
            This is the amplitude of the operating point's change in velocity. Its
            units are meters per second and its default value is 0 meters per second.
        :param velocity_period: float, optional
            This is the period of the operating point's change in its velocity. Its
            units are seconds and its default value is 0 seconds.
        :param velocity_spacing: string, optional
            This value determines the spacing of the operating point's change in its
            velocity. The options are "sine", and "uniform". The default value is
            "sine".
        """
```

### feature/improved_geometry_definitions (pterasoftware\movements\operating_point_movement.py)
```python
        """
        :param base_operating_point: OperatingPoint

            This is the base OperatingPoint, from which the OperatingPoint at each
            time step will be created.

        :param ampVCg__E: number, optional

            The amplitude of the OperatingPointMovement's changes in its
            OperatingPoints' vCg__E parameters. Must be a non-negative number (int or
            float), and is converted to a float internally. The default value is 0.0.
            The units are in meters per second.

        :param periodVCg__E: number, optional

            The period of the OperatingPointMovement's changes in its
            OperatingPoints' vCg__E parameter. Must be a non-negative number (int or
            float), and is converted to a float internally. The default value is 0.0.
            It must be 0.0 if ampVCg__E 0.0 and non-zero if not. The units are in
            seconds.

        :param spacingVCg__E: string, optional

            The value determines the spacing of the OperatingPointMovement's change
            in its OperatingPoints' Cgi_E_I parameters. Must be either "sine" or
            "uniform". The default value is "sine".

        :param phaseVCg__E: number optional

            The phase offsets of the first time step's OperatingPoint's vCg__E
            parameter relative to the base OperatingPoint's vCg__E parameter. Must be
            a number (int or float) in the range [0.0, 360.0), and is converted to a
            float internally. The default value is 0.0. It must be 0.0 if ampVCg__E
            is 0.0 and non-zero if not. The units are in degrees.
        """
```

## Panel

### main (pterasoftware\panel.py)
```python
        """
        :param front_right_vertex: (3,) array
            This is an array containing the x, y, and z coordinates of the panel's
            front right vertex.
        :param front_left_vertex: (3,) array
            This is an array containing the x, y, and z coordinates of the panel's
            front left vertex.
        :param back_left_vertex: (3,) array
            This is an array containing the x, y, and z coordinates of the panel's
            back left vertex.
        :param back_right_vertex: (3,) array
            This is an array containing the x, y, and z coordinates of the panel's
            back right vertex.
        :param is_leading_edge: bool
            This is true if the panel is a leading edge panel on a wing, and false
            otherwise.
        :param is_trailing_edge: bool
            This is true if the panel is a trailing edge panel on a wing, and false
            otherwise.
        """
```

### feature/improved_geometry_definitions (pterasoftware\geometry\panel.py)
```python
        """
        :param Frpp_G_Cg: array-like of 3 numbers

            Position [x, y, z] of the Panel's front right vertex (in geometry axes,
            relative to the CG). Can be a list, tuple, or numpy array of numbers (int
            or float). Values are converted to floats internally. Front, back, left,
            and right are defined with respect to this Panel's position on its wing
            section, with front meaning towards the leading edge, right meaning
            towards the next wing section or the Wing's tip, etc. The units are in
            meters.

        :param Flpp_G_Cg: array-like of 3 numbers

            Position [x, y, z] of the Panel's front left vertex (in geometry axes,
            relative to the CG). Can be a list, tuple, or numpy array of numbers (int
            or float). Values are converted to floats internally. Front, back, left,
            and right are defined with respect to this Panel's position on its wing
            section, with front meaning towards the leading edge, right meaning
            towards the next wing section or the Wing's tip, etc. The units are in
            meters.

        :param Blpp_G_Cg: array-like of 3 numbers

            Position [x, y, z] of the Panel's back left vertex (in geometry axes,
            relative to the CG). Can be a list, tuple, or numpy array of numbers (int
            or float). Values are converted to floats internally. Front, back, left,
            and right are defined with respect to this Panel's position on its wing
            section, with front meaning towards the leading edge, right meaning
            towards the next wing section or the Wing's tip, etc. The units are in
            meters.

        :param Brpp_G_Cg: array-like of 3 numbers

            Position [x, y, z] of the Panel's back right vertex (in geometry axes,
            relative to the CG). Can be a list, tuple, or numpy array of numbers (int
            or float). Values are converted to floats internally. Front, back, left,
            and right are defined with respect to this Panel's position on its wing
            section, with front meaning towards the leading edge, right meaning
            towards the next wing section or the Wing's tip, etc. The units are in
            meters.

        :param is_leading_edge: boolLike

            This is True if the Panel is a leading edge Panel on a Wing, and False
            otherwise. It can be a boolean or a NumPy boolean and will be converted
            internally to a boolean.

        :param is_trailing_edge: boolLike

            This is True if the Panel is a trailing edge Panel on a Wing, and False
            otherwise. It can be a boolean or a NumPy boolean and will be converted
            internally to a boolean.
        """
```

## RingVortex

### main (pterasoftware\aerodynamics.py)
```python
        """
        :param front_left_vertex: 1D array
            This is a vector containing the x, y, and z coordinates of the vortex's
            front left point. It's a (,3) array with units of meters.
        :param front_right_vertex: 1D array
            This is a vector containing the x, y, and z coordinates of the vortex's
            front right point. It's a (,3) array with units of meters.
        :param back_left_vertex: 1D array
            This is a vector containing the x, y, and z coordinates of the vortex's
            back left point. It's a (,3) array with units of meters.
        :param back_right_vertex: 1D array
            This is a vector containing the x, y, and z coordinates of the vortex's
            back right point. It's a (,3) array with units of meters.
        :param strength: float
            This is the strength of the vortex in meters squared per second.
        """
```

### feature/improved_geometry_definitions (pterasoftware\aerodynamics.py)
```python
        """
        :param Frrvp_G_Cg: array-like of 3 numbers

            Position [x, y, z] of the RingVortex's front-right point (in geometry
            axes, relative to the CG). The front-right point is defined as the end
            point of the RingVortex's right leg and the start point of its front leg.
            Can be a list, tuple, or numpy array of numbers (int or float). Values
            are converted to floats internally. The units are in meters.

        :param Flrvp_G_Cg: array-like of 3 numbers

            Position [x, y, z] of the RingVortex's front-left point (in geometry
            axes, relative to the CG). The front-left point is defined as the end
            point of the RingVortex's front leg and the start point of its left leg.
            Can be a list, tuple, or numpy array of numbers (int or float). Values
            are converted to floats internally. The units are in meters.

        :param Blrvp_G_Cg: array-like of 3 numbers

            Position [x, y, z] of the RingVortex's back-left point (in geometry axes,
            relative to the CG). The back-left point is defined as the end point of
            the RingVortex's left leg and the start point of its back leg. Can be a
            list, tuple, or numpy array of numbers (int or float). Values are
            converted to floats internally. The units are in meters.

        :param Brrvp_G_Cg: array-like of 3 numbers

            Position [x, y, z] of the RingVortex's back-right point (in geometry
            axes, relative to the CG). The back-right point is defined as the end
            point of the RingVortex's back leg and the start point of its right leg.
            Can be a list, tuple, or numpy array of numbers (int or float). Values
            are converted to floats internally. The units are in meters.

        :param strength: number

            This is the strength of the RingVortex It must be a positive number and
            will be converted internally to a float. Its units are in meters squared
            per second.
        """
```

## SteadyHorseshoeVortexLatticeMethodSolver

### main (pterasoftware\steady_horseshoe_vortex_lattice_method.py)
```python
        """
        :param steady_problem: SteadyProblem
            This is the SteadyProblem to be solved.
        :return: None
        """
```

### feature/improved_geometry_definitions (pterasoftware\steady_horseshoe_vortex_lattice_method.py)
```python
        """
        :param steady_problem: SteadyProblem
            This is the SteadyProblem to be solved.
        :return: None
        """
```

## SteadyProblem

### main (pterasoftware\problems.py)
```python
        """
        :param airplanes: list of Airplanes

            This is the list of the Airplanes for this SteadyProblem.

        :param operating_point: OperatingPoint

            This is the OperatingPoint for this SteadyProblem.
        """
```

### feature/improved_geometry_definitions (pterasoftware\problems.py)
```python
        """
        :param airplanes: list of Airplanes

            This is the list of the Airplanes for this SteadyProblem.

        :param operating_point: OperatingPoint

            This is the OperatingPoint for this SteadyProblem.
        """
```

## SteadyRingVortexLatticeMethodSolver

### main (pterasoftware\steady_ring_vortex_lattice_method.py)
```python
        """
        :param steady_problem: SteadyProblem
            This is the SteadyProblem to be solved.
        :return: None
        """
```

### feature/improved_geometry_definitions (pterasoftware\steady_ring_vortex_lattice_method.py)
```python
        """
        :param steady_problem: SteadyProblem
            This is the SteadyProblem to be solved.
        :return: None
        """
```

## UnsteadyProblem

### main (pterasoftware\problems.py)
```python
        """
        :param movement: Movement

            This is the Movement that contains this UnsteadyProblem's
            OperatingPointMovement and AirplaneMovements.

        :param only_final_results: boolLike, optional

            If set to True, the Solver will only calculate forces, moments,
            and pressures for the final complete cycle (of the Movement's
            sub-Movement with the longest period), which increases simulation speed.
            The default value is False.
        """
```

### feature/improved_geometry_definitions (pterasoftware\problems.py)
```python
        """
        :param movement: Movement

            This is the Movement that contains this UnsteadyProblem's
            OperatingPointMovement and AirplaneMovements.

        :param only_final_results: boolLike, optional

            If set to True, the Solver will only calculate forces, moments,
            and pressures for the final complete cycle (of the Movement's
            sub-Movement with the longest period), which increases simulation speed.
            The default value is False.
        """
```

## UnsteadyRingVortexLatticeMethodSolver

### main (pterasoftware\unsteady_ring_vortex_lattice_method.py)
```python
        """
        :param unsteady_problem: UnsteadyProblem
            This is the UnsteadyProblem to be solved.
        :return: None
        """
```

### feature/improved_geometry_definitions (pterasoftware\unsteady_ring_vortex_lattice_method.py)
```python
        """
        :param unsteady_problem: UnsteadyProblem
            This is the UnsteadyProblem to be solved.
        :return: None
        """
```

## Wing

### main (pterasoftware\geometry.py)
```python
        """
        :param wing_cross_sections: list of WingCrossSection objects
            This is a list of WingCrossSection objects, that represent the wing's
            cross sections.
        :param name: str, optional
            This is a sensible name for the wing. The default is "Untitled Wing".
        :param x_le: float, optional
            This is the x coordinate of the leading edge of the wing, relative to the
            airplane's reference point. The default is 0.0.
        :param y_le: float, optional
            This is the y coordinate of the leading edge of the wing, relative to the
            airplane's reference point. The default is 0.0.
        :param z_le: float, optional
            This is the z coordinate of the leading edge of the wing, relative to the
            airplane's reference point. The default is 0.0.
        :param unit_normal_vector: array, optional
            This is an (3,) array of floats that represents the unit normal vector
            of the wing's symmetry plane. It is also the direction vector that the
            wing's span will be assessed relative to. Additionally, this vector
            crossed with the "unit_chordwise_vector" defines the normal vector of the
            plane that the wing's projected area will reference. It must be
            equivalent to this wing's root wing cross section's "unit_normal_vector"
            attribute. The default is np.array([0.0, 1.0, 0.0]), which is the XZ
            plane's unit normal vector.
        :param symmetric: bool, optional
            Set this to true if the wing is across the xz plane. Set it to false if
            not. The default is false.
        :param unit_chordwise_vector: array, optional
            This is an (3,) array of floats that represents the unit vector that
            defines the wing's chordwise direction. This vector crossed with the
            "unit_normal_vector" defines the normal vector of the plane that
            the wing's projected area will reference. This vector must be parallel to
            the intersection of the wing's symmetry plane with each of its wing cross
            section's planes. The default is np.array([1.0, 0.0, 0.0]), which is the
            X unit vector.
        :param num_chordwise_panels: int, optional
            This is the number of chordwise panels to be used on this wing. The
            default is 8.
        :param chordwise_spacing: str, optional
            This is the type of spacing between the wing's chordwise panels. It can
            be set to "cosine" or "uniform". Using a cosine spacing is highly
            recommended for steady simulations and a uniform spacing is highly
            recommended for unsteady simulations. The default is "cosine".
        """
```

### feature/improved_geometry_definitions (pterasoftware\geometry\wing.py)
```python
        """
        :param wing_cross_sections: list of WingCrossSections

            This is a list of WingCrossSections that represent the wing's cross
            sections in order from root to tip. It must contain at least two
            WingCrossSections.

        :param name: str, optional

            This is a sensible name for the Wing. The default is "Untitled Wing".

        :param prelimLer_G_Cg: array-like of 3 numbers, optional

            This is the position [x, y, z] of the origin of this Wing's axes (in
            geometry axes, relative to the starting point) before any symmetry or
            mirror has been applied. Can be a tuple, list, or numpy array of numbers
            (int or float). Values are converted to floats internally. It may differ
            from the actual position as explained in the class docstring. The units
            are meters. The default is (0.0, 0.0, 0.0).

        :param angles_G_to_prelimWn_izyx: array-like of 3 numbers, optional

            This is the rotation angles [angleX, angleY, angleZ] in degrees that
            define the orientation of this Wing's axes relative to the geometry axes
            before any symmetry or mirroring has been applied. Can be a tuple, list,
            or numpy array of numbers (int or float). Values are converted to floats
            internally. All angles must be in the range (-90, 90) degrees. Rotations
            are intrinsic, and proceed in the z-y'-x'' order conventional for Euler
            angles. It may differ from the actual position as explained in the class
            docstring. The units are degrees. The default is (0.0, 0.0, 0.0).

        :param symmetric: boolLike, optional

            Set this to True if the Wing's geometry should be mirrored across the
            symmetry plane while retaining the non-mirrored side. If mirror_only is
            True, symmetric must be False. If symmetric is true, then neither
            symmetry_normal_G nor symmetry_point_G_Cg can be None. If the symmetry
            plane is coincident with this Wing's wing axes' xz-plane, the mirrored
            and non-mirrored geometry will be meshed as a single wing. If not,
            this Wing's Airplane will automatically create another Wing with the
            mirrored geometry, modify both Wings' parameters, and add the reflected
            Wing to its list of wings immediately following this one. For more
            details on how that process, and how this parameter interacts with
            symmetry_normal_G, symmetry_point_G_Cg, and mirror_only, see the class
            docstring. It can be a boolean or a NumPy boolean and will be converted
            internally to a boolean. The default is False.

        :param mirror_only: boolLike, optional

            Set this to True if the Wing's geometry should be reflected about the
            symmetry plane without retaining the non-reflected geometry. If symmetric
            is True, mirror_only must be False. If mirror_only is true, then neither
            symmetry_normal_G nor symmetry_point_G_Cg can be None. For more
            details on how this parameter interacts with symmetry_normal_G,
            symmetry_point_G_Cg, and symmetric, see the class docstring. It can be
            a boolean or a NumPy boolean and will be converted internally to a
            boolean. The default is False.

        :param symmetry_normal_G: array-like of 3 numbers or None, optional

            The unit normal vector (in preliminary wing axes) that, together with
            symmetry_point_G_Cg, defines the plane used for symmetry or mirroring.
            Can be a tuple, list, or numpy array of numbers (int or float), or None.
            Values are converted to floats and normalized internally. Note that
            reversing the normal direction (using the antiparallel vector) defines
            the same plane and produces the same result. This value must be None if
            both symmetric and mirror_only are False, and cannot be None if either
            are True. For more details on how this parameter interacts with
            symmetry_point_G_Cg, symmetric, and mirror_only, see the class
            docstring. The default is None. In mirror_only cases, the wing axes
            themselves are reflected, so this vector is identical in preliminary wing
            axes and wing axes.

        :param symmetry_point_G_Cg: array-like of 3 numbers or None,
        optional

            A point [x, y, z] (in preliminary wing axes, relative to the preliminary
            leading edge root point) that, along with symmetry_normal_G, defines the
            location of the plane about which symmetry or mirroring is applied. Can
            be a list, tuple, or numpy array of numbers (int or float), or None.
            Values are converted to floats internally. This value must be None if
            both symmetric and mirror_only are False, and cannot be None if either
            are True. For more details on how this parameter interacts with
            symmetry_normal_G, symmetric, and mirror_only, see the class docstring.
            The units are meters. The default is None. In mirror_only cases, the wing
            axes themselves are reflected, so this vector is identical in preliminary
            wing axes relative to the preliminary leading edge root point and in wing
            axes relative to the leading edge root point.

        :param num_chordwise_panels: int, optional

            This is the number of chordwise panels to be used on this wing,
            which must be set to a positive integer. The default is 8.

        :param chordwise_spacing: str, optional

            This is the type of spacing between the wing's chordwise panels. It can
            be "cosine" or "uniform". Using cosine spacing is highly recommended for
            steady simulations and uniform spacing is highly recommended for unsteady
            simulations. The default is "cosine".

        """
```

## WingCrossSection

### main (pterasoftware\geometry.py)
```python
        """
        :param airfoil: Airfoil
            This is the airfoil to be used at this wing cross section.
        :param x_le: float, optional
            This is the x coordinate of the leading edge of the wing cross section
            relative to the wing's datum. The default value is 0.0.
        :param y_le: float, optional
            This is the y coordinate of the leading edge of the wing cross section
            relative to the wing's leading edge. The default value is 0.0.
        :param z_le: float, optional
            This is the z coordinate of the leading edge of the wing cross section
            relative to the wing's datum. The default value is 0.0.
        :param chord: float, optional
            This is the chord of the wing at this wing cross section. The default
            value is 1.0.
        :param unit_normal_vector: array, optional
            This is an (3,) array of floats that represents the unit normal vector
            of the plane this wing cross section lies on. If this wing cross section
            is a wing's root, this vector must be equal to the wing's
            unit_normal_vector attribute. Also, every wing cross section
            must have a plane that intersects its parent wing's symmetry plane at a
            line parallel to the parent wing's "unit_chordwise_vector". The default
            is np.array([ 0.0, 1.0, 0.0]), which is the XZ plane's unit normal vector.
        :param twist: float, optional
            This is the twist of the cross section about the leading edge in degrees.
            The default value is 0.0.
        :param control_surface_type: str, optional
            This is type of control surfaces for this wing cross section. It can be
            "symmetric" or "asymmetric". An example of symmetric control surfaces are
            flaps. An example of asymmetric control surfaces are ailerons. The
            default value is "symmetric".
        :param control_surface_hinge_point: float, optional
            This is the location of the control surface hinge from the leading edge
            as a fraction of chord. The default value is 0.75.
        :param control_surface_deflection: float, optional
            This is the control deflection in degrees. Deflection downwards is
            positive. The default value is 0.0 degrees.
        :param num_spanwise_panels: int, optional
            This is the number of spanwise panels to be used between this wing cross
            section and the next one. The default value is 8.
        :param spanwise_spacing: str, optional
            This can be "cosine" or "uniform". Using cosine spacing is highly
            recommended. The default value is "cosine".
        """
```

### feature/improved_geometry_definitions (pterasoftware\geometry\wing_cross_section.py)
```python
        """
        :param airfoil: Airfoil

            This is the Airfoil to be used at this WingCrossSection.

        :param num_spanwise_panels: int or None

            This is the number of spanwise panels to be used between this
            WingCrossSection and the next one. For tip WingCrossSections, this must
            be None. For all other WingCrossSections, this must be a positive integer.

        :param chord: number, optional

            This is the chord of the wing at this WingCrossSection. The units are
            meters. It must be greater than 0.0 and a number (int or float). The
            default value is 1.0.

        :param Lp_Wcsp_Lpp: array-like of 3 numbers, optional

            This is the position [x, y, z] in meters of this WingCrossSection's
            leading edge in parent wing cross section axes, relative to the parent
            leading edge point. Can be a tuple, list, or numpy array of numbers (int
            or float). Values are converted to floats internally. If this is the root
            WingCrossSection, the parent wing cross section axes are the wing axes
            and the parent leading point is the Wing's leading edge root point. If
            not, the parent axes and point are those of the previous
            WingCrossSection. If this is the root WingCrossSection, it must be a zero
            vector. The second component must be non-negative. The default is (0.0,
            0.0, 0.0).

        :param angles_Wcsp_to_Wcs_izyx: array-like of 3 numbers, optional

            This is the angle vector of rotation angles [roll, pitch, yaw] in degrees
            that define the orientation of this WingCrossSection's axes relative to
            the parent wing cross section axes. Can be a tuple, list, or numpy array
            of numbers (int or float). Values are converted to floats internally. If
            this is a root WingCrossSection, these are the wing axes. If not,
            the parent axes are the previous WingCrossSection's axes. For the root
            WingCrossSection, this must be a zero vector. For other
            WingCrossSections, all angles must be in the range (-90, 90) degrees.
            Roll is rotation about the x-axis, pitch is rotation about the y-axis,
            and yaw is rotation about the z-axis. Rotations are intrinsic,
            and proceed in the z-y'-x'' order conventional for Euler angles. The
            units are degrees. The default is (0.0, 0.0, 0.0).

        :param control_surface_symmetry_type: str or None, optional

            This determines how control surfaces behave when the Wing has symmetry.
            It can be "symmetric", "asymmetric", or None. With "symmetric", mirrored
            control surfaces have the same deflection (like flaps). With
            "asymmetric", mirrored control surfaces have opposite deflections (like
            ailerons). The default value is None. For Wings with type 4 or 5 symmetry,
            this parameter must be specified. For Wings with type 1, 2, or 3 symmetry,
            this parameter must be None. This validation is performed by the parent
            Airplane during Wing processing.

        :param control_surface_hinge_point: float, optional

            This is the location of the control surface hinge from the leading edge
            as a fraction of chord. It must be a float the range (0.0, 1.0). The
            default value is 0.75.

        :param control_surface_deflection: number, optional

            This is the control deflection in degrees. Deflection downwards is
            positive. It must be a number (int or float) in the range [-5.0,
            5.0] degrees. The default value is 0.0 degrees.

        :param spanwise_spacing: str or None, optional

            For non-tip WingCrossSections, this can be "cosine" or "uniform". Using
            cosine spacing is highly recommended. For tip WingCrossSections it must
            be None.
        """
```

## WingCrossSectionMovement

### main (pterasoftware\movement.py)
```python
        """
        :param base_wing_cross_section: WingCrossSection
            This is the first wing cross section object, from which the others will
            be created.
        :param sweeping_amplitude: float, optional
            This is the amplitude of the cross section's change in its sweep,
            relative to the vehicle's body axes. Its units are degrees and its
            default value is 0.0 degrees.
        :param sweeping_period: float, optional
            This is the period of the cross section's change in its sweep. Its units
            are seconds and its default value is 0.0 seconds.
        :param sweeping_spacing: string, optional
            This value determines the spacing of the cross section's change in its
            sweep. The options are "sine", "uniform", and "custom". The default value
            is "sine". If "custom", then the value of custom_sweep_function must not
            be none. If both sweeping_spacing and custom_sweep_function are not none,
            then the value of sweeping_spacing will take precedence.
        :param custom_sweep_function: function, optional
            This is a function that describes the motion of the sweeping. For
            example, it could be np.cos or np.sinh (assuming numpy had previously
            been imported as np). It will be horizontally scaled by the
            sweeping_period, vertically scaled by the sweeping_amplitude. For
            example, say the function has an amplitude of 2 units, a period of 3
            units, sweeping_amplitude is set to 4 units and sweeping_period is set to
            5 units. The sweeping motion will have a net amplitude of 8 units and a
            net period of 15 units.
        :param pitching_amplitude: float, optional
            This is the amplitude of the cross section's change in its pitch,
            relative to the vehicle's body axes. Its units are degrees and its
            default value is 0.0 degrees.
        :param pitching_period: float, optional
            This is the period of the cross section's change in its pitch. Its units
            are seconds and its default value is 0.0 seconds.
        :param pitching_spacing: string, optional
            This value determines the spacing of the cross section's change in its
            pitch. The options are "sine", "uniform", and "custom". The default value
            is "sine". If "custom", then the value of custom_pitch_function must not
            be none. If both pitching_spacing and custom_pitch_function are not none,
            then the value of pitching_spacing will take precedence.
        :param custom_pitch_function: function, optional
            This is a function that describes the motion of the pitching. For
            example, it could be np.cos or np.sinh (assuming numpy had previously
            been imported as np). It will be horizontally scaled by the
            pitching_period, vertically scaled by the pitching_amplitude. For
            example, say the function has an amplitude of 2 units, a period of 3
            units, pitching_amplitude is set to 4 units and pitching_period is set to
            5 units. The pitching motion will have a net amplitude of 8 units and a
            net period of 15 units.
        :param heaving_amplitude: float, optional
            This is the amplitude of the cross section's change in its heave,
            relative to the vehicle's body axes. Its units are degrees and its
            default value is 0.0 degrees.
        :param heaving_period: float, optional
            This is the period of the cross section's change in its heave. Its units
            are seconds and its default value is 0.0 seconds.
        :param heaving_spacing: string, optional
            This value determines the spacing of the cross section's change in its
            heave. The options are "sine", "uniform", and "custom". The default value
            is "sine". If "custom", then the value of custom_heave_function must not
            be none. If both heaving_spacing and custom_heave_function are not none,
            then the value of heaving_spacing will take precedence.
        :param custom_heave_function: function, optional
            This is a function that describes the motion of the heaving. For example,
            it could be np.cos or np.sinh (assuming numpy had previously been
        """
```

### feature/improved_geometry_definitions (pterasoftware\movements\wing_cross_section_movement.py)
```python
        """
        :param base_wing_cross_section: WingCrossSection

            This is the base WingCrossSection, from which the WingCrossSection at
            each time step will be created.

        :param ampLp_Wcsp_Lpp: array-like of 3 numbers, optional

            The amplitudes of the WingCrossSectionMovement's changes in its
            WingCrossSections' Lp_Wcsp_Lpp parameters. Can be a tuple, list, or numpy
            array of non-negative numbers (int or float). Values are converted to
            floats internally. The default value is (0.0, 0.0, 0.0). The units are in
            meters.

        :param periodLp_Wcsp_Lpp: array-like of 3 numbers, optional

            The periods of the WingCrossSectionMovement's changes in its
            WingCrossSections' Lp_Wcsp_Lpp parameters. Can be a tuple, list, or numpy
            array of non-negative numbers (int or float). Values are converted to
            floats internally. The default value is (0.0, 0.0, 0.0). Each element
            must be 0.0 if the corresponding element in ampLp_Wcsp_Lpp is 0.0 and
            non-zero if not. The units are in seconds.

        :param spacingLp_Wcsp_Lpp: array-like of 3 strs or callables, optional

            The value determines the spacing of the WingCrossSectionMovement's change
            in its WingCrossSections' Lp_Wcsp_Lpp parameters. Can be a tuple, list,
            or numpy array. Each element can be the string "sine", the string
            "uniform", or a callable custom spacing function. Custom spacing functions
            are for advanced users and must start at 0, return to 0 after one period
            of 2*pi radians, have zero mean, have amplitude of 1, be periodic, return
            finite values only, and accept a ndarray as input and return a ndarray of
            the same shape. The custom function is scaled by ampLp_Wcsp_Lpp, shifted
            by phaseLp_Wcsp_Lpp, and centered around the base value, with the period
            controlled by periodLp_Wcsp_Lpp. The default value is ("sine", "sine",
            "sine").

        :param phaseLp_Wcsp_Lpp: array-like of 3 numbers, optional

            The phase offsets of the elements in the first time step's
            WingCrossSection's Lp_Wcsp_Lpp parameter relative to the base
            WingCrossSection's Lp_Wcsp_Lpp parameter. Can be a tuple, list, or numpy
            array of non-negative numbers (int or float) in the range [0.0, 360.0).
            Values are converted to floats internally. The default value is (0.0,
            0.0, 0.0). Each element must be 0.0 if the corresponding element in
            ampLp_Wcsp_Lpp is 0.0 and non-zero if not. The units are in degrees.

        :param ampAngles_Wcsp_to_Wcs_izyx: array-like of 3 numbers, optional

            The amplitudes of the WingCrossSectionMovement's changes in its
            WingCrossSections' angles_Wcsp_to_Wcs_izyx parameters. Can be a tuple,
            list, or numpy array of numbers (int or float) in the range [0.0,
            180.0). Values are converted to floats internally. The default value is (
            0.0, 0.0, 0.0). The units are in degrees.

        :param periodAngles_Wcsp_to_Wcs_izyx: array-like of 3 numbers, optional

            The periods of the WingCrossSectionMovement's changes in its
            WingCrossSections' angles_Wcsp_to_Wcs_izyx parameters. Can be a tuple,
            list, or numpy array of non-negative numbers (int or float). Values are
            converted to floats internally. The default value is (0.0, 0.0,
            0.0). Each element must be 0.0 if the corresponding element in
            ampAngles_Wcsp_to_Wcs_izyx is 0.0 and non-zero if not. The units are in
            seconds.

        :param spacingAngles_Wcsp_to_Wcs_izyx: array-like of 3 strs or callables, optional

            The value determines the spacing of the WingCrossSectionMovement's change
            in its WingCrossSections' angles_Wcsp_to_Wcs_izyx parameters. Can be a
            tuple, list, or numpy array. Each element can be the string "sine", the
            string "uniform", or a callable custom spacing function. Custom spacing
            functions are for advanced users and must start at 0, return to 0 after
            one period of 2*pi radians, have zero mean, have amplitude of 1, be
            periodic, return finite values only, and accept a ndarray as input and
            return a ndarray of the same shape. The custom function is scaled by
            ampAngles_Wcsp_to_Wcs_izyx, shifted by phaseAngles_Wcsp_to_Wcs_izyx, and
            centered around the base value, with the period controlled by
            periodAngles_Wcsp_to_Wcs_izyx. The default value is ("sine", "sine",
            "sine").

        :param phaseAngles_Wcsp_to_Wcs_izyx: array-like of 3 numbers, optional

            The phase offsets of the elements in the first time step's
            WingCrossSection's angles_Wcsp_to_Wcs_izyx parameter relative to the base
            WingCrossSection's angles_Wcsp_to_Wcs_izyx parameter. Can be a tuple,
            list, or numpy array of numbers (int or float) in the range [0.0,
            360.0). Values are converted to floats internally. The default value is (
            0.0, 0.0, 0.0). Each element must be 0.0 if the corresponding element in
            ampAngles_Wcsp_to_Wcs_izyx is 0.0 and non-zero if not. The units are in
            degrees.
        """
```

## WingMovement

### main (pterasoftware\movement.py)
```python
        """
        :param base_wing: Wing
            This is the first wing object, from which the others will be created.
        :param wing_cross_sections_movements: list of WingCrossSectionMovement objects
            This is a list of the WingCrossSectionMovement objects associated with
            each of the base wing's cross sections.
        :param x_le_amplitude: float, optional
            This is the amplitude of the wing's change in its x reference point. Its
            units are meters and its default value is 0 meters.
        :param x_le_period: float, optional
            This is the period of the wing's change in its x reference point. Its
            units are seconds and its default value is 0 seconds.
        :param x_le_spacing: string, optional
            This value determines the spacing of the wing's change in its x reference
            point. The options are "sine", and "uniform". The default value is "sine".
        :param y_le_amplitude: float, optional
            This is the amplitude of the wing's change in its y reference point. Its
            units are meters and its default value is 0 meters.
        :param y_le_period: float, optional
            This is the period of the wing's change in its y reference point. Its
            units are seconds and its default value is 0 seconds.
        :param y_le_spacing: string, optional
            This value determines the spacing of the wing's change in its y reference
            point. The options are "sine", and "uniform". The default value is "sine".
        :param z_le_amplitude: float, optional
            This is the amplitude of the wing's change in its z reference point. Its
            units are meters and its default value is 0 meters.
        :param z_le_period: float, optional
            This is the period of the wing's change in its z reference point. Its
            units are seconds and its default value is 0 seconds.
        :param z_le_spacing: string, optional
            This value determines the spacing of the wing's change in its z reference
            point. The options are "sine", and "uniform". The default value is "sine".
        """
```

### feature/improved_geometry_definitions (pterasoftware\movements\wing_movement.py)
```python
        """
        :param base_wing: Wing

            This is the base Wing, from which the Wing at each time step will be
            created.

        :param wing_cross_section_movements: list of WingCrossSectionMovements

            This is a list of the WingCrossSectionMovements associated with each of
            the base Wing's WingCrossSections. It must have the same length as the
            base Wing's list of WingCrossSections.

        :param ampPrelimLer_G_Cg: array-like of 3 numbers, optional

            The amplitudes of the WingMovement's changes in its Wings' prelimLer_G_Cg
            parameters. Can be a tuple, list, or numpy array of non-negative numbers
            (int or float). Values are converted to floats internally. The default
            value is (0.0, 0.0, 0.0). The units are in meters.

        :param periodPrelimLer_G_Cg: array-like of 3 numbers, optional

            The periods of the WingMovement's changes in its Wings' prelimLer_G_Cg
            parameters. Can be a tuple, list, or numpy array of non-negative numbers
            (int or float). Values are converted to floats internally. The default
            value is (0.0, 0.0, 0.0). Each element must be 0.0 if the corresponding
            element in ampPrelimLer_G_Cg is 0.0 and non-zero if not. The units are in
            seconds.

        :param spacingPrelimLer_G_Cg: array-like of 3 strs or callables, optional

            The value determines the spacing of the WingMovement's change in its
            Wings' prelimLer_G_Cg parameters. Can be a tuple, list, or numpy array.
            Each element can be the string "sine", the string "uniform", or a callable
            custom spacing function. Custom spacing functions are for advanced users
            and must start at 0, return to 0 after one period of 2*pi radians, have
            zero mean, have amplitude of 1, be periodic, return finite values only,
            and accept a ndarray as input and return a ndarray of the same shape. The
            custom function is scaled by ampPrelimLer_G_Cg, shifted by
            phasePrelimLer_G_Cg, and centered around the base value, with the period
            controlled by periodPrelimLer_G_Cg. The default value is ("sine", "sine",
            "sine").

        :param phasePrelimLer_G_Cg: array-like of 3 numbers, optional

            The phase offsets of the elements in the first time step's Wing's
            prelimLer_G_Cg parameter relative to the base Wing's prelimLer_G_Cg
            parameter. Can be a tuple, list, or numpy array of non-negative numbers (
            int or float) in the range [0.0, 360.0). Values are converted to floats
            internally. The default value is (0.0, 0.0, 0.0). Each element must be
            0.0 if the corresponding element in ampPrelimLer_G_Cg is 0.0 and non-zero
            if not. The units are in degrees.

        :param ampAngles_G_to_prelimWn_izyx: array-like of 3 numbers, optional

            The amplitudes of the WingMovement's changes in its Wings'
            angles_G_to_prelimWn_izyx parameters. Can be a tuple, list, or numpy
            array of numbers (int or float) in the range [0.0, 180.0). Values are
            converted to floats internally. The default value is (0.0, 0.0, 0.0). The
            units are in degrees.

        :param periodAngles_G_to_prelimWn_izyx: array-like of 3 numbers, optional

            The periods of the WingMovement's changes in its Wings'
            angles_G_to_prelimWn_izyx parameters. Can be a tuple, list, or numpy
            array of non-negative numbers (int or float). Values are converted to
            floats internally. The default value is (0.0, 0.0, 0.0). Each element
            must be 0.0 if the corresponding element in ampAngles_G_to_prelimWn_izyx
            is 0.0 and non-zero if not. The units are in seconds.

        :param spacingAngles_G_to_prelimWn_izyx: array-like of 3 strs or callables, optional

            The value determines the spacing of the WingMovement's change in its
            Wings' angles_G_to_prelimWn_izyx parameters. Can be a tuple, list,
            or numpy array. Each element can be the string "sine", the string
            "uniform", or a callable custom spacing function. Custom spacing functions
            are for advanced users and must start at 0, return to 0 after one period
            of 2*pi radians, have zero mean, have amplitude of 1, be periodic, return
            finite values only, and accept a ndarray as input and return a ndarray of
            the same shape. The custom function is scaled by
            ampAngles_G_to_prelimWn_izyx, shifted by phaseAngles_G_to_prelimWn_izyx,
            and centered around the base value, with the period controlled by
            periodAngles_G_to_prelimWn_izyx. The default value is ("sine", "sine",
            "sine").

        :param phaseAngles_G_to_prelimWn_izyx: array-like of 3 numbers, optional

            The phase offsets of the elements in the first time step's Wing's
            angles_G_to_prelimWn_izyx parameter relative to the base Wing's
            angles_G_to_prelimWn_izyx parameter. Can be a tuple, list, or numpy array
            of numbers (int or float) in the range [0.0, 360.0). Values are converted
            to floats internally. The default value is (0.0, 0.0, 0.0). Each element
            must be 0.0 if the corresponding element in ampAngles_G_to_prelimWn_izyx
            is 0.0 and non-zero if not. The units are in degrees.
        """
```