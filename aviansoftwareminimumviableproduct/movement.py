"""

"""

import numpy as np
import aviansoftwareminimumviableproduct as asmvp


class Movement:
    """

    """

    def __init__(self, airplane_movement, operating_point_movement, num_steps, delta_time):
        """

        :param airplane_movement:
        :param operating_point_movement:
        :param num_steps:
        :param delta_time:
        """

        self.airplane_movement = airplane_movement
        self.operating_point_movement = operating_point_movement

        self.airplanes = self.airplane_movement.generate_airplanes(
            num_steps=num_steps,
            delta_time=delta_time
        )
        self.operating_points = self.operating_point_movement.generate_operating_points(
            num_steps=num_steps,
            delta_time=delta_time
        )


class AirplaneMovement:
    """

    """

    def __init__(self, base_airplane, wing_movements):
        """

        :param base_airplane:
        :param wing_movements:
        """

        self.base_airplane = base_airplane
        self.wing_movements = wing_movements

    def generate_airplanes(self, num_steps, delta_time):
        """

        :param num_steps:
        :param delta_time:
        :return:
        """

        wings = np.empty((len(self.wing_movements), num_steps), dtype=object)

        for wing_movement_location, wing_movement in np.ndenumerate(self.wing_movements):
            wing_movements_list_of_wings = wing_movement.generate_wings(num_steps=num_steps, delta_time=delta_time)
            wing_movements_array_of_wings = np.array(wing_movements_list_of_wings)
            wings[wing_movement_location, :] = wing_movements_array_of_wings

        airplanes = []

        name = self.base_airplane.name
        x_ref = self.base_airplane.x_ref
        y_ref = self.base_airplane.y_ref
        z_ref = self.base_airplane.z_ref

        for step in range(num_steps):
            wings_at_time_step = wings[:, step]

            this_airplane = asmvp.geometry.Airplane(
                name=name,
                x_ref=x_ref,
                y_ref=y_ref,
                z_ref=z_ref,
                wings=wings_at_time_step
            )

            airplanes.append(this_airplane)

        return airplanes


class WingMovement:
    """
    
    """
    
    def __init__(self, base_wing,
                 sweeping_amplitude, sweeping_period, sweeping_spacing, sweeping_base,
                 heaving_amplitude, heaving_period, heaving_spacing, heaving_base,
                 pitching_amplitude, pitching_period, pitching_spacing, pitching_base):
        """
        
        :param base_wing: 
        :param sweeping_amplitude: 
        :param sweeping_period: 
        :param sweeping_spacing:
        :param sweeping_base: 
        :param heaving_amplitude: 
        :param heaving_period: 
        :param heaving_spacing:
        :param heaving_base: 
        :param pitching_amplitude: 
        :param pitching_period: 
        :param pitching_spacing:
        :param pitching_base: 
        """
        
        self.base_wing = base_wing

        self.sweeping_amplitude = sweeping_amplitude
        self.sweeping_period = sweeping_period
        self.sweeping_spacing = sweeping_spacing
        self.sweeping_base = sweeping_base
        
        self.heaving_amplitude = heaving_amplitude
        self.heaving_period = heaving_period
        self.heaving_spacing = heaving_spacing
        self.heaving_base = heaving_base
        
        self.pitching_amplitude = pitching_amplitude
        self.pitching_period = pitching_period
        self.pitching_spacing = pitching_spacing
        self.pitching_base = pitching_base

    def generate_wings(self, num_steps, delta_time):
        """

        :param num_steps:
        :param delta_time:
        :return:
        """

        if self.sweeping_spacing == 'sine':
            sweeping_values = asmvp.performance.oscillating_sinspace(
                amplitude=self.sweeping_amplitude,
                period=self.sweeping_period,
                base_value=self.sweeping_base,
                num_steps=num_steps,
                delta_time=delta_time
            )
        elif self.sweeping_spacing == 'uniform':
            sweeping_values = asmvp.performance.oscillating_linspace(
                amplitude=self.sweeping_amplitude,
                period=self.sweeping_period,
                base_value=self.sweeping_base,
                num_steps=num_steps,
                delta_time=delta_time
            )
        else:
            raise Exception("Bad value of sweeping_spacing!")

        if self.heaving_spacing == 'sine':
            heaving_values = asmvp.performance.oscillating_sinspace(
                amplitude=self.heaving_amplitude,
                period=self.heaving_period,
                base_value=self.heaving_base,
                num_steps=num_steps,
                delta_time=delta_time
            )
        elif self.heaving_spacing == 'uniform':
            heaving_values = asmvp.performance.oscillating_linspace(
                amplitude=self.heaving_amplitude,
                period=self.heaving_period,
                base_value=self.heaving_base,
                num_steps=num_steps,
                delta_time=delta_time
            )
        else:
            raise Exception("Bad value of heaving_spacing!")

        if self.pitching_spacing == 'sine':
            pitching_values = asmvp.performance.oscillating_sinspace(
                amplitude=self.pitching_amplitude,
                period=self.pitching_period,
                base_value=self.pitching_base,
                num_steps=num_steps,
                delta_time=delta_time
            )
        elif self.pitching_spacing == 'uniform':
            pitching_values = asmvp.performance.oscillating_linspace(
                amplitude=self.pitching_amplitude,
                period=self.pitching_period,
                base_value=self.pitching_base,
                num_steps=num_steps,
                delta_time=delta_time
            )
        else:
            raise Exception("Bad value of pitching_spacing!")

        for panel in self.base_wing.panels:
            front_left_vertex = 

        wing_cross_sections = np.empty((len(self.wing_cross_section_movements), num_steps), dtype=object)

        for wing_cross_section_movement in self.wing_cross_section_movements:


        for wing_cross_section_movement_location in range(len(self.wing_cross_section_movements)):
            wing_cross_section_movement = self.wing_cross_section_movements[wing_cross_section_movement_location]
            x = np.array(
                wing_cross_section_movement.generate_wing_cross_sections(num_steps=num_steps, delta_time=delta_time))
            wing_cross_sections[wing_cross_section_movement_location, :] = x

        # Create an empty list of wings.
        wings = []

        name = self.base_wing.name
        symmetric = self.base_wing.symmetric
        num_chordwise_panels = self.base_wing.num_chordwise_panels
        chordwise_spacing = self.base_wing.chordwise_spacing

        for step in range(num_steps):
            x_le = x_le_list[step]
            y_le = y_le_list[step]
            z_le = z_le_list[step]
            cross_sections = wing_cross_sections[:, step]

            this_wing = asmvp.geometry.Wing(
                name=name,
                x_le=x_le,
                y_le=y_le,
                z_le=z_le,
                wing_cross_sections=cross_sections,
                symmetric=symmetric,
                num_chordwise_panels=num_chordwise_panels,
                chordwise_spacing=chordwise_spacing
            )

            wings.append(this_wing)

        # Return the ndarray of wing cross sections.
        return wings


class OperatingPointMovement:
    """
    
    """
    
    def __init__(self, base_operating_point,
                 velocity_amplitude, velocity_period, velocity_spacing,
                 alpha_amplitude, alpha_period, alpha_spacing,
                 beta_amplitude, beta_period, beta_spacing):
        """
        
        :param base_operating_point: 
        :param velocity_amplitude: 
        :param velocity_period: 
        :param velocity_spacing: 
        :param alpha_amplitude: 
        :param alpha_period: 
        :param alpha_spacing: 
        :param beta_amplitude: 
        :param beta_period: 
        :param beta_spacing: 
        """
        
        self.base_operating_point = base_operating_point
        
        self.base_velocity = self.base_operating_point.velocity
        self.velocity_amplitude = velocity_amplitude
        self.velocity_period = velocity_period
        self.velocity_spacing = velocity_spacing
        
        self.base_alpha = self.base_operating_point.alpha
        self.alpha_amplitude = alpha_amplitude
        self.alpha_period = alpha_period
        self.alpha_spacing = alpha_spacing
        
        self.base_beta = self.base_operating_point.beta
        self.beta_amplitude = beta_amplitude
        self.beta_period = beta_period
        self.beta_spacing = beta_spacing

    def generate_operating_points(self, num_steps, delta_time):
        """

        :param num_steps:
        :param delta_time:
        :return:
        """

        pass
