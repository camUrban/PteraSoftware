"""This is script is an example of how to run Ptera Software's unsteady ring vortex
lattice method solver on a simple airplane object."""

import pterasoftware as ps

# import numpy as np

size = 4

simple_airplane = ps.geometry.Airplane(
    wings=[
        ps.geometry.Wing(
            wing_cross_sections=[
                ps.geometry.WingCrossSection(
                    airfoil=ps.geometry.Airfoil(
                        name="naca0012",
                        coordinates=None,
                        repanel=True,
                        n_points_per_side=400,
                    ),
                    x_le=0.0,
                    y_le=0.0,
                    z_le=0.0,
                    chord=2.0,
                    # unit_normal_vector=np.array([0.0, 1.0, 0.0]),
                    twist=0.0,
                    control_surface_type="symmetric",
                    control_surface_hinge_point=0.75,
                    control_surface_deflection=0.0,
                    num_spanwise_panels=int(size / 2),
                    spanwise_spacing="cosine",
                ),
                ps.geometry.WingCrossSection(
                    airfoil=ps.geometry.Airfoil(
                        name="naca0012",
                        coordinates=None,
                        repanel=True,
                        n_points_per_side=400,
                    ),
                    x_le=0.0,
                    y_le=1.0,
                    z_le=0.0,
                    chord=2.0,
                    # unit_normal_vector=np.array([0.0, 1.0, 0.0]),
                    twist=0.0,
                    control_surface_type="symmetric",
                    control_surface_hinge_point=0.75,
                    control_surface_deflection=0.0,
                    num_spanwise_panels=int(size / 2),
                    spanwise_spacing="cosine",
                ),
            ],
            name="Main Wing",
            x_le=-1.0,
            y_le=0.0,
            z_le=0.0,
            # symmetry_unit_normal_vector=np.array([0.0, 1.0, 0.0]),
            symmetric=True,
            # unit_chordwise_vector=np.array([1.0, 0.0, 0.0]),
            num_chordwise_panels=size,
            chordwise_spacing="cosine",
        ),
    ],
    name="Simple Airplane",
    x_ref=0.0,
    y_ref=0.0,
    z_ref=0.0,
    weight=0.0,
    s_ref=None,
    c_ref=None,
    b_ref=None,
)

wing = simple_airplane.wings[0]
panels = wing.panels
show_z = False

[num_rows, num_cols] = panels.shape

p_text = ""
row_num = 1
for row_of_panels in panels:
    col_num = 1
    for panel in row_of_panels:
        p_text += "│"
        col_num += 1

        [this_fl_x, this_fl_y, this_fl_z] = panel.front_left_vertex
        [this_fr_x, this_fr_y, this_fr_z] = panel.front_right_vertex

        f_text = "({: .5e}, {: .5e}"
        if show_z:
            f_text += ", {: .5e}"
        f_text += ") "
        f_text += "({: .5e}, {: .5e}"
        if show_z:
            f_text += ", {: .5e}"
        f_text += ")"

        if show_z:
            this_f_text = f_text.format(
                this_fl_x,
                this_fl_y,
                this_fl_z,
                this_fr_x,
                this_fr_y,
                this_fr_z,
            )
        else:
            this_f_text = f_text.format(
                this_fl_x,
                this_fl_y,
                this_fr_x,
                this_fr_y,
            )

        p_text += this_f_text

    p_text += "│\n"
    for panel in row_of_panels:
        p_text += "│"

        [this_bl_x, this_bl_y, this_bl_z] = panel.back_left_vertex
        [this_br_x, this_br_y, this_br_z] = panel.back_right_vertex

        b_text = "({: .5e}, {: .5e}"
        if show_z:
            b_text += ", {: .5e}"
        b_text += ") "
        b_text += "({: .5e}, {: .5e}"
        if show_z:
            b_text += ", {: .5e}"
        b_text += ")"

        if show_z:
            this_b_text = b_text.format(
                this_bl_x,
                this_bl_y,
                this_bl_z,
                this_br_x,
                this_br_y,
                this_br_z,
            )
        else:
            this_b_text = b_text.format(
                this_bl_x,
                this_bl_y,
                this_br_x,
                this_br_y,
            )

        p_text += this_b_text
    p_text += "│"
    last_line = p_text.splitlines()[-1]

    col_len = int((len(last_line) - 1) / num_cols)

    if row_num != num_rows:
        col_text = "─" * (col_len - 1) + "┼"
    else:
        col_text = "─" * (col_len - 1) + "┴"

    all_col_text = col_text * num_cols
    all_col_text = all_col_text[:-1]

    if row_num != num_rows:
        all_col_text = "├" + all_col_text + "┤"
    else:
        all_col_text = "└" + all_col_text + "┘"

    p_text = p_text + "\n" + all_col_text + "\n"

    row_num += 1

p_text = p_text[:-1]

last_line = p_text.splitlines()[-1]
col_len = int((len(last_line) - 1) / num_cols)
col_text = "─" * (col_len - 1) + "┬"
all_col_text = col_text * num_cols
all_col_text = all_col_text[:-1]
all_col_text = "┌" + all_col_text + "┐\n"

p_text = all_col_text + p_text

print(p_text)
