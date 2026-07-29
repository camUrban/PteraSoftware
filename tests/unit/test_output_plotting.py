"""This module contains classes to test the output plotting functions."""

import csv
import tempfile
import unittest
from pathlib import Path

import matplotlib.colors
import matplotlib.image
import matplotlib.pyplot as plt
import numpy as np
import numpy.testing as npt

# noinspection PyProtectedMember
from pterasoftware import _output_plotting, _output_rendering, _transformations
from tests.unit.fixtures import operating_point_fixtures, output_plotting_fixtures


class TestTextColorNormalized(unittest.TestCase):
    """This class contains methods for testing the plots' normalized text color."""

    def test_matches_the_rendered_visualizations_text_color(self) -> None:
        """Test that the plots' text color is the rendered visualizations' text color.

        The two kinds of output are meant to look related, which they only do while the
        plots' normalized color stays the rendered color divided by 255.0.
        """
        npt.assert_allclose(
            _output_plotting._TEXT_COLOR_NORMALIZED,
            np.array(_output_rendering.TEXT_COLOR, dtype=float) / 255.0,
        )


class TestGetOperatingPointVelocity(unittest.TestCase):
    """This class contains methods for testing
    _output_plotting.get_operating_point_velocity."""

    def test_returns_a_three_component_vector(self) -> None:
        """Test that the CG velocity is returned as a (3,) ndarray of floats.

        The CG velocity is the first Airplane's CG velocity (in Earth axes, observed
        from the Earth frame), which every test in this class refers to by that shorter
        name.
        """
        operating_point = operating_point_fixtures.make_basic_operating_point_fixture()
        velocity_E__E = _output_plotting.get_operating_point_velocity(operating_point)
        self.assertEqual(velocity_E__E.shape, (3,))
        self.assertEqual(velocity_E__E.dtype, float)

    def test_alpha_alone_does_not_tilt_wind_aligned_earth_axes(self) -> None:
        """Test that a large alpha leaves the CG velocity along the Earth x axis.

        An OperatingPoint given no attitude angles derives them from alpha and beta,
        which aligns Earth axes with wind axes. A large alpha then tilts Earth axes
        along with the first Airplane, which travels along the Earth x axis at its
        speed.
        """
        operating_point = (
            operating_point_fixtures.make_high_alpha_operating_point_fixture()
        )
        velocity_E__E = _output_plotting.get_operating_point_velocity(operating_point)
        npt.assert_allclose(velocity_E__E, [10.0, 0.0, 0.0], atol=1e-12)

    def test_alpha_tilts_the_velocity_of_a_level_airplane(self) -> None:
        """Test that alpha tilts the CG velocity of an airplane held level.

        With zero attitude angles, Earth axes and the first Airplane's body axes
        coincide, so an alpha of 5.0 degrees puts the first Airplane's nose 5.0 degrees
        above its flight path. Earth axes have +z pointing down, so the CG velocity has
        a positive z component.
        """
        operating_point = operating_point_fixtures.make_basic_operating_point_fixture()
        velocity_E__E = _output_plotting.get_operating_point_velocity(operating_point)
        alpha_rad = np.deg2rad(5.0)
        npt.assert_allclose(
            velocity_E__E,
            [10.0 * np.cos(alpha_rad), 0.0, 10.0 * np.sin(alpha_rad)],
            atol=1e-12,
        )

    def test_is_the_negative_of_the_freestream_velocity(self) -> None:
        """Test that the CG velocity is the negative of the freestream velocity.

        Rotating the CG velocity back into the first Airplane's geometry axes must
        recover the freestream velocity the OperatingPoint stores, negated.
        """
        operating_point = (
            operating_point_fixtures.make_with_attitude_angles_operating_point_fixture()
        )
        velocity_E__E = _output_plotting.get_operating_point_velocity(operating_point)
        velocity_GP1__E = _transformations.apply_T_to_vectors(
            operating_point.T_pas_E_CgP1_to_GP1_CgP1,
            velocity_E__E,
            is_position=False,
        )
        npt.assert_allclose(velocity_GP1__E, -operating_point.vInf_GP1__E, atol=1e-12)

    def test_scales_with_the_operating_points_speed(self) -> None:
        """Test that a faster OperatingPoint gives a faster CG velocity."""
        operating_point = (
            operating_point_fixtures.make_high_speed_operating_point_fixture()
        )
        velocity_E__E = _output_plotting.get_operating_point_velocity(operating_point)
        npt.assert_allclose(velocity_E__E, [100.0, 0.0, 0.0], atol=1e-12)

    def test_does_not_share_memory_with_the_operating_point(self) -> None:
        """Test that the CG velocity is a new array rather than a view of a cached one.

        The OperatingPoint's cached freestream velocity is read-only, so a caller that
        writes to the returned array must not reach it.
        """
        operating_point = operating_point_fixtures.make_basic_operating_point_fixture()
        vInf_GP1__E = operating_point.vInf_GP1__E.copy()
        velocity_E__E = _output_plotting.get_operating_point_velocity(operating_point)
        velocity_E__E[0] = 0.0
        npt.assert_array_equal(operating_point.vInf_GP1__E, vInf_GP1__E)


class TestCsvHeaders(unittest.TestCase):
    """This class contains methods for testing _output_plotting.csv_headers."""

    def test_appends_the_context_and_the_unit(self) -> None:
        """Test that a header is its label, then its context, then its unit."""
        headers = _output_plotting.csv_headers(
            ["Induced Drag", "Side Force", "Lift"], "(in Wind Axes)", "Force (N)"
        )
        self.assertEqual(
            headers,
            [
                "Induced Drag in Wind Axes (N)",
                "Side Force in Wind Axes (N)",
                "Lift in Wind Axes (N)",
            ],
        )

    def test_drops_the_subtitles_parentheses_and_commas(self) -> None:
        """Test that a subtitle loses its parentheses and commas on its way to a header.

        A subtitle is parenthesized and comma separated, and a header is neither.
        """
        headers = _output_plotting.csv_headers(
            ["Rolling Moment"], "(in Wind Axes, Relative to the CG)", "Moment (N m)"
        )
        self.assertEqual(
            headers, ["Rolling Moment in Wind Axes Relative to the CG (N*m)"]
        )

    def test_omits_the_unit_when_the_y_label_has_none(self) -> None:
        """Test that a header carries no unit when the y-axis label names none."""
        headers = _output_plotting.csv_headers(
            ["Induced Drag Coefficient"], "(in Wind Axes)", "Force Coefficient"
        )
        self.assertEqual(headers, ["Induced Drag Coefficient in Wind Axes"])

    def test_omits_the_context_when_the_subtitle_is_empty(self) -> None:
        """Test that a header carries no context when the figure has no subtitle."""
        headers = _output_plotting.csv_headers(
            ["Angle of Attack", "Sideslip Angle"], "", "Angle (deg)"
        )
        self.assertEqual(headers, ["Angle of Attack (deg)", "Sideslip Angle (deg)"])

    def test_substitutes_the_quantity_for_a_component_label(self) -> None:
        """Test that a component label takes its quantity from the y-axis label.

        A figure that labels its series by component alone names the quantity in its
        title, which a column has no room for.
        """
        headers = _output_plotting.csv_headers(
            ["X Component", "Y Component", "Z Component"],
            "(of the First Airplane's CG, in Earth Axes, Relative to the Earth Origin)",
            "Position (m)",
        )
        self.assertEqual(
            headers,
            [
                "X Position of the First Airplane's CG in Earth Axes Relative to the "
                "Earth Origin (m)",
                "Y Position of the First Airplane's CG in Earth Axes Relative to the "
                "Earth Origin (m)",
                "Z Position of the First Airplane's CG in Earth Axes Relative to the "
                "Earth Origin (m)",
            ],
        )

    def test_leaves_a_label_that_is_only_component_alone(self) -> None:
        """Test that only a label ending in " Component" takes on the quantity."""
        headers = _output_plotting.csv_headers(["Component"], "", "Position (m)")
        self.assertEqual(headers, ["Component (m)"])

    def test_returns_no_headers_for_no_labels(self) -> None:
        """Test that a figure with no series composes no headers."""
        self.assertEqual(
            _output_plotting.csv_headers([], "(in Wind Axes)", "Force (N)"), []
        )

    def test_does_not_modify_the_labels(self) -> None:
        """Test that composing headers leaves the figure's legend labels alone."""
        labels = ["X Component", "Y Component", "Z Component"]
        _output_plotting.csv_headers(labels, "(in Earth Axes)", "Position (m)")
        self.assertEqual(labels, ["X Component", "Y Component", "Z Component"])


class TestWriteTimeHistoryCsv(unittest.TestCase):
    """This class contains methods for testing
    _output_plotting.write_time_history_csv."""

    def setUp(self) -> None:
        """Set up a temporary directory to write each test's CSV file to."""
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary_directory.cleanup)
        self.save_path = Path(self.temporary_directory.name) / "time_history.csv"

    def _read_rows(self) -> list[list[str]]:
        """Reads back the rows of the CSV file that a test wrote.

        :return: A list of the file's rows, each a list of its cells' text.
        """
        with open(self.save_path, "r", newline="", encoding="utf-8") as csv_file:
            return list(csv.reader(csv_file))

    def test_writes_a_header_row_led_by_time(self) -> None:
        """Test that the first row is the time column's header, then the others."""
        times = output_plotting_fixtures.make_times_fixture()
        headers = output_plotting_fixtures.make_three_headers_fixture()
        columns = output_plotting_fixtures.make_three_series_fixture()
        _output_plotting.write_time_history_csv(times, headers, columns, self.save_path)
        self.assertEqual(self._read_rows()[0], ["Time (s)"] + headers)

    def test_writes_each_time_step_in_order(self) -> None:
        """Test that each row leads with its time step's time, in order."""
        times = output_plotting_fixtures.make_times_fixture()
        headers = output_plotting_fixtures.make_three_headers_fixture()
        columns = output_plotting_fixtures.make_three_series_fixture()
        _output_plotting.write_time_history_csv(times, headers, columns, self.save_path)
        written_times = [float(row[0]) for row in self._read_rows()[1:]]
        npt.assert_allclose(written_times, times)

    def test_writes_each_column_in_the_order_it_was_given(self) -> None:
        """Test that each column's values land in that column's cells."""
        times = output_plotting_fixtures.make_times_fixture()
        headers = output_plotting_fixtures.make_three_headers_fixture()
        columns = output_plotting_fixtures.make_three_series_fixture()
        _output_plotting.write_time_history_csv(times, headers, columns, self.save_path)
        rows = self._read_rows()[1:]
        for column_id, column in enumerate(columns):
            written_column = [float(row[column_id + 1]) for row in rows]
            npt.assert_allclose(written_column, column)

    def test_writes_only_the_time_column_for_no_columns(self) -> None:
        """Test that a time history with no columns still writes its times."""
        times = output_plotting_fixtures.make_times_fixture()
        _output_plotting.write_time_history_csv(times, [], [], self.save_path)
        rows = self._read_rows()
        self.assertEqual(rows[0], ["Time (s)"])
        self.assertEqual(len(rows), len(times) + 1)

    def test_writes_only_the_header_row_for_no_time_steps(self) -> None:
        """Test that a time history with no time steps writes its headers alone."""
        headers = output_plotting_fixtures.make_three_headers_fixture()
        columns = [np.zeros(0, dtype=float) for _ in headers]
        _output_plotting.write_time_history_csv(
            np.zeros(0, dtype=float), headers, columns, self.save_path
        )
        self.assertEqual(self._read_rows(), [["Time (s)"] + headers])

    def test_raises_for_more_columns_than_headers(self) -> None:
        """Test that more columns than headers is rejected."""
        times = output_plotting_fixtures.make_times_fixture()
        headers = output_plotting_fixtures.make_three_headers_fixture()[:1]
        columns = output_plotting_fixtures.make_three_series_fixture()
        with self.assertRaises(ValueError):
            _output_plotting.write_time_history_csv(
                times, headers, columns, self.save_path
            )

    def test_writes_nothing_when_it_rejects_its_input(self) -> None:
        """Test that a rejected time history leaves no partly written file behind."""
        times = output_plotting_fixtures.make_times_fixture()
        headers = output_plotting_fixtures.make_three_headers_fixture()
        columns = output_plotting_fixtures.make_single_series_fixture()
        with self.assertRaises(ValueError):
            _output_plotting.write_time_history_csv(
                times, headers, columns, self.save_path
            )
        self.assertFalse(self.save_path.exists())

    def test_ends_its_lines_with_a_line_feed(self) -> None:
        """Test that the file's lines end with a line feed alone.

        The csv module writes a carriage return and a line feed unless it is told
        otherwise. No other text file this project writes ends its lines that way.
        """
        times = output_plotting_fixtures.make_times_fixture()
        headers = output_plotting_fixtures.make_three_headers_fixture()
        columns = output_plotting_fixtures.make_three_series_fixture()
        _output_plotting.write_time_history_csv(times, headers, columns, self.save_path)
        contents = self.save_path.read_bytes()
        self.assertNotIn(b"\r", contents)
        self.assertTrue(contents.endswith(b"\n"))

    def test_overwrites_an_existing_file(self) -> None:
        """Test that writing over an earlier run's file leaves none of it behind."""
        times = output_plotting_fixtures.make_times_fixture()
        headers = output_plotting_fixtures.make_three_headers_fixture()
        columns = output_plotting_fixtures.make_three_series_fixture()
        _output_plotting.write_time_history_csv(times, headers, columns, self.save_path)
        _output_plotting.write_time_history_csv(
            times[:2], headers[:1], columns[:1], self.save_path
        )
        rows = self._read_rows()
        self.assertEqual(rows[0], ["Time (s)", headers[0]])
        self.assertEqual(len(rows), 3)


class TestPlotTimeHistory(unittest.TestCase):
    """This class contains methods for testing _output_plotting.plot_time_history."""

    def setUp(self) -> None:
        """Set up a temporary directory to save each test's figure to."""
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary_directory.cleanup)
        self.save_path = Path(self.temporary_directory.name) / "time_history.png"

    def tearDown(self) -> None:
        """Close the figures that the tested function leaves open."""
        plt.close("all")

    def test_draws_one_line_per_series(self) -> None:
        """Test that each series is drawn as its own line, in the order it was given."""
        _output_plotting.plot_time_history(
            output_plotting_fixtures.make_times_fixture(),
            output_plotting_fixtures.make_three_series_fixture(),
            output_plotting_fixtures.make_three_labels_fixture(),
            output_plotting_fixtures.make_three_colors_fixture(),
            "Example Airplane Forces",
            "(in Wind Axes)",
            "Force (N)",
            output_plotting_fixtures.make_figure_size_fixture(),
            False,
            self.save_path,
            300.0,
        )
        axes = plt.gcf().axes[0]
        series = output_plotting_fixtures.make_three_series_fixture()
        self.assertEqual(len(axes.lines), len(series))
        for line, this_series in zip(axes.lines, series):
            npt.assert_allclose(
                np.asarray(line.get_xdata(), dtype=float),
                output_plotting_fixtures.make_times_fixture(),
            )
            npt.assert_allclose(np.asarray(line.get_ydata(), dtype=float), this_series)

    def test_labels_and_colors_each_line(self) -> None:
        """Test that each line takes the label and color given for its series."""
        labels = output_plotting_fixtures.make_three_labels_fixture()
        colors = output_plotting_fixtures.make_three_colors_fixture()
        _output_plotting.plot_time_history(
            output_plotting_fixtures.make_times_fixture(),
            output_plotting_fixtures.make_three_series_fixture(),
            labels,
            colors,
            "Example Airplane Forces",
            "(in Wind Axes)",
            "Force (N)",
            output_plotting_fixtures.make_figure_size_fixture(),
            False,
            self.save_path,
            300.0,
        )
        axes = plt.gcf().axes[0]
        self.assertEqual([line.get_label() for line in axes.lines], labels)
        self.assertEqual([line.get_color() for line in axes.lines], colors)

    def test_draws_the_lines_from_thickest_to_thinnest(self) -> None:
        """Test that the lines thin evenly from the first series to the last.

        Drawing them in this order is what keeps every line visible where the curves
        overlap.
        """
        _output_plotting.plot_time_history(
            output_plotting_fixtures.make_times_fixture(),
            output_plotting_fixtures.make_three_series_fixture(),
            output_plotting_fixtures.make_three_labels_fixture(),
            output_plotting_fixtures.make_three_colors_fixture(),
            "Example Airplane Forces",
            "(in Wind Axes)",
            "Force (N)",
            output_plotting_fixtures.make_figure_size_fixture(),
            False,
            self.save_path,
            300.0,
        )
        axes = plt.gcf().axes[0]
        self.assertEqual([line.get_linewidth() for line in axes.lines], [3.5, 2.5, 1.5])

    def test_draws_a_lone_series_at_the_thickest_line_width(self) -> None:
        """Test that a figure with one series draws it at the thickest line width."""
        _output_plotting.plot_time_history(
            output_plotting_fixtures.make_times_fixture(),
            output_plotting_fixtures.make_single_series_fixture(),
            ["Lift"],
            ["tab:blue"],
            "Example Airplane Forces",
            "(in Wind Axes)",
            "Force (N)",
            output_plotting_fixtures.make_figure_size_fixture(),
            False,
            self.save_path,
            300.0,
        )
        axes = plt.gcf().axes[0]
        self.assertEqual([line.get_linewidth() for line in axes.lines], [3.5])

    def test_names_the_title_subtitle_and_axis_labels(self) -> None:
        """Test that the figure carries its title, subtitle, and axis labels."""
        _output_plotting.plot_time_history(
            output_plotting_fixtures.make_times_fixture(),
            output_plotting_fixtures.make_three_series_fixture(),
            output_plotting_fixtures.make_three_labels_fixture(),
            output_plotting_fixtures.make_three_colors_fixture(),
            "Example Airplane Forces",
            "(in Wind Axes)",
            "Force (N)",
            output_plotting_fixtures.make_figure_size_fixture(),
            False,
            self.save_path,
            300.0,
        )
        figure = plt.gcf()
        axes = figure.axes[0]
        self.assertEqual(figure.get_suptitle(), "Example Airplane Forces")
        self.assertEqual(axes.get_title(), "(in Wind Axes)")
        self.assertEqual(axes.get_xlabel(), "Time (s)")
        self.assertEqual(axes.get_ylabel(), "Force (N)")

    def test_omits_the_subtitle_when_it_is_empty(self) -> None:
        """Test that a figure given no subtitle draws none."""
        _output_plotting.plot_time_history(
            output_plotting_fixtures.make_times_fixture(),
            output_plotting_fixtures.make_three_series_fixture(),
            output_plotting_fixtures.make_three_labels_fixture(),
            output_plotting_fixtures.make_three_colors_fixture(),
            "Example Airplane Aerodynamic Angles",
            "",
            "Angle (deg)",
            output_plotting_fixtures.make_figure_size_fixture(),
            False,
            self.save_path,
            300.0,
        )
        self.assertEqual(plt.gcf().axes[0].get_title(), "")

    def test_sizes_the_figure_as_it_was_asked_to(self) -> None:
        """Test that the figure takes the width and height it was given."""
        figure_size_in = output_plotting_fixtures.make_figure_size_fixture()
        _output_plotting.plot_time_history(
            output_plotting_fixtures.make_times_fixture(),
            output_plotting_fixtures.make_three_series_fixture(),
            output_plotting_fixtures.make_three_labels_fixture(),
            output_plotting_fixtures.make_three_colors_fixture(),
            "Example Airplane Forces",
            "(in Wind Axes)",
            "Force (N)",
            figure_size_in,
            False,
            self.save_path,
            300.0,
        )
        npt.assert_allclose(plt.gcf().get_size_inches(), figure_size_in)

    def test_hides_the_top_and_right_spines(self) -> None:
        """Test that the plot keeps only the spines its axis labels sit beside."""
        _output_plotting.plot_time_history(
            output_plotting_fixtures.make_times_fixture(),
            output_plotting_fixtures.make_three_series_fixture(),
            output_plotting_fixtures.make_three_labels_fixture(),
            output_plotting_fixtures.make_three_colors_fixture(),
            "Example Airplane Forces",
            "(in Wind Axes)",
            "Force (N)",
            output_plotting_fixtures.make_figure_size_fixture(),
            False,
            self.save_path,
            300.0,
        )
        axes = plt.gcf().axes[0]
        self.assertFalse(axes.spines.top.get_visible())
        self.assertFalse(axes.spines.right.get_visible())
        self.assertTrue(axes.spines.bottom.get_visible())
        self.assertTrue(axes.spines.left.get_visible())

    def test_colors_the_axis_labels_with_the_text_color(self) -> None:
        """Test that the axis labels take the color the visualizations' text takes."""
        _output_plotting.plot_time_history(
            output_plotting_fixtures.make_times_fixture(),
            output_plotting_fixtures.make_three_series_fixture(),
            output_plotting_fixtures.make_three_labels_fixture(),
            output_plotting_fixtures.make_three_colors_fixture(),
            "Example Airplane Forces",
            "(in Wind Axes)",
            "Force (N)",
            output_plotting_fixtures.make_figure_size_fixture(),
            False,
            self.save_path,
            300.0,
        )
        axes = plt.gcf().axes[0]
        text_color = matplotlib.colors.to_rgba(_output_plotting._TEXT_COLOR_NORMALIZED)
        self.assertEqual(
            matplotlib.colors.to_rgba(axes.xaxis.label.get_color()), text_color
        )
        self.assertEqual(
            matplotlib.colors.to_rgba(axes.yaxis.label.get_color()), text_color
        )

    def test_leaves_the_backgrounds_transparent(self) -> None:
        """Test that neither the figure nor the plot draws an opaque background.

        A transparent background is what lets a saved PNG sit on a page of any color.
        """
        _output_plotting.plot_time_history(
            output_plotting_fixtures.make_times_fixture(),
            output_plotting_fixtures.make_three_series_fixture(),
            output_plotting_fixtures.make_three_labels_fixture(),
            output_plotting_fixtures.make_three_colors_fixture(),
            "Example Airplane Forces",
            "(in Wind Axes)",
            "Force (N)",
            output_plotting_fixtures.make_figure_size_fixture(),
            False,
            self.save_path,
            300.0,
        )
        figure = plt.gcf()
        self.assertEqual(matplotlib.colors.to_rgba(figure.get_facecolor())[3], 0.0)
        self.assertEqual(
            matplotlib.colors.to_rgba(figure.axes[0].get_facecolor())[3], 0.0
        )

    def test_legend_lists_every_series_at_one_line_width(self) -> None:
        """Test that the legend names every series and draws its lines one width.

        The lines thin across the plot to keep them all visible, which the legend has no
        reason to repeat.
        """
        labels = output_plotting_fixtures.make_three_labels_fixture()
        _output_plotting.plot_time_history(
            output_plotting_fixtures.make_times_fixture(),
            output_plotting_fixtures.make_three_series_fixture(),
            labels,
            output_plotting_fixtures.make_three_colors_fixture(),
            "Example Airplane Forces",
            "(in Wind Axes)",
            "Force (N)",
            output_plotting_fixtures.make_figure_size_fixture(),
            False,
            self.save_path,
            300.0,
        )
        legend = plt.gcf().axes[0].get_legend()
        self.assertIsNotNone(legend)
        assert legend is not None
        self.assertEqual([text.get_text() for text in legend.get_texts()], labels)
        self.assertEqual(
            [line.get_linewidth() for line in legend.get_lines()], [2.5, 2.5, 2.5]
        )

    def test_does_not_save_when_it_is_not_asked_to(self) -> None:
        """Test that a figure the caller only wants shown is not written to disk."""
        _output_plotting.plot_time_history(
            output_plotting_fixtures.make_times_fixture(),
            output_plotting_fixtures.make_three_series_fixture(),
            output_plotting_fixtures.make_three_labels_fixture(),
            output_plotting_fixtures.make_three_colors_fixture(),
            "Example Airplane Forces",
            "(in Wind Axes)",
            "Force (N)",
            output_plotting_fixtures.make_figure_size_fixture(),
            False,
            self.save_path,
            100.0,
        )
        self.assertFalse(self.save_path.exists())

    def test_saves_at_the_resolution_it_was_given(self) -> None:
        """Test that the saved PNG's size is its size in inches times its resolution."""
        _output_plotting.plot_time_history(
            output_plotting_fixtures.make_times_fixture(),
            output_plotting_fixtures.make_three_series_fixture(),
            output_plotting_fixtures.make_three_labels_fixture(),
            output_plotting_fixtures.make_three_colors_fixture(),
            "Example Airplane Forces",
            "(in Wind Axes)",
            "Force (N)",
            (6.0, 4.0),
            True,
            self.save_path,
            100.0,
        )
        image = matplotlib.image.imread(self.save_path)
        self.assertEqual(image.shape[0], 400)
        self.assertEqual(image.shape[1], 600)
