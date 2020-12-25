try:
    from books.cct.cctpy.cctpy import *
except ModuleNotFoundError:
    pass

from cctpy import *


if True:
    DL2 = 2.1162209
    GAP3 = 0.1978111
    QS3_LEN = 0.2382791
    QS3_GRADIENT = -7.3733
    QS3_SECOND_GRADIENT = -45.31 * 2.0
    QS3_APERTURE = 60 * MM
    big_r_part2 = 0.95
    CCT_APERTURE = 80 * MM

    small_r_gap = 15 * MM
    small_r_innerest = 83 * MM
    agcct_small_r_in = small_r_innerest
    agcct_small_r_out = small_r_innerest + small_r_gap
    dipole_cct_small_r_in = small_r_innerest + small_r_gap * 2
    dipole_cct_small_r_out = small_r_innerest + small_r_gap * 3

    dipole_cct_winding_num: int = 128
    agcct_winding_nums: List[int] = [21, 50, 50]

    dipole_cct_bending_angle = 67.5
    dipole_cct_bending_rad = BaseUtils.angle_to_radian(dipole_cct_bending_angle)
    agcct_bending_angles: List[float] = [8 + 3.716404, 8 + 19.93897, 8 + 19.844626]
    agcct_bending_angles_rad: List[float] = BaseUtils.angle_to_radian(
        agcct_bending_angles
    )

    dipole_cct_tilt_angles = numpy.array([30.0, 80.0, 90.0, 90.0])
    agcct_tilt_angles = numpy.array([90.0, 30.0, 90.0, 90.0])

    dipole_cct_current = 9664.0
    agcct_current = -6000.0

    disperse_number_per_winding: int = 120

    trajectory_part2 = (
        Trajectory.set_start_point(P2(3.703795764767297, 1.5341624380266456))
        .first_line(P2(1, 1), DL2)
        .add_arc_line(0.95, True, dipole_cct_bending_angle)
        .add_strait_line(GAP3 * 2 + QS3_LEN)
    )

    beamline = Beamline(None)
    beamline.magnets.append(
        CCT.create_cct_along(
            trajectory=trajectory_part2,
            s=DL2,
            big_r=big_r_part2,
            small_r=dipole_cct_small_r_in,
            bending_angle=dipole_cct_bending_angle,
            tilt_angles=dipole_cct_tilt_angles,
            winding_number=dipole_cct_winding_num,
            current=dipole_cct_current,
            starting_point_in_ksi_phi_coordinate=P2.origin(),
            end_point_in_ksi_phi_coordinate=P2(
                2 * math.pi * dipole_cct_winding_num, -dipole_cct_bending_rad
            ),
            disperse_number_per_winding=disperse_number_per_winding,
        )
    )

    beamline.magnets.append(
        CCT.create_cct_along(
            trajectory=trajectory_part2,
            s=DL2,
            big_r=big_r_part2,
            small_r=dipole_cct_small_r_out,  # diff
            bending_angle=dipole_cct_bending_angle,
            tilt_angles=-dipole_cct_tilt_angles,  # diff ⭐
            winding_number=dipole_cct_winding_num,
            current=dipole_cct_current,
            starting_point_in_ksi_phi_coordinate=P2.origin(),
            end_point_in_ksi_phi_coordinate=P2(
                -2 * math.pi * dipole_cct_winding_num, -dipole_cct_bending_rad
            ),
            disperse_number_per_winding=disperse_number_per_winding,
        )
    )

    agcct_index = 0
    agcct_start_in = P2.origin()
    agcct_start_out = P2.origin()
    agcct_end_in = P2(
        ((-1.0) ** agcct_index) * 2 * math.pi * agcct_winding_nums[agcct_index],
        -agcct_bending_angles_rad[agcct_index],
    )
    agcct_end_out = P2(
        ((-1.0) ** (agcct_index + 1)) * 2 * math.pi * agcct_winding_nums[agcct_index],
        -agcct_bending_angles_rad[agcct_index],
    )
    beamline.magnets.append(
        CCT.create_cct_along(
            trajectory=trajectory_part2,
            s=DL2,
            big_r=big_r_part2,
            small_r=agcct_small_r_in,
            bending_angle=agcct_bending_angles[agcct_index],
            tilt_angles=-agcct_tilt_angles,
            winding_number=agcct_winding_nums[agcct_index],
            current=agcct_current,
            starting_point_in_ksi_phi_coordinate=agcct_start_in,
            end_point_in_ksi_phi_coordinate=agcct_end_in,
            disperse_number_per_winding=disperse_number_per_winding,
        )
    )

    beamline.magnets.append(
        CCT.create_cct_along(
            trajectory=trajectory_part2,
            s=DL2,
            big_r=big_r_part2,
            small_r=agcct_small_r_out,
            bending_angle=agcct_bending_angles[agcct_index],
            tilt_angles=agcct_tilt_angles,
            winding_number=agcct_winding_nums[agcct_index],
            current=agcct_current,
            starting_point_in_ksi_phi_coordinate=agcct_start_out,
            end_point_in_ksi_phi_coordinate=agcct_end_out,
            disperse_number_per_winding=disperse_number_per_winding,
        )
    )

    for ignore in range(len(agcct_bending_angles) - 1):
        agcct_index += 1
        agcct_start_in = agcct_end_in + P2(
            0,
            -agcct_bending_angles_rad[agcct_index - 1]
            / agcct_winding_nums[agcct_index - 1],
        )
        agcct_start_out = agcct_end_out + P2(
            0,
            -agcct_bending_angles_rad[agcct_index - 1]
            / agcct_winding_nums[agcct_index - 1],
        )
        agcct_end_in = agcct_start_in + P2(
            ((-1) ** agcct_index) * 2 * math.pi * agcct_winding_nums[agcct_index],
            -agcct_bending_angles_rad[agcct_index],
        )
        agcct_end_out = agcct_start_out + P2(
            ((-1) ** (agcct_index + 1)) * 2 * math.pi * agcct_winding_nums[agcct_index],
            -agcct_bending_angles_rad[agcct_index],
        )
        beamline.magnets.append(
            CCT.create_cct_along(
                trajectory=trajectory_part2,
                s=DL2,
                big_r=big_r_part2,
                small_r=agcct_small_r_in,
                bending_angle=agcct_bending_angles[agcct_index],
                tilt_angles=-agcct_tilt_angles,
                winding_number=agcct_winding_nums[agcct_index],
                current=agcct_current,
                starting_point_in_ksi_phi_coordinate=agcct_start_in,
                end_point_in_ksi_phi_coordinate=agcct_end_in,
                disperse_number_per_winding=disperse_number_per_winding,
            )
        )

        beamline.magnets.append(
            CCT.create_cct_along(
                trajectory=trajectory_part2,
                s=DL2,
                big_r=big_r_part2,
                small_r=agcct_small_r_out,
                bending_angle=agcct_bending_angles[agcct_index],
                tilt_angles=agcct_tilt_angles,
                winding_number=agcct_winding_nums[agcct_index],
                current=agcct_current,
                starting_point_in_ksi_phi_coordinate=agcct_start_out,
                end_point_in_ksi_phi_coordinate=agcct_end_out,
                disperse_number_per_winding=disperse_number_per_winding,
            )
        )

    # Plot3.plot_line2(trajectory_part2)
    # Plot3.plot_beamline(beamline, ["r", "y-"] * 3)
    # Plot3.set_center(trajectory_part2.point_at(DL2).to_p3(), 5)
    # Plot3.show()

    print(beamline.magnetic_field_at(trajectory_part2.point_at(DL2).to_p3()))
    print(beamline.magnetic_field_at(trajectory_part2.point_at(DL2+0.5).to_p3()))

if True:
    DL2 = 2.1162209
    GAP3 = 0.1978111
    QS3_LEN = 0.2382791
    QS3_GRADIENT = -7.3733
    QS3_SECOND_GRADIENT = -45.31 * 2.0
    QS3_APERTURE = 60 * MM
    big_r_part2 = 0.95
    CCT_APERTURE = 80 * MM

    small_r_gap = 15 * MM
    small_r_innerest = 83 * MM
    agcct_small_r_in = small_r_innerest
    agcct_small_r_out = small_r_innerest + small_r_gap
    dipole_cct_small_r_in = small_r_innerest + small_r_gap * 2
    dipole_cct_small_r_out = small_r_innerest + small_r_gap * 3

    dipole_cct_winding_num: int = 128
    agcct_winding_nums: List[int] = [21, 50, 50]

    dipole_cct_bending_angle = 67.5
    dipole_cct_bending_rad = BaseUtils.angle_to_radian(dipole_cct_bending_angle)
    agcct_bending_angles: List[float] = [8 + 3.716404, 8 + 19.93897, 8 + 19.844626]

    dipole_cct_tilt_angles = numpy.array([30.0, 80.0, 90.0, 90.0])
    agcct_tilt_angles = numpy.array([90.0, 30.0, 90.0, 90.0])

    dipole_cct_current = 9664.0
    agcct_current = -6000.0

    disperse_number_per_winding: int = 120

    # trajectory_part2 = (
    #     Trajectory.set_start_point(P2(3.703795764767297, 1.5341624380266456))
    #     .first_line(P2(1, 1), DL2)
    #     .add_arc_line(0.95, True, dipole_cct_bending_angle)
    #     .add_strait_line(GAP3*2+QS3_LEN)
    # )

    #############################################
    beamline = (
        Beamline.set_start_point(P2(3.703795764767297, 1.5341624380266456))
        .first_drift(P2(1, 1), DL2)
        .append_agcct(
            big_r=0.95,
            small_rs=[
                dipole_cct_small_r_out,
                dipole_cct_small_r_in,
                agcct_small_r_out,
                agcct_small_r_in,
            ],
            bending_angles=BaseUtils.list_multiply(agcct_bending_angles, -1),
            tilt_angles=[dipole_cct_tilt_angles, agcct_tilt_angles],
            winding_numbers=[[dipole_cct_winding_num], agcct_winding_nums],
            currents=[dipole_cct_current, agcct_current],
            disperse_number_per_winding=120,
        )
        .append_drift(GAP3 * 2 + QS3_LEN)
    )

    # Plot3.plot_line2(beamline.trajectory)
    # Plot3.plot_beamline(beamline, ["r", "y-"] * 3)
    # Plot3.set_center(beamline.trajectory.point_at(DL2).to_p3(), 5)
    # Plot3.show()

    print(beamline.magnetic_field_at(beamline.trajectory.point_at(DL2).to_p3()))
    print(beamline.magnetic_field_at(beamline.trajectory.point_at(DL2+0.5).to_p3()))
