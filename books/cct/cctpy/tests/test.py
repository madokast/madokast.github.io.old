import unittest
import random

import numpy as np

from cctpy import (
    GPU_ON,
    BaseUtils,
    P2,
    P3,
    StraightLine2,
    Trajectory,
    Plot3,
    CCT,
    LocalCoordinateSystem,
    MM,
    M,
    QS,
)


class AssertCase(unittest.TestCase):
    def test_equal(self):
        BaseUtils.equal(1, 1, msg="不相等")
        with self.assertRaises(AssertionError):
            BaseUtils.equal(1, 2, msg="不相等")

    def test_equal2(self):
        BaseUtils.equal(P2(), P2(1e-10), msg="不相等")
        with self.assertRaises(AssertionError):
            BaseUtils.equal(P2(), P2(1e-5), msg="不相等")


class P2Case(unittest.TestCase):
    def test_p2_len(self):
        p = P2(1, 1)
        self.assertTrue(BaseUtils.equal(p.length(), np.sqrt(2)))

    def test_p2_norm(self):
        p = P2(1, 1)
        self.assertTrue(
            BaseUtils.equal(p.normalize(), P2(np.sqrt(2) / 2, np.sqrt(2) / 2))
        )

    def test_p2_change_len(self):
        for i in range(1, 100):
            a = random.random() + i
            b = random.random() + i
            p = P2(a, b)
            self.assertTrue(p.change_length(1).normalize() == P2(a, b).normalize())

    def test_p2_copy(self):
        for i in range(100):
            p = P2(random.random() + i, random.random() + i)
            self.assertTrue(BaseUtils.equal(p, p.copy()))
            self.assertFalse(p is p.copy())

    def test_p2_add(self):
        for i in range(100):
            a = random.random() + i
            b = random.random() + i
            c = random.random() + i
            d = random.random() + i
            self.assertEqual(P2(a, b) + P2(c, d), P2(a + c, b + d))

    def test_p2_neg(self):
        for i in range(100):
            a = random.random() + i
            b = random.random() + i
            self.assertEqual(P2(-a, -b), -P2(a, b))

    def test_p2_sub(self):
        for i in range(100):
            a = random.random() + i
            b = random.random() + i
            c = random.random() + i
            d = random.random() + i
            self.assertEqual(P2(a, b) - P2(c, d), P2(a - c, b - d))

    def test_p2_iadd(self):
        for i in range(100):
            a = random.random() + i
            b = random.random() + i
            c = random.random() + i
            d = random.random() + i
            p1 = P2(a, b)
            p2 = P2(c, d)
            p1_copy = p1.copy()
            p1 += p2
            self.assertEqual(p1, p1_copy + p2)

    def test_p2_isub(self):
        for i in range(100):
            a = random.random() + i
            b = random.random() + i
            c = random.random() + i
            d = random.random() + i
            p1 = P2(a, b)
            p2 = P2(c, d)
            p1_copy = p1.copy()
            p1 -= p2
            self.assertEqual(p1, p1_copy - p2)

    def test_p2_rotate(self):
        for phi in np.linspace(0.1, 2 * np.pi * 0.99, 100):
            p = P2(1.0)
            p = p.rotate(phi)
            self.assertTrue(BaseUtils.equal(p.angle_to_x_axis(), phi))

    def test_angle_to_x_axis(self):
        for i in range(100):
            phi = random.random() * np.pi * 1.9 + 1e-6
            p = P2(1.0)
            p = p.rotate(phi)
            self.assertTrue(BaseUtils.equal(p.angle_to_x_axis(), phi))

    def test_mul1(self):
        for i in range(100):
            a = random.random() + i + 1
            b = random.random() + i + 1
            c = random.random() + i + 1
            self.assertTrue(P2(a, b) * c, P2(a * c, b * c))

    def test_mul2(self):
        for i in range(100):
            a = random.random() + i + 1
            b = random.random() + i + 1
            c = random.random() + i + 1
            d = random.random() + i + 1
            self.assertTrue(
                BaseUtils.equal(P2(a, b) * P2(c, d), a * c + b * d, msg="error")
            )

    def test_angle_to(self):
        for i in range(100):
            a0 = random.random() * np.pi * 1.9 + 1e-6
            b0 = random.random() * np.pi * 1.9 + 1e-6
            a = max(a0, b0)
            b = min(a0, b0)
            p1 = P2(1.0)
            p2 = P2(1.0)

            p1 = p1.rotate(a)
            p2 = p2.rotate(b)

            diff = a - b

            # print(p1.angle_to(p2),diff)

            self.assertTrue(
                BaseUtils.equal(
                    p2.angle_to(p1), diff #if diff >0 else 2 * np.pi - diff
                )
            )

    def test_to_p3(self):
        for i in range(100):
            a = random.random() + i
            b = random.random() + i
            self.assertTrue(P2(a, b).to_p3() == P3(a, b, 0.0))


class P3Case(unittest.TestCase):
    def test_length(self):
        for i in range(100):
            a = random.random() + i
            b = random.random() + i
            c = random.random() + i
            self.assertTrue(
                BaseUtils.equal(
                    P3(a, b, c).length(), float(np.sqrt(a ** 2 + b ** 2 + c ** 2))
                )
            )

    def test_normalize(self):
        for i in range(100):
            a = random.random() + i
            b = random.random() + i
            c = random.random() + i
            self.assertTrue(BaseUtils.equal(P3(a, b, c).normalize().length(), 1.0))

    def test_change_length(self):
        for i in range(100):
            a = random.random() + i
            b = random.random() + i
            c = random.random() + i
            d = random.random() + i
            self.assertTrue(BaseUtils.equal(P3(a, b, c).change_length(d).length(), d))

    def test_copy(self):
        for i in range(100):
            a = random.random() + i
            b = random.random() + i
            c = random.random() + i
            self.assertTrue(P3(a, b, c) == P3(a, b, c).copy())

    def test_p3_add(self):
        for i in range(100):
            a = random.random() + i
            b = random.random() + i
            c = random.random() + i
            d = random.random() + i
            e = random.random() + i
            f = random.random() + i
            self.assertEqual(P3(a, b, e) + P3(c, d, f), P3(a + c, b + d, e + f))

    def test_p3_neg(self):
        for i in range(100):
            a = random.random() + i
            b = random.random() + i
            c = random.random() + i
            self.assertEqual(P3(-a, -b, -c), -P3(a, b, c))

    def test_p3_sub(self):
        for i in range(100):
            a = random.random() + i
            b = random.random() + i
            c = random.random() + i
            d = random.random() + i
            e = random.random() + i
            f = random.random() + i
            self.assertEqual(P3(a, b, e) - P3(c, d, f), P3(a - c, b - d, e - f))

    def test_p2_iadd(self):
        for i in range(100):
            a = random.random() + i
            b = random.random() + i
            c = random.random() + i
            d = random.random() + i
            e = random.random() + i
            f = random.random() + i
            p1 = P3(a, b, e)
            p2 = P3(c, d, f)
            p1_copy = p1.copy()
            p1 += p2
            self.assertEqual(p1, p1_copy + p2)

    def test_p3_isub(self):
        for i in range(100):
            a = random.random() + i
            b = random.random() + i
            c = random.random() + i
            d = random.random() + i
            e = random.random() + i
            f = random.random() + i
            p1 = P3(a, b, e)
            p2 = P3(c, d, f)
            p1_copy = p1.copy()
            p1 -= p2
            self.assertEqual(p1, p1_copy - p2)

    def test_mul1(self):
        for i in range(100):
            a = random.random() + i + 1
            b = random.random() + i + 1
            c = random.random() + i + 1
            d = random.random() + i + 1
            self.assertTrue(P3(a, b, d) * c, P3(a * c, b * c, d * c))

    def test_mul2(self):
        for i in range(100):
            a = random.random() + i + 1
            b = random.random() + i + 1
            c = random.random() + i + 1
            d = random.random() + i + 1
            e = random.random() + i + 1
            f = random.random() + i + 1
            self.assertTrue(
                BaseUtils.equal(
                    P3(a, b, e) * P3(c, d, f), a * c + b * d + e * f, msg="error"
                )
            )

    def test__matmul__(self):
        x = P3.x_direct()
        y = P3.y_direct()
        z = P3.z_direct()

        self.assertEqual(x @ y, z)
        self.assertEqual(y @ z, x)
        self.assertEqual(z @ x, y)

        self.assertEqual(x @ -y, -z)
        self.assertEqual(y @ -z, -x)
        self.assertEqual(z @ -x, -y)


class TrajectoryTestCase(unittest.TestCase):
    def test_straight(self):
        s = StraightLine2(1, np.array([1.0, 1.0]), np.array([0.0, 0.0]))
        self.assertTrue(BaseUtils.equal(s.get_length(), 1))

    def test_tr_01(self):
        t = (
            Trajectory(StraightLine2(2.0, P2(1, 0), P2(0, 0)))
            .add_arc_line(0.95, False, 22.5)
            .add_strait_line(1.5)
            .add_arc_line(0.95, False, 22.5)
            .add_strait_line(2.0 + 2.2)
            .add_arc_line(0.95, True, 67.5)
            .add_strait_line(1.5)
            .add_arc_line(0.95, True, 67.5)
            .add_strait_line(2.2)
        )

        self.assertTrue(BaseUtils.equal(t.direct_at_end(), P2(0, -1)))

        # Plot3.plot3d(t.line_and_color())
        # Plot3.show()

    def test_tr_02(self):
        """
        彩蛋，把绘图代码注释取消即可
        Returns
        -------

        """
        c1 = (
            Trajectory(StraightLine2(0.01, P2(0, 1), P2(0, 0)))
            .add_arc_line(1, True, 135)
            .add_arc_line(0.01, True, 90)
            .add_strait_line(0.1)
            .add_arc_line(0.01, True, 90)
            .add_arc_line(0.9, False, 360 - 90)
            .add_arc_line(0.01, True, 90)
            .add_strait_line(0.1)
            .add_arc_line(0.01, True, 90)
            .add_arc_line(1, True, 135)
        )

        c2 = c1 + P2(3, 0)

        t = Trajectory(StraightLine2(0.8, P2(1, 0), P2(6, 1))).add_arc_line(
            0.01, True, 90
        ).add_strait_line(0.2).add_arc_line(0.01, True, 90).add_strait_line(
            0.7
        ).add_arc_line(
            0.01, False, 90
        ).add_strait_line(
            1.7
        ).add_arc_line(
            0.01, True, 90
        ).add_strait_line(
            0.2
        ).add_arc_line(
            0.01, True, 90
        ).add_strait_line(
            1.7
        ).add_arc_line(
            0.01, False, 90
        ).add_strait_line(
            0.7
        ).add_arc_line(
            0.01, True, 90
        ).add_strait_line(
            0.2
        ).add_arc_line(
            0.01, True, 90
        ).add_strait_line(
            0.8
        ) + P2(
            0.5, 0
        )

        Plot3.plot_line2(c1, describe='r')
        Plot3.plot_line2(c2, describe='b')
        Plot3.plot_line2(t, describe='g')
        Plot3.set_center()
        
        Plot3.show()

        self.assertTrue(True)


class CCTTestCase(unittest.TestCase):
    def test_magnet_at(self):
        cct = CCT(
            LocalCoordinateSystem.global_coordinate_system(),
            0.95,
            83 * MM + 15 * MM * 2,
            67.5,
            [30.0, 80.0, 90.0, 90.0],
            128,
            -9664,
            P2(0, 0),
            P2(128 * np.pi * 2, 67.5 / 180.0 * np.pi),
        )

        m = cct.magnetic_field_at(P3.origin())
        self.assertEqual(
            m, P3(0.0031436355039083964, -0.00470478301086915, 0.00888627084434009)
        )


class QsTest(unittest.TestCase):
    def plot_qs(self):
        length = 0.2 * M
        aper = 30 * MM
        g = 10.0
        L = 0
        lc = LocalCoordinateSystem()
        qs = QS(lc, length, g, L, aper)

        Plot3.plot_qs(qs)

    def test_quad_0(self):
        """
        测试 qs 四极场
        Returns
        -------

        """
        length = 0.2 * M
        aper = 30 * MM
        g = 10.0
        L = 0
        lc = LocalCoordinateSystem(P3(),-P3.x_direct(),P3.y_direct())
        qs = QS(lc, length, g, L, aper)

        m = qs.magnetic_field_at(P3(10 * MM, 0.1, 0.0))
        self.assertTrue(m == P3(0.0, 0.0, -0.1))

        m = qs.magnetic_field_at(P3(15 * MM, 0.1, 0.0))
        self.assertTrue(m == P3(0.0, 0.0, -0.15))

        m = qs.magnetic_field_at(P3(15 * MM, 0.1, 5 * MM))
        self.assertTrue(m == P3(-0.05, -3.061616997868383e-18, -0.15))

    def test_quad_1(self):
        """
        测试 qs 四极场
        Returns
        -------

        """
        length = 0.2 * M
        aper = 30 * MM
        g = -45.7
        L = 0
        lc = LocalCoordinateSystem(P3(),-P3.x_direct(),P3.y_direct())
        qs = QS(lc, length, g, L, aper)

        m = qs.magnetic_field_at(P3(10 * MM, 0.1, 0))
        self.assertTrue(m == P3(0.0, 0.0, 0.457))

        m = qs.magnetic_field_at(P3(15 * MM, 0.1, 0))
        self.assertTrue(m == P3(0.0, 0.0, 0.6855))

        m = qs.magnetic_field_at(P3(15 * MM, 0.1, 5 * MM))
        self.assertTrue(m == P3(0.2285, 1.399158968025851e-17, 0.6855))


if __name__ == "__main__":
    unittest.main()
