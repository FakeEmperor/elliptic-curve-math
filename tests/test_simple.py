import unittest
from kily.math.ec import *


class TestSmoke(unittest.TestCase):

    def test_check_ok(self):
        ec1 = EllipticCurveModulo(2, 3, 97)
        self.assertTrue(ec1.check())

        ec2 = EllipticCurveModulo(4, 5, 97)
        self.assertTrue(ec2.check())

        ec3 = EllipticCurveModulo(1, 6, 97)
        self.assertTrue(ec3.check())

    def test_check_fail(self, raise_error=False):
        ec_wrong_modulo = EllipticCurveModulo(2, 3, 91)
        try:
            self.assertFalse(ec_wrong_modulo.check(raise_error=raise_error))
        except ValueError as e:
            if not raise_error:
                raise e

        ec_wrong_crypto_params = EllipticCurveModulo(1, 5, 97)
        try:
            self.assertFalse(ec_wrong_crypto_params.check(raise_error=raise_error))
        except ValueError as e:
            if not raise_error:
                raise e

        ec_wrong_modulo_2 = EllipticCurveModulo(2, 3, 2)
        try:
            self.assertFalse(ec_wrong_modulo_2.check(raise_error=raise_error))
        except ValueError as e:
            if not raise_error:
                raise e

    def test_check_raise(self):
        self.test_check_fail(True)

    def test_point_ok(self):
        ec1 = EllipticCurveModulo(2, 3, 97)
        print("\nField size: {}".format(ec1.size))
        for p in ec1.point_begin():
            print("Point: {}".format(p))
            self.assertTrue(ec1.valid(p))

    def test_point_gen_unique(self):
        ec1 = EllipticCurveModulo(2, 3, 97)
        points = set()
        flat_points = []
        print("\nField size: {}".format(ec1.size))
        s = 0
        for p in ec1.point_begin():
            print("Point: {}".format(p))
            points.add(p)
            s += 1
            if p in flat_points:
                self.fail("Generated point {} is not unique".format(points))
        self.assertTrue(s, len(points))

    def test_hasse_check(self):
        ec1 = EllipticCurveModulo(2, 3, 97)
        self.assertTrue(ec1.size_hasse_check())

    def test_point_fail(self):
        ec1 = EllipticCurveModulo(2, 3, 4)
        self.assertFalse(ec1.valid(ec1.point(2, 0)))

    def test_repr(self):
        # sanity check
        ec1 = EllipticCurveModulo(2, 3, 97)
        self.assertEquals(repr(ec1), "EllipticCurveModulo(2, 3, 97) = 'y^2 = x^3 + 2x + 3 mod 97'")


class TestSimple(unittest.TestCase):

    def test_two_point_add(self):
        pass