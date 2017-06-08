import math


def is_prime(n):
    """
    from an answer
    # https://stackoverflow.com/questions/18833759/python-prime-number-checker
    # bruteforce check
    """
    if n <= 2:
        return True
    if n % 2 == 0:
        return False
    return all(n % i for i in range(3, int(math.sqrt(n)) + 1, 2))

# from an answer
# https://stackoverflow.com/questions/31074172/elliptic-curve-point-addition-over-a-finite-field-in-python


class ECMPoint(object):
    O = POINT_INFINITY = "inf"

    def __init__(self, x: int, y: int, ec: 'EllipticCurveModulo'):
        self.x = x
        self.y = y
        self.ec = ec

    def __mul__(self, mult: int):
        # dumb mult
        result = self.ec.point(self.x, self.y)
        for i in range(mult):
            result = self.ec.ec_add(result, self)
        return result

    def __add__(self, p: 'ECMPoint'):
        return self.ec.ec_add(self, p)

    def __sub__(self, p: 'ECMPoint'):
        return self.ec.ec_add(self, self.ec.ec_inv(p))

    def __invert__(self):
        return self.ec.ec_inv(self)

    def __repr__(self):
        return "ECMPoint({}, {})".format(self.x, self.y)

    def __eq__(self, p: 'ECMPoint'):
        return p is not None and p != ECMPoint.O and self.x == p.x and self.y == p.y and self.ec == p.ec

    def __hash__(self):
        return hash((self.x, self.y))


# elliptic curve like y^2 = x^3 + a*x + b
class EllipticCurveModulo(object):
    class PointIter(object):
        def __init__(self, gen_point: ECMPoint, ec: 'EllipticCurveModulo'):
            self.gen_point = gen_point
            self.inter_point = None  # type: ECMPoint
            self.ec = ec
            self.max_gen = 1

        def next(self):
            if self.inter_point != ECMPoint.O:
                if self.inter_point is not None:
                    self.inter_point = self.gen_point + self.inter_point
                else:
                    self.inter_point = self.ec.point(self.gen_point.x, self.gen_point.y)
                return self.inter_point if self.inter_point == ECMPoint.O else self.ec.point(self.inter_point.x, self.inter_point.y)
            else:
                raise StopIteration()

        def __iter__(self):
            return self

        # Python 3 compatibility
        def __next__(self):
            return self.next()

    def __init__(self, a: int, b: int, mod: int):
        self._mod = mod
        self._a = a
        self._b = b
        self.__size__ = None  # type: int

    def point_begin(self) -> 'PointIter':
        for x in range(1, self._mod):
            try:
                p = self.eval(x)
                return EllipticCurveModulo.PointIter(p, self)
            except ValueError:
                pass

    def eval_x(self, x: int) -> int:
        x = x % self._mod
        return (x ** 3 + self._a * x + self._b) % self._mod

    def eval(self, x) -> ECMPoint:
        y_square = self.eval_x(x)
        # find a square root by brute force
        for y in range(self._mod):
            if y**2 % self._mod == y_square:
                return self.point(x, y)
        raise ValueError("given x cannot evaluate a point")

    def check(self, crypto_check: bool = True, raise_error=True) -> bool:
        if self._mod <= 2:
            if raise_error:
                raise ValueError("_mod value '{}' is incorrect (<=2)!".format(self._mod))
            return False
        if not is_prime(self._mod):
            if raise_error:
                raise ValueError("_mod value '{}' is incorrect (not prime)!".format(self._mod))
            return False
        if crypto_check and ((4*self._a**3 % self._mod) + (27*self._b**2 % self._mod) % self._mod) % self._mod == 0:
            if raise_error:
                raise ValueError("EC is singular: 4*a^3 + 27*b^2 is zero mod p (crypto check) "
                                 "with a={},b={},p={}".format(self._a, self._b, self._mod))
            return False
        return True

    def valid(self, p: ECMPoint):
        """
        Determine whether we have a valid representation of a point
        on our curve.  We assume that the x and y coordinates
        are always reduced modulo p, so that we can compare
        two points for equality with a simple ==.
        """
        if p == ECMPoint.O:
            return True
        elif p.ec != self:
            return False
        else:
            return ((p.y ** 2 - (p.x ** 3 + self._a * p.x + self._b)) % self._mod == 0 and
                    0 <= p.x < self._mod and 0 <= p.y < self._mod)

    def point(self, x, y) -> ECMPoint:
        return ECMPoint(x, y, self)

    @property
    def size(self) -> int:
        if self.__size__ is None:
            self.__size__ = sum(1 for _ in self.point_begin())  # type: int
        return self.__size__

    def size_hasse_check(self) -> bool:
        s = self.size
        return abs(s - self._mod - 1) < 2*math.sqrt(self._mod)

    def inv_mod(self, x):
        """
        Compute an inverse for x modulo p, assuming that x
        is not divisible by p.
        """
        if x % self._mod == 0:
            raise ZeroDivisionError("Impossible inverse")
        # three-args pow - ( x^p-2 mod p )
        # if mod is prime - result is always an inverse (because x^(p-1) = 1 mod p)
        return pow(x, self._mod-2, self._mod)

    def ec_inv(self, p: ECMPoint, check=True) -> ECMPoint:
        """
        Inverse of the point P on the elliptic curve y^2 = x^3 + ax + b.
        """
        if check and not self.valid(p):
            raise ValueError("Point '{}' is invalid".format(p))
        if p == ECMPoint.O:
            return p
        return ECMPoint(p.x, (-p.y) % self._mod, self)

    def ec_add(self, p: ECMPoint, q: ECMPoint, check=True):
        """
        Sum of the points P and Q on the elliptic curve y^2 = x^3 + ax + b.
        """
        if check and not self.valid(p):
            raise ValueError("Point '{}' is invalid".format(p))

        # Deal with the special cases where either P, Q, or P + Q is
        # the origin.
        if p == ECMPoint.O:
            result = q
        elif q == ECMPoint.O:
            result = p
        elif q == self.ec_inv(p, check=False):
            result = ECMPoint.O
        else:
            # Cases not involving the origin.
            if p == q:
                dydx = (3 * p.x ** 2 + self._a) * self.inv_mod(2 * p.y)
            else:
                dydx = (q.y - p.y) * self.inv_mod(q.x - p.x)

            x = int((dydx ** 2 - p.x - q.x) % self._mod)
            y = int((dydx * (p.x - x) - p.y) % self._mod)
            result = ECMPoint(x, y, self)

        # The above computations *should* have given us another point
        # on the curve.
        assert self.valid(result)
        return result

    def __repr__(self):
        return "EllipticCurveModulo({}, {}, {}) =" \
               " 'y^2 = x^3 + {}x + {} mod {}'".format(*((self._a, self._b, self._mod)*2))

