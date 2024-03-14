import numpy as np
from fractions import Fraction


def Schwinger(j: int = 2):
    """
    Only works for integer j, too lazy to consider half-integer j's
    """
    print("\\begin{align*}")
    m_values = np.arange(j, -j-1, -1)
    for m in m_values:
        for mp in m_values:
            k_max = np.min([j+m, j-mp])
            k_min = np.max([0, m-mp])
            k_values = np.arange(k_max, k_min-1, -1)
            for k in k_values:
                if k == k_values[0]:
                    print(f"    \\mathcal{{D}}^{{(2)}}_{{{mp}, {m}}}", end="=&")
                sign_power = k-m+mp

                if (-1)**sign_power < 0:
                    print("-", end="")
                else:
                    if k != k_values[0]:
                        print("+", end="")

                coeff = Fraction(
                    np.math.factorial(j + m) * np.math.factorial(j - m)
                    * np.math.factorial(j + mp) * np.math.factorial(j - mp)
                    /
                    (np.math.factorial(j + m - k) * np.math.factorial(k)
                        * np.math.factorial(j - k - mp) * np.math.factorial(k-m+mp))**2
                )

                if np.sqrt(coeff.numerator/coeff.denominator).is_integer() and coeff.numerator/coeff.denominator != 1:
                    print(int(np.sqrt(coeff.numerator/coeff.denominator)), end="")
                elif coeff.numerator/coeff.denominator != 1:
                    print(f"\\sqrt{{{int(coeff.numerator / coeff.denominator)}}}", end="")
                if not ((m == 0) and (mp == 0)):
                    print("e^{{-i(", end="")

                    if mp == 1:
                        print(f"\\alpha", end="")
                    elif mp == -1:
                        print(f"-\\alpha", end="")
                    elif mp != 0:
                        print(f"{mp}\\alpha", end="")

                    if m == 1:
                        if mp != 0:
                            print("+", end="")
                        print(f"\\gamma", end="")
                    elif m == -1:
                        print(f"-\\gamma", end="")
                    elif m != 0:
                        if mp != 0 and m > 1:
                            print("+", end="")
                        print(f"{m}\\gamma", end="")

                    print(")}}", end="")
                cos_power = 2*j-2*k+m-mp
                if cos_power == 1:
                    print(f"\\cos(\\frac\\beta 2)", end="")
                elif cos_power != 0:
                    print(f"\\cos[{cos_power}](\\frac{{\\beta}}{{2}})", end="")

                sin_power = 2*k-m+mp
                if sin_power == 1:
                    print(f"\\sin(\\frac\\beta 2)", end="")
                elif sin_power != 0:
                    print(f"\\sin[{sin_power}](\\frac{{\\beta}}{{2}})", end="")

                if k == k_values[-1] and not (m == m_values[-1] and mp == m_values[-1]):
                    print("\\\\")
                else:
                    print("", end="")
    print("\n\\end{align*}")


def J_minus(m1, m2, j=2, j1=1, j2=1):
    m = m1+m2+1
    LHS = (j+m)*(j-m+1)
    RHS1 = (j1+m1+1)*(j1-m1)
    RHS2 = (j2+m2+1)*(j2-m2)

    def simplify_root(X):
        if X != 0:
            if np.sqrt(X).is_integer():
                print(f"{int(np.sqrt(X))}", end="")
            else:
                print(f"\\sqrt{{{X}}}", end="")
    simplify_root(LHS)
    print(f"\\braket{{{m1},{m2}|{j},{m-1}}}", end="=&")
    simplify_root(RHS1)
    if RHS1 != 0:
        print(f"\\braket{{{m1+1},{m2}|{j},{m}}}", end="")
    if RHS1 != 0 and RHS2 != 0:
        print("+", end="")
    simplify_root(RHS2)
    if RHS2 != 0:
        print(f"\\braket{{{m1},{m2+1}|{j},{m}}}", end="")
    print(f"\\\\")


def number_1():
    Schwinger()


def number_4():
    m1m2_values = (
        *[(i, 1) for i in (0, -1)],
        *[(i, 0) for i in (1, 0, -1)],
        *[(i, -1) for i in (1, 0, -1)]
    )
    for m1m2 in m1m2_values:
        J_minus(*m1m2)


if __name__ == '__main__':
    number_1()
