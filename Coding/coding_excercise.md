# Documentation for Numerical Polynomial Algebra Course

Aim: Download sympy from github, create a branch and modify and test the code.

1. Clone sympy from github

git clone git://github.com/sympy/sympy.git

2. install mpmath

conda install mpmath

3. cd sympy

run python, then import sympy ...

4. create your experimental branch

git branch mybranch
git checkout mybranch

5. install numpy

conda install numpy


6. run tests for groebnertools

 ... start python, then

>>> import sympy
>>> sympy.test("sympy/polys/tests/test_groebnertools.py")

check out the testing code

7. change some code

* in sympy/polys/groebnertools.py change spoly, then run the tests again ...

def spoly(p1, p2, ring):
    """
    Compute LCM(LM(p1), LM(p2))/LM(p1)*p1 - LCM(LM(p1), LM(p2))/LM(p2)*p2
    This is the S-poly provided p1 and p2 are monic
    """

    print("*** this is my own spoly code ***")

    LM1 = p1.LM
    LM2 = p2.LM
    LCM12 = ring.monomial_lcm(LM1, LM2)
    m1 = ring.monomial_div(LCM12, LM1)
    m2 = ring.monomial_div(LCM12, LM2)
    s1 = p1.mul_monom(m1)
    s2 = p2.mul_monom(m2)
    s = s1 - s2
    return s
