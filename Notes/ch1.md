---
title: Special topics course Numerical Polynomial Algebra
author: Anand Deopurkar and Markus Hegland
date: ANU winter semester 2018
header-includes:
  $$\newcommand{\N}{\mathbb{N}}$$
  $$\newcommand{\Q}{\mathbb{Q}}$$
  $$\newcommand{\R}{\mathbb{R}}$$
  $$\newcommand{\C}{\mathbb{C}}$$
  $$\newcommand{\T}{\mathcal{T}}$$
  $$\newcommand{\P}{\mathcal{P}}$$
  $$\newcommand{\range}{\operatorname{range}}$$
  $$\newcommand{\null}{\operatorname{null}}$$
  $$\newcommand{\x}{\underline{x}}$$
...
# Chapter 1: Numerical introduction

# introductory section -- numerical aspects

* here we consider polynomial systems of equations for the real and complex fields
* polynomial systems of equations are a generalisation of linear systems of equations which include powers $x_i^k$ and products of powers of the unknowns
* applications of polynomial systems of equations:
    * computing Gaussian quadrature points
    * computing coefficients of Runge-Kutta methods
    * eigenvalue problems
    * rank-k matrix approximation
    * tensor approximations and decompositions
    * polynomial optimisation
    * semidefinite programming
    * machine learning
    * robotics
    * cryptography
    * ...
* for some cases there are stable and well-studied methods using specialised algorithms
* for some cases very large systems can be solved using general (Buchberger-based) algorithms
* in some cases even small system are practically unsolvable
* algorithms in the rational field and finite fields are often effective
* algorithms for the real and complex fields based on floating point numbers are often unstable
* simple example is Euclid's algorithm which is related to Gaussian elimination (without pivoting)
    * problems to solve are Toeplitz systems
    * some numerical problems can be illustrated there
    * algorithms based on stable (QR) algorithms have been investigated

## 1.1 Linear spaces of polynomials

This section is based on section 1.1 in Stetter's book. Here we define the linear space of multivariate polynomials. We first define the monomial basis. A univariate monomial is a power $x^j$ for $x\in\R$ and $j\in\N_0$ where $\N_0$ is the set of non-negative integers. A multivariate monomial is a product of the powers of the variables, i.e.,
$$x^j = x_1^{j_1}\cdots x_s^{j_s}$$ where $x\in\R^s$ and $j\in\N_0^s$. The *set of (s-variate) monomials* $\T^s$ is
$$\T^s = \{x^j \mid x\in\R^s,\; j\in\N_0^s\}.$$
The set of monomials is *closed under multiplication*, i.e.,
if $x^j, x^k\in \T^s$ then $x^j*x^k = x^{j+k}\in \T^s$. As the multiplication is associative we have

**Proposition:** The set $(\T^s,*)$ of monomials is a symmetric semigroup with unit $1$.

*proof:* use map $x^j \rightarrow j$

**Definition:**
A *semigroup* $(S,*)$ is a set $S$ with a binary operation $*$ satisfying the associative law
$$a*(b*c) = (a*b)*c$$

**Examples:**

* non-negative integers $(\N_0,+)$ $$\N_0 = \{0,1,2,\ldots \}$$

    * is symmetric and contains unit 0
    * does not contain negative integers and thus is not a group

* vectors of non-negative integers $(\N_0^s,+)$

**Proposition:** The semigroup $(\T^s,*)$ is isomorph to $(\N_0^s,+)$.

**Definition:** The set of $s$-variate polynomials is $$\P^s = \{p(x) = \sum_{j\in J} c_j x^j \mid J\subset \N_0^s, |J| < \infty, c_j\in\R\}.$$

Here $|J|$ denotes the cardinality (size) of the set $J$. A shorthand notation for this is $$\P^s = \R(\T^s).$$

From the definition it follows that $\P^s$ is a *commutative algebra*, i.e.,

* a real vector space with component-wise addition and multiplication by reals
* a ring with the multiplication defined by $$p(x)q(x) = \sum_{j,k}c_j d_k x^{j+k}$$
where $p(x) = \sum_j c_j x^j$ and $q(x)= \sum_k d_k x^k$.

In practice, the polynomials occurring are sparse, i.e., the size $|J|$ in the sums are small.

-------------------------------------------------

In computations, the degree of the polynomials will be limited. One uses, for example

**Definition:** (monomials with finite degree)

$$\T_d^s = \{x^j \in \T^s \mid |j| < d\}$$

where degree of $x^j$ is $|j| = j_1+\cdots+j_s$

* $\T_d^s$ is not a semigroup

Generating a vector space from this set gives polynomials of degree $d$
$$\P_d^s = \R(\T_d^s).$$

The size $\T_d^s = \binom{d+s}{s}$ is equal to the dimension of the vector space $\P_d^s$.

**Examples:**

* $\T_d^1 = \{1,x,x^2,\ldots,x^d\}$ which generates
$$\P_d^1 = \R(\T_d^1),$$ the univariate polynomials of degree up to $d$ and
* $\T_1^s = \{1,x_1,x_2,\ldots,x_s\}$ which generates $$\P_1^s = \R(\T_1^s),$$ the multivariate polynomials of degree up to one.

In some applications one uses the tensor product of univariate polynomials of a fixed degree and has
$$\T_d^{\otimes s} = \{x^j \in \T^s \mid j_i = 0,\ldots,d, \; i=1,2,\ldots,s\}$$ which generates
$$\P_d^{\otimes s} = \R(\T_d^{\otimes s}) = \P_d^1\otimes \cdots \otimes \P_d^1.$$

---------------------------------------
For a given finite subset of $\T^s$ one defines the vector of basis functions for the corresponding polynomial vector space as array $\x = (x^j)_{j\in J}$. A polynomial in the space $\P_J$ generated by $x^j$ for $j\in J$ then takes the form $$p(x) = c^T \x$$ where $c\in\R^{|J|}$.

For example, if $J=\{0,\ldots,d\}$ then $$\x = (x^d,x^{d-1},\ldots,1)^T$$ (we will order highest degree first). Then $c^T\x$ generates all the elements of $\P^1_d$, i.e. $$\P^1_d = \{c^T\x \mid c\in \R^{d+1}\}.$$ If $J=\{e_1,e_2,\ldots,e_s,0\}$ where $e_i$ is the i-th standard basis vector and 0 the zero vector one has $$\x=(x,1)^T=(x_1,x_2,\ldots,x_s,1)^T$$ which generates the elements of $\P^s_1$, i.e., $$\P_1^s = \{c^T \x \mid c\in \R^{s+1}\}$$ in this case.


**Definition:**
A polynomial system $P(x)$ is an array of polynomials $p_i(x)\in \P^s$ such that $$P(x) = \begin{bmatrix}p_1(x)\\ \vdots \\ p_n(x)\end{bmatrix}.$$

A polynomial system (or system of polynomials) can be expressed as the matrix vector product $$P(x) = C \x$$
for some matrix $C\in \R^{n\times |J|}$.

## 1.2 Solutions of polynomial systems of equations

This section includes basic material plus topics from Section 1.2 of Stetter's book.

**Definition:** A real (complex) solution or zero of a polynomial system is a vector $x\in\R^s$ ($x\in\C^s$) which satisfies $$P(x) = 0.$$

From a mathematical perspective, we are interested in two questions:

* Does a particular polynomial system have a solution?
* What is the nature of the set of all solutions of a polynomial system? (Eg, is it a smooth curve or manifold, a point or a line or vector space?)

### Solutions of univariate polynomial systems

We remember first that zero degree polynomials have no zeros (unless the are the zero polynomial which has value identical zero). In the following we implicitly assume that any polynomial we are talking about is not the zero polynomial.

Then we also know that any polynomial of degree one has exactly one zero. This is valid for any field, including $\R$ and $\C$.

The following is an application of a theorem for continuous functions as polynomials are continuous.

**Proposition:** Any real, univariate polynomial which has both positive and negative function values has at least one real zero.

A direct consequence is

**Corollary:** Any real univariate polynomial of odd degree has at least one real zero.

However, not all real univariate polynomials have real zeros, a simple example is $p(x) = x^2 +1$. A fundamental result of algebra, however, is the following.

**Theorem:** Any complex polynomial of degree larger than zero has at least one (complex) zero.

Note, however, that a system with $n>1$ univariate polynomial equations may not have a solution. For example, the system
$$\begin{align*} x^2 - 1 &= 0 \\ x^3-1 &= 0\end{align*}$$
has the solution $x=1$ while the system
$$\begin{align*} x^2 - 4 &= 0\\ x^3-1 &= 0\end{align*}$$
does not have a solution.

The nature of the set of solutions is provided now:

**Proposition:**
The set of all solutions of a system of real or complex univariate polynomials is finite and may have at most $d$ elements where $d$ is the degree of the system.

Summarising, for the case $s=1$ one has a good understanding of the set of solutions using results from algebra and calculus.

#### Solutions of linear systems of equations
  In this case of first degree polynomials in $\P^s_1$ and the polynomial system is given by
  $$P(x) = \begin{bmatrix}A & -b \end{bmatrix} \begin{bmatrix}x \\ 1\end{bmatrix} = Ax - b.$$ The matrix is $A\in\R^{n,s}$. Now let $$\range(A) = \{Ax \mid x \in \R^s\}$$ be the *range* of the matrix $A$. By definition one has

*A linear system $Ax-b=0$ has a solution $x$ if $$b\in \range(A).$$*

A general $b\in\R^n$ can be decomposed as
$$b = b_1 + b_2$$
where $b_1\in\range(A)$ and $b_2\in\range(A)^\perp$, i.e., $b_2$ is orthogonal to the range. It follows that we have a solution only if $b_2=0$.

Now for any $y\in \range(A)^\perp$ one has
$$(Ax)^T y = 0$$
for all $x\in\R^s$ and thus $y$ has to satisfy the homogeneous
system of equations
$$A^T y = 0.$$ The set of all such $y$ is called the *null-space* of $A$
$$\null(A^T) = \{y \mid A^T y = 0\}$$
and thus
$$\range(A)^\perp = \null(A^T).$$

If the nullspace of $A$ is nontrivial (not equal $\{0\}$) then the determinant of $A$ is zero $$\det(A)=0.$$ Note that the determinant is a polynomial in the matrix elements and thus the nullspace is a solution of a particular polynomial equation.

If $b\in\range(A)$ it has to be orthogonal to $\range(A)^\perp$ and one gets the major result of linear algebra:

**Proposition:** The linear system $Ax=b$ *has a solution* if $b$ is orthogonal to any vector $y\in \null(A^T)$.

The null space $\null(A)$ can be shown to be a linear space with dimension $s-r$ where $r$ is the rank of the matrix $A$ (which basically is defined by the dimension of the null space). One can see that if $x$ is a solution of $Ax=b$ and $y\in\null(A)$ then $x+y$ is also a solution as $A(x+y) = b$. One then has a characterisation of all the solutions of $Ax=b$:

**Proposition:** The set of solutions of $Ax=b$ is an *affine set* given by $$\{x+y \mid Ay = 0\}$$ for any solution $x$ of $Ax=b$.

The original characterisation $Ax+b=0$ of the solutions and the proposition above give an *implicit* description of the set of solutions. As the null space of $A$ is a linear or vector space, there exists a matrix $Y\in\R^{s,m}$ such that $AY=0$ and where the columns of $Y$ are just the basis of the null space of $A$. It then follows that the set of solutions of $Ax -b=0$ can be represented as $$\{x+Yt \mid t\in \R^m\}.$$ Note that $m$ is the dimension of the null space of $A$. We call this representation *explicit*.

As in the case of univariate polynomials one also has a good idea of what the solutions are for linear systems of equations. The situation is not as simple for higher degree multivariate polynomials.

#### Inverse function theorem

The next theorem reduces the theory of nonlinear problems to linear problems for the case $s=n$. Here we consider continuously differentiable functions $F$. We will use the *Jacobi matrix* which we denote by $$J_F(x)=F^\prime(x)=\begin{bmatrix}\partial F_1(x)/\partial x_1 & \cdots & \partial F_1(x)/\partial x_s\\ \vdots & & \vdots \\ \partial F_n(x)/\partial x_1 & \cdots & \partial F_n(x) / \partial x_s\end{bmatrix}.$$

**Inverse Function Theorem:**
*Let $F: \R^s \rightarrow \R^s$ be $C^1$ in a neighborhood of some point $x_0\in \R^s.$ If the null space $$\null(F^\prime(x_0))=0$$ then $F$ is invertible in a neighborhood of $x_0$.*

This means that the equation $$F(x) = F(x_0) + y$$ has a unique solution $x$ for all $y$ in a neigborhood of $F(x_0)$.

#### Solutions of polynomial systems of equations

This is related to the section 1.2 of Stetter's book.

We will denote the Jacobian of a polynomial by $P^\prime(x) = J_P(x)$.

More generally, we introduce differential operators related to polynomials. Let the gradient be
$$\partial = \nabla = (\partial_1,\ldots,\partial_s).$$
(A row vector.) Then let
$$\partial_j = \nabla^j = \frac{1}{j!}\frac{\partial^{|j|}}{\partial x^j}$$
where

* $j! = j_1! \cdots j_s!$
* $|j|=j_1+\cdots+j_s$
* $\partial x^j = \partial x_1^{j_1}\cdots \partial x_s^{j_s}$

Furthermore, let $q(x) = \sum_j b_j x_j$ then we define $$q(\partial) = \sum_j b_j \partial_j.$$

Then one has

**Proposition:** (Leibniz rule)
$$\partial_j (p\cdot q) = \
\sum_{k\leq j} \partial_{j-k}p\, \partial_k q$$

**Proposition:** (derivative of polynomials)

*If $p\in\P_d^s$ then*

* *$\partial_j p \in \P_{d-|j|}^s$ for $|j| \leq d$*
* *$\partial_j p = 0$ for $|j|>d$*

**Proposition:** (Taylor)

$$p(x+h) = \sum_{k=0}^d \sum_{|j|=k} \partial_j p(x)\, h^j =: p(h; x)$$


**Proposition:** (zero polynomial)
*A polynomial $p(x)$ defines the zero mapping if and only if it is the zero polynomial (i.e. has only zero coefficients).*

**proof**: induction over $s$:

* show for $s=1$
* $p(x) = p(x_1,...,x_{s-1},x_s)$ is univariate in $x_s$ ...

**Proposition:** (existence in $\C^s$)
*Any polynomial $p\in \P^s$ has at least one zero in $\C^s$ if $s>0$.*

**proof:** induction as before

**Definition:** A polynomial system $P(x)$ is *essentially linear* if $\det(P^\prime(x))$ is constant.

**Proposition:**
*If $n=s$, a polynomial system $P(x)$ has a non-invertible Jacobian at least in one point $x\in \C^s$ unless $P(x)$ is essentially linear.*

**proof:** induction as before and note that determinant of Jacobian is a polynomial

## 1.3 Floating point

* relative error model

## 1.4 Random polynomials