import numpy as np
from numpy import sin, cos, exp, pi, sqrt, inf
from scipy.integrate import quad
import math 
import scipy as scipy


def integrator(function, a, b):
    step = 0.001
    t = np.linspace(a, b, int((b - a) / step))
    return sum([function(i) * step for i in t])


def factorial(n):
    if n < 0:
        return 0
    if n == 0:
        return 1
    if n <= 2:
        return n

    ans = 1
    for i in range(2, n+1):
        ans *= i
    return ans


def F(x):
    return quad(lambda t: 1/sqrt(2 * pi) * exp(-(t**2)/2), -inf, x)[0]


def Moivre_Laplace_local_approx(n, p, q, k):
    x = (k - n * p)/sqrt(n * p * q)
    try:
        ans = 1/sqrt(n * p * q * 2 * pi) * exp((-x**2)/2)
    except:
        return "Too big nums!"
    if ans < math.e**(-14):
        ans = 0.0
    return ans


def Moivre_Laplace_integral_approx(n, p, q, a, b):
    x_1 = (b - n * p) / sqrt(n * p * q)
    x_2 = (a - n * p) / sqrt(n * p * q)
    try:
        ans = F(x_1) - F(x_2)
    except:
        return "Too big nums!"
    if ans < math.e**(-14):
        return 0.0
    return ans


def Pois_approx(n, p, q, k):
    lambda_pois = n * p 
    try:
        ans = exp(-lambda_pois) * (lambda_pois**k)/factorial(k)
    except:
        return "Too big nums!"
    if (ans < math.e**-14):
        return 0.00
    else:
        return ans

def classic_calculation(n, p, q, k):
    try:
        comb = factorial(n)/(factorial(k) * factorial(n - k))
        ans = comb * p**(k) * q**(n - k)
    except:
        return "Too big nums!"
    if (ans < math.e**-14):
        return 0.00
    else:
        return ans


def S_n_btw_classical(n, p, q):
    a = int(n/2 - sqrt(n * p * q))
    b = int(n/2 + sqrt(n * p * q))
    ans = 0
    for k in range(a, b+1):
        try:
            ans += classic_calculation(n, p, q, k)
        except:
            return "Too big nums"
    if (ans < math.e**(-14)):
        return 0
    return ans


def S_n_btw_local(n, p, q):
    a = int(n/2 - sqrt(n * p * q))
    b = int(n/2 + sqrt(n * p * q))
    ans = 0
    for i in range(a, b+1):
        try:
            ans += Moivre_Laplace_local_approx(n, p, q, i)
        except:
            return "Too big nums"

    return ans


def S_n_btw_pois(n, p, q):
    a = int(n/2 - sqrt(n * p * q))
    b = int(n/2 + sqrt(n * p * q))
    ans = 0
    for i in range(a, b+1):
        try:
            ans += Pois_approx(n, p, q, i)
        except:
            return "Too big nums"
    return ans


def S_n_btw_integral(n, p, q):
    a = n/2 - sqrt(n * p * q)
    b = n/2 + sqrt(n * p * q)
    return Moivre_Laplace_integral_approx(n, p, q, a, b)


def S_n_less_classical(n, p, q, maximum_k):
    ans = 0
    for k in range(0, maximum_k + 1):
        comb = factorial(n) / (factorial(k) * factorial(n-k))
        try:
            ans += comb * p**k * q**(n-k)
        except:
            return "Too big nums"
    return ans 


def S_n_less_local(n, p, q, maximum_k):
    ans = 0
    for k in range(0, maximum_k + 1):
        try:
            ans += Moivre_Laplace_local_approx(n, p, q, k)
        except:
            return "Too big nums."
    return ans


def S_n_less_pois(n, p, q, maximum_k):
    ans = 0
    for k in range(0, maximum_k + 1):
        try:
            ans += Pois_approx(n, p, q, k)
        except:
            return "Too big numbers!!"

    return ans


def S_n_less_integral(n, p, q, maximum_k):
    return Moivre_Laplace_integral_approx(n, p, q, 0, maximum_k)


def S_n_likely(n, p, q):
    guess = (n + 1) * p 
    if (guess - int(guess) != 0):
        return math.ceil(guess)
    return guess - 1


if __name__ == "__main__":
    p_s = [0.001, 0.01, 0.1, 0.25, 0.5]
    q_s = [1 - i for i in p_s]
    n_s = [100, 1000, 10000]
    print("====================-")
    for n in n_s:
        for p in p_s:
            print(f"n = {n} p = {p} q = {1-p}")
            print(f'P(x1 <= S <= x2): ')
            print(f'Классическая формула | теорема Пуассона | Локальная теорема | Интегральная теорема')
            print(S_n_btw_classical(n, p, 1 - p), "|", S_n_btw_pois(n, p, 1 - p), "|", S_n_btw_local(n, p, 1 - p), "|", S_n_btw_integral(n, p, 1 - p))
            print(f'P(S <= 5 ): ')
            print(f'Классическая формула | теорема Пуассона | Локальная теорема | Интегральная теорема')
            print(S_n_less_classical(n, p, 1 - p, 5), "|", S_n_less_pois(n, p, 1 - p, 5), "|", S_n_less_local(n, p, 1 - p, 5), "|", S_n_less_integral(n, p, 1 - p, 5))
            most_likely_shot = S_n_likely(n, p, 1 - p) 
            print("k* =", S_n_likely(n, p, 1 - p))
            print("P(S = k*)")
            print(f"Классическая формула | теорема Пуассона | Локальная теорема | Интегральная теорема")
            print(classic_calculation(n, p, 1-p, most_likely_shot), Pois_approx(n, p, 1 - p, most_likely_shot), Moivre_Laplace_local_approx(n, p, 1-p, most_likely_shot), Moivre_Laplace_integral_approx(n, p, 1-p, 0, most_likely_shot))
            print("=====================")