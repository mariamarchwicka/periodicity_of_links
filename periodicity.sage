#!/usr/bin/env python

# Copyright (c) 2018: Maria Marchwicka, Wojciech Politarczyk.
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>



import sys
import os
import numpy as np
import warnings



class MySettings(object):

    def __init__(self):

        self.f_pd_knot_11_15 = os.path.join(os.getcwd(), "knots_11_15.txt")
        self.f_knot_up_to_10 = os.path.join(os.getcwd(), "knots_3_10.txt")

        self.f_homfly_lm_in = os.path.join(os.getcwd(), "homflypt.input")

        self.f_results_out = os.path.join(os.getcwd(), "results.out")

        self.periods = [3, 5, 7, 9, 11]
        self.set_to_check = self.get_set()

        # check only knots from defined set
        self.only_chosen = True
        self.only_chosen = False


        self.print_results = False
        self.print_results = True

        # HOMFLYPT polynomials from file
        self.input_file_with_homflypt = True
        # self.input_file_with_homflypt = False


        if self.input_file_with_homflypt:
            if not os.path.isfile(self.f_homfly_lm_in):
                warnings.warn("No input file with HOMFLYPT polynomials")
                self.input_file_with_homflypt = False

    def get_set(self):

        set_to_check = set()
        return set_to_check


class PeriodicityTester(object):

    def __init__(self, name, pd_code, A=None, f_homfly_in=None):

        self.results = []
        '''
        To results for each period q a list in following form will be appended:
        [q, murasugi, naik_1, naik_2, borodzik, przytycki]
        Crierion is set to be:
        -1 if it is not applicable (details in check_naik_2, check_przytycki,
        1 if criterion doesn't exclude periodic,
        0 if criterion excludes periodicity.
        murasugi, naik_1, naik_2 or borodzik is also set to be:
        2 if alexander_polynomial == 1.
        0 if previous criterion in the list is 0.
        '''

        self.name = name
        self.pd_code = pd_code

        self.smith = None
        self.reset_results()

        if pd_code is not None:
            self.K = Link(pd_code)
            self.seifert = self.K.seifert_matrix()
        else:
            self.seifert = A
        # delta := Alexander polynomial
        delta = (self.seifert.transpose() - t * self.seifert).determinant()
        self.delta = delta.shift(-delta.exponents()[0])
        self.delta_factors = self.set_delta_factors()
        self.przytycki_tester = self.get_przytycki_tester(f_homfly_in)

    def reset_results(self):
        self.murasugi = 0
        self.naik_1 = 0
        self.naik_2 = 0
        self.borodzik = 0
        self.przytycki = 0
        self.murasugi_fulfilling = set()
        self.naik_1_fulfilling = []
        self.naik_2_fulfilling = []

    def set_smith(self):
        symetric_from_seifert = self.seifert + self.seifert.transpose()
        assert symetric_from_seifert.determinant() != 0, \
               "The determinant of A + A^T is zero."
        self.smith = symetric_from_seifert.smith_form()
        D, U, V = self.smith
        self.diagonal = D.diagonal()
        self.maximum_in_diagonal = max(self.diagonal)
        C = U.inverse()
        E_inverse = V
        self.C_tran_E_inv_D_inv = C.transpose() * E_inverse * D.inverse()
        self.matrix_C = C
        self.matrix_E_inverse = E_inverse

    def get_przytycki_tester(self, f_homfly_in):
        if self.pd_code is not None:
            try:
                return PrzytyckiTester(self.K, self.name, f_homfly_in)
            except ImportError as e:
                pass
        return None

    def get_C_tran_E_inv_D_inv(self):
        if self.smith is None:
            self.set_smith()
        return self.C_tran_E_inv_D_inv

    def get_maximum_in_diagonal(self):
        if self.smith is None:
            self.set_smith()
        return self.maximum_in_diagonal

    def set_delta_factors(self):
        # find all delta (alexander polynomial) factors
        lst_of_factors = [[f[0]] * f[1] for f in self.delta.factor()]
        # flattening a list
        lst_of_factors = [el for sublist in lst_of_factors for el in sublist]
        delta_candidates = set()
        for s in get_subsets(lst_of_factors):
            d = t^0
            for el in s:
                d *= el
            delta_candidates.add(d)
        return delta_candidates

    def check_criteria_for_period(self, q):

        self.reset_results()
        self.przytycki = self.check_przytycki(q)

        if self.delta == 1:
            self.murasugi = 2
            self.naik_1 = 2
            self.naik_2 = 2
            self.borodzik = 2
            return 2

        self.murasugi = self.check_murasugi(q)
        self.naik_1 = self.check_naik_1(q)
        self.naik_2 = self.check_naik_2(q)
        self.borodzik = self.check_borodzik(q)

        return self.borodzik * self.przytycki

    def check_murasugi(self, q):
        '''
        Select these delta factors and natural number r such that:
        delta = delta_prime^q * (1 + t^1 + ... + t^(r-1))^(q-1) mod q
        where "delta_prime" is a delta factor.
        '''
        quotient_delta = self.delta.change_ring(GF(q))
        # Underlying polynomial of quotient_delta:
        quotient_delta = quotient_delta.polynomial_construction()[0]
        delta_degree = quotient_delta.degree()

        for candidate in self.delta_factors:
            quotient_candidate = candidate.change_ring(GF(q))
            power_candidate = quotient_candidate^q
            power_candidate = power_candidate.polynomial_construction()[0]
            # (r - 1) - possible t-polynomial degree
            r = (delta_degree - power_candidate.degree()) / (q - 1) + 1
            if r < 1 or not r.is_integer():
                continue
            t_polynomial = get_t_polynomial(q, r)
            right_side = t_polynomial * power_candidate
            if quotient_delta != right_side and -quotient_delta != right_side:
                continue
            self.murasugi_fulfilling.add((candidate, r))

        return int(bool(self.murasugi_fulfilling))

    def check_naik_1_candidate(self, delta_prime, delta_evaluated, q):

        t_delta = delta_evaluated / delta_prime(-1)
        t_delta_dict = {f[0]: f[1] for f in factor(t_delta)}
        t_delta_factors = [f for f in t_delta_dict.keys()
                           if f != 2 and gcd(q, f) == 1]
        for f in t_delta_factors:
            f_q = naik_number_dict.setdefault((f, q), get_naik_number(f, q))
            if not (t_delta_dict[f] / (2 * f_q)).is_integer():
                return None
        return t_delta_factors

    def check_naik_1(self, q):
        '''
        For each delta' find a set P of prime numbers p such that:
        gcd(p, q) == 1, p != 2 and p| t_delta, t_delta = delta(-1)/delta'(-1).
        Check if all p factors of t_delta has multiplicity divisible by 2*[p|q].
        If it holds for at least one delta' candidate, set naik_1 = True.
        '''
        delta_evaluated = self.delta(-1)

        for delta_prime, _ in self.murasugi_fulfilling:
            t_delta_factors = self.check_naik_1_candidate(delta_prime,
                                                          delta_evaluated, q)
            if t_delta_factors is not None:
                self.naik_1_fulfilling.append((delta_prime, t_delta_factors))

        return int(bool(self.naik_1_fulfilling))

    def check_naik_2_candidate(self, q, p_list):
        delta_prime_bases = []
        maximum_in_diagonal = self.get_maximum_in_diagonal()
        for p in p_list:
            p_q = naik_number_dict[(p, q)]
            bases_for_p_torsion = []
            factor_power = p
            # find all p^k torsion parts
            while (maximum_in_diagonal / factor_power).is_integer():
                basis_for_p_k_part = []
                for el in self.diagonal:
                    to_be_append = el / factor_power
                    is_int = (to_be_append / p).is_integer()
                    if to_be_append.is_integer() and not is_int:
                        basis_for_p_k_part.append(to_be_append)
                    else:
                        basis_for_p_k_part.append(0)
                len_non_zero = sum(x != 0 for x in basis_for_p_k_part)
                # check if dimension is multiple of 2 * naik_number
                if not (len_non_zero / (2 * p_q)).is_integer():
                    return None
                factor_power *= p
                bases_for_p_torsion.append(basis_for_p_k_part)
            delta_prime_bases.append((p, bases_for_p_torsion))
        return delta_prime_bases

    def check_naik_2(self, q):
        '''
        For each delta' consider a set P of primes p such that: gcd(p, q) == 1,
        p != 2, p| delta(-1)/delta'(-1) (self.naik_1_fulfilling) and p is not
        a factor of delta'(-1). Check if dimension of p^k torsion part
        is divisible by 2*[p|q] for all k and all p from P.
        If it holds for at least one delta' candidate, we set naik_2 to be True.
        In particular naik_2 is set to be -1 if the criterion passes,
        but only in cases where P is an empty set.
        '''
        for delta_prime, p_list in self.naik_1_fulfilling:
            delta_prime_factors = set([d[0] for d in factor(delta_prime(-1))])
            p_list = [p for p in p_list if p not in delta_prime_factors]

            if not p_list:
                self.naik_2 = -1
                self.borodzik = -1
                continue

            delta_prime_bases = self.check_naik_2_candidate(q, p_list)
            if delta_prime_bases is not None:
                self.naik_2_fulfilling.append((delta_prime,
                                               delta_prime_bases))
        if self.naik_2_fulfilling:
            return 1
        return self.naik_2

    def check_borodzik(self, q):
        '''
        Consider all delta' that meet criterion Naik 2.
        For all p from a set P (defined as in check_naik_2)
        and all k consider p^k torsion part.
        For each p^k torsion check if eta == epsilon_1 * epsilon_2
        (see check_borodzik_candidate()).
        If it holds for at least one delta' candidate, set borodzik to be True.
        In particular borodzik is set to be -1 if the criterion passes,
        but only in cases where P is an empty set.
        '''

        for delta_prime, delta_prime_bases in self.naik_2_fulfilling:
            borodzik_pass = True
            for p, bases_for_p in delta_prime_bases:
                # if len(bases_for_p) > 1:
                #     print "HURA"  # more than one p^k part - not found yet
                if not self.check_borodzik_candidate(q, p, bases_for_p):
                    borodzik_pass = False
                    break
            if borodzik_pass:
                return 1
        return self.borodzik

    def check_borodzik_candidate(self, q, p, bases):
        '''
        For each p^k torsion check if eta == epsilon_1 * epsilon_2.
        If determinant of corsesponding matrix P is square modulo p, then:
        episilon_1 = 1, else: episilon_1 = -1.
        If p == 3 mod(4) and a rank of p^k torsion part n == 2 mod(4), then:
        epsilon_2 = -1, else: epsilon_2 = 1.
        eta = naik_sign ^ d, where  d = n / (2 * [p, q]).
        If p^([p, q]) % q == 1, then: naik_sign = 1, else: naik_sign = -1.
        '''
        for k, p_k_basis in enumerate(bases):
            X = np.diagflat(p_k_basis)
            # columns that up to zero (element in diagonal is zero):
            zero_columns = np.nonzero(X.sum(axis=0) == 0)
            X = np.delete(X, zero_columns, axis=1)
            n = X.shape[1]
            X = matrix(X)
            P = p^(k + 1) * X.transpose() * self.get_C_tran_E_inv_D_inv() * X
            P_det = P.determinant()
            if P_det % p == 0:
                raise ValueError("P determinant is 0 modulo p.")

            if p % 4 == 3 and n % 4 == 2:  # epsilon_1
                epsilon = -1
            else:
                epsilon = 1

            if not mod(P_det, p).is_square():
                epsilon *= -1  # epsilon = epsilon_1 * epsilon_2

            p_q = naik_number_dict[(p, q)]
            d = n / (2 * p_q)
            # sign(p_q) - whether rest is -1 or 1
            if sign(p_q)^d != epsilon:
                return False
        return True

    def check_przytycki(self, q):
        if self.przytycki_tester is not None and q in prime_numbers:
            try:
                return self.przytycki_tester.check_congruence(q)
            except (AttributeError, OverflowError) as e:
                pass
        return -1

    def save_results(self, f_out):

        for result in self.results:
            line_to_write = self.name + "," + ",".join(map(str, result))
            f_out.writelines(line_to_write + "\n")

    def print_results(self):

        print "\n" + "#" * 15 + " " + str(self.name) + " " + "#" * 15

        for result in self.results:

            q = result[0]
            print
            self.print_przytycki_result(q, result[5])

            if result[1] == 2:
                print "Alexander polynomial is 1"
                continue

            if not result[1]:
                print "\t\tMurasugi: fail, q = " + str(q)
                continue

            print "Murasugi: pass, q = " + str(q)

            if not result[2]:
                print "\t\tNaik 1: fail, q = " + str(q)
                continue

            print "Naik 1: pass, q = " + str(q)

            if not result[3]:
                print "\t\tNaik 2: fail, q = " + str(q)
                continue

            if result[3] == -1:
                print "Naik 2: not applicable, q = " + str(q)
                continue

            print "Naik 2: pass, q = " + str(q)

            if not result[4]:
                print ("\t\tBorodzik: fail, q = " + str(q))
                continue

            if result[4] == -1:
                print ("Borodzik: not applicable, q = " + str(q))
                continue

            print ("Borodzik: pass, q = " + str(q))

    def print_przytycki_result(self, q, result):
        if not result:
            print "\t\tPrzytycki: fail, q = " + str(q)
        elif result == -1:
            print "Przytycki: not applicable, q = " + str(q)
        else:
            print "Przytycki: pass, q = " + str(q)


class PrzytyckiTester(object):

    def __init__(self, K, name, f_homfly_in=None):
        homflypt = self.get_homflypt_polynomial(K, name, f_homfly_in)
        homfly_difference = homflypt(a, -z) - homflypt(a^-1, -z)
        self.homfly_difference = z * homfly_difference
        self.homflypt_polynomial = homflypt

    def get_homflypt_polynomial(self, K, name, f_homfly_in=None):
        if f_homfly_in is not None:
            try:
                current_name, homflypt = f_homfly_in.readline().split(',')
                while current_name != name:
                    current_name, homflypt = f_homfly_in.readline().split(',')
                homflypt = sage_eval(homflypt, locals={'a': a, 'z': z})
                return homflypt
            except (AttributeError, ValueError) as e:
                pass
        return K.homfly_polynomial('a', 'z', 'lm')

    def check_congruence(self, q):
        for i in range(q + 1):
            z_coefficient = self.homfly_difference.coefficient(z^(i+1))
            ideal = (a + a^-1)^(q - i)  # for i == q will be 1
            coefficient_modulo_ideal = z_coefficient.quo_rem(ideal)[1]
            coefficient_modulo_q = coefficient_modulo_ideal.change_ring(GF(q))
            if coefficient_modulo_q != 0:
                return 0
        return 1


def check_criteria(name, pd_code, f_homfly_in=None):


    tester = PeriodicityTester(name, pd_code, None, f_homfly_in)

    for i, q in enumerate(settings.periods):

        tester.check_criteria_for_period(q)
        tester.results.append([q, tester.murasugi, tester.naik_1,
                               tester.naik_2, tester.borodzik,
                               tester.przytycki])
    if settings.print_results:
        tester.print_results()

    return tester


def get_naik_number(p, q):
    '''
    Calculate the smallest integer i such that p^i == +/-1 mod q.
    Signum of i shows whether rest is -1 or 1
    '''
    if gcd(q, p) > 1:
        return 0
    p_power = p
    for i in xrange(1, sys.maxint):
        pq = p_power % q
        if pq == 1:
            return i
        if pq == q - 1:
            return -i
        p_power *= p


def get_t_polynomial(q, r):  # for check_murasugi(), r coresponds to l in paper
    t_polynomial = sum([t^i for i in range(r)])
    t_polynomial = t_polynomial.change_ring(GF(q))
    t_polynomial ^= (q - 1)
    return t_polynomial


def get_subsets(myset):
    return reduce(lambda z, x: z + [y + [x] for y in z], myset, [[]])


def parse_pd_code(pd_code_from_file):
    set = '0987654321[],'
    pd_code = ''.join([c for c in pd_code_from_file if c in set])
    return eval(pd_code)


def parse_knot_name(name):
    data = name[5: -2].split(',')
    name = data[0].strip() + data[1].strip().lower()[:1] + data[2].strip()
    return name


def check_11_to_15(f_out, f_homfly_out=None, f_homfly_in=None):
    with open(settings.f_pd_knot_11_15, 'r') as f:
        line = f.readline()
        while line:
            name = parse_knot_name(line)
            pd_code = parse_pd_code(f.readline())
            line = f.readline()
            tester = check_criteria(name, pd_code, f_homfly_in)
            if tester is None:
                continue
            tester.save_results(f_out)


def check_up_to_10(f_out, f_homfly_in=None):
    with open(settings.f_knot_up_to_10, 'r') as f:
        line = f.readline()
        while line:
            line = line.split(" = ")
            name = str(line[0])[5:]
            pd_code = parse_pd_code(str(line[1]))
            line = f.readline()
            tester = check_criteria(name, pd_code, f_homfly_in)
            if tester is None:
                continue
            tester.save_results(f_out)


def test_all(f_out, f_homfly_in=None):
    check_up_to_10(f_out, f_homfly_in)
    if f_out is not None:
        f_out.flush()
    check_11_to_15(f_out, f_homfly_in)


if __name__ == '__main__':

    settings = MySettings()
    S.<a, z> = LaurentPolynomialRing(ZZ)
    R.<t> = LaurentPolynomialRing(ZZ)
    prime_numbers = Primes()
    naik_number_dict = {}
    if settings.input_file_with_homflypt:
        with open(settings.f_results_out, 'w') as f_out,\
             open(settings.f_homfly_lm_in, 'r') as f_homfly_in:
            test_all(f_out, f_homfly_in)
    else:
        with open(settings.f_results_out, 'w') as f_out:
            test_all(f_out)
