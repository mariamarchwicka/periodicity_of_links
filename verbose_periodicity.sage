#!/usr/bin/python

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


if not os.path.isfile('knots_periodicity.py'):
    os.system('sage --preparse periodicity.sage')
    os.system('mv periodicity.sage.py knots_periodicity.py')
from knots_periodicity import get_t_polynomial, \
                        parse_pd_code, parse_knot_name
import knots_periodicity



class MySettings(object):

    def __init__(self):

        self.f_pd_knot_11_15 = os.path.join(os.getcwd(), "knots_11_15.txt")
        self.f_knot_up_to_10 = os.path.join(os.getcwd(), "knots_3_10.txt")

        self.f_homfly_lm_out = os.path.join(os.getcwd(), "homflypt.out")
        self.f_homfly_lm_in = os.path.join(os.getcwd(), "homflypt.input")

        self.f_results_out = os.path.join(os.getcwd(), "results.out")
        self.f_old_results = os.path.join(os.getcwd(), "old_results.input")

        # self.periods = [9]

        self.periods = [3, 5, 7, 9, 11]
        self.set_to_check = self.get_set()

        # check only knots from defined set
        self.only_chosen = True
        # self.only_chosen = False

        self.debugging = True
        # self.debugging = False

        # only if debugging
        self.print_matrices = True
        self.print_matrices = False

        # only if only_chosen
        self.only_periods_where_borodzik = True
        self.only_periods_where_borodzik = False

        # only if only_chosen
        self.only_periods = True
        self.only_periods = False

        self.print_results = True
        # self.print_results = False

        # saving HOMFLYPT polynomials into self.f_homfly_lm_out
        self.save_homfly = True
        self.save_homfly = False

        # reuse HOMFLYPT polynomials previously saved
        self.input_file_with_homflypt = True
        # self.input_file_with_homflypt = False

        self.check_old_results = True
        self.check_old_results = False

        if self.input_file_with_homflypt:
            if not os.path.isfile(self.f_homfly_lm_in):
                warnings.warn("No input file with HOMFLYPT polynomials")
                self.input_file_with_homflypt = False

    def get_set(self):
        set_to_check = set()

        periodic_burde = set(["3_1", "5_1", "7_1", "8_19", "9_1",
                              "9_35", "9_40", "9_41", "9_47", "9_49",
                              "10_3", "10_123", "10_124"])
        set_to_check |= periodic_burde

        # knots that fail Borodzik criterion
        self.fails_dict = {
                        "12a100": 3,
                        "12a348": 3,
                        "13a4648": 3,
                        "13n3659": 3,
                        "14a7583": 3,
                        "14a7948": 3,
                        "14a8670": 3,
                        "14a9356": 3,
                        "14a14971": 3,
                        "14a16311": 3,
                        "14a17173": 3,
                        "14a17260": 3,
                        "14a18647": 3,
                        "14n908": 3,
                        "14n913": 3,
                        "14n2451": 3,
                        "14n2458": 3,
                        "14n6565": 3,
                        "14n9035": 3,
                        "14n11989": 3,
                        "14n14577": 3,
                        "14n23051": 3,
                        "14n24618": 3,
                        "15a6030": 3,
                        "15a6066": 3,
                        "15a10622": 3,
                        "15a15077": 3,
                        "15a33910": 3,
                        "15a36983": 3,
                        "15a46768": 3,
                        "15a72333": 3,
                        "15a82771": 3,
                        "15n3147": 3,
                        "15n3369": 3,
                        "15n3372": 3,
                        "15n4496": 3,
                        "15n4514": 3,
                        "15n4517": 3,
                        "15n11293": 3,
                        "15n11533": 3,
                        "15n14173": 3,
                        "15n15251": 3,
                        "15n19351": 3,
                        "15n19989": 3,
                        "15n20691": 3,
                        "15n33684": 3,
                        "15n34725": 3,
                        "15n36715": 3,
                        "15n45612": 3,
                        "15n49287": 3,
                        "15n55026": 3,
                        "15n58771": 3,
                        "15n59908": 3,
                        "15n61622": 3,
                        "15n61790": 3,
                        "15n61833": 3,
                        "15n63397": 3,
                        "15n67585": 3,
                        "15n69848": 3,
                        "15n90233": 3,
                        "15n90525": 3,
                        "15n112198": 3,
                        "15n115648": 3,
                        "15n116414": 3,
                        "15n120198": 3,
                        "15n120375": 3,
                        "15n133302": 3,
                        "15n135864": 3,
                        "15n135918": 3,
                        "15n148509": 3,
                        "15n155150": 3,
                        "15n158831": 3,
                        "15n162066": 3,
                        "15n162237": 3,
                        "15n163140": 3,
                        "15n165092": 3,
                        "15n165622": 3,
                        "15n167650": 3,
                        "14n26993": 5,
                        "15a80526": 5,
                        "15n83514": 5,
                        "15n95792": 5,
                        }
        set_to_check |= set(self.fails_dict.keys())

        # knots that pass Borodzik criterion
        self.success_dict = {
                        "9_40": 3,
                        "9_41": 3,
                        "9_49": 3,
                        "11a297": 3,
                        "11a321": 3,
                        "11n133": 3,
                        "12a561": 3,
                        "12a780": 3,
                        "12a1019": 3,
                        "12a1202": 3,
                        "12a1206": 3,
                        "12n706": 3,
                        "12n837": 3,
                        "12n839": 3,
                        "12n843": 3,
                        "12n844": 3,
                        "12n847": 3,
                        "12n881": 3,
                        "13n2694": 3,
                        "14a2160": 3,
                        "14a7206": 3,
                        "14a10416": 3,
                        "14a12671": 3,
                        "14a15296": 3,
                        "14a16592": 3,
                        "14a18362": 3,
                        "14n945": 3,
                        "14n3276": 3,
                        "14n3888": 3,
                        "14n4912": 3,
                        "14n9288": 3,
                        "14n10327": 3,
                        "14n11309": 3,
                        "14n11898": 3,
                        "14n13447": 3,
                        "14n13863": 3,
                        "14n14083": 3,
                        "14n14183": 3,
                        "14n14497": 3,
                        "14n16414": 3,
                        "14n16415": 3,
                        "14n16428": 3,
                        "14n16682": 3,
                        "14n17032": 3,
                        "14n17183": 3,
                        "14n17871": 3,
                        "14n17959": 3,
                        "14n21996": 3,
                        "14n23568": 3,
                        "14n24905": 3,
                        "15a8033": 3,
                        "15a15545": 3,
                        "15a20833": 3,
                        "15a22423": 3,
                        "15a23751": 3,
                        "15a24566": 3,
                        "15a24687": 3,
                        "15a33565": 3,
                        "15a35274": 3,
                        "15a39992": 3,
                        "15a40971": 3,
                        "15a54610": 3,
                        "15a74206": 3,
                        "15a74381": 3,
                        "15a77993": 3,
                        "15a81135": 3,
                        "15a81151": 3,
                        "15a81179": 3,
                        "15a81370": 3,
                        "15a81477": 3,
                        "15a81796": 3,
                        "15a82451": 3,
                        "15a82698": 3,
                        "15a83361": 3,
                        "15a83814": 3,
                        "15a85128": 3,
                        "15a85145": 3,
                        "15a85169": 3,
                        "15a85223": 3,
                        "15a85254": 3,
                        "15a85257": 3,
                        "15n15450": 3,
                        "15n15810": 3,
                        "15n17487": 3,
                        "15n17658": 3,
                        "15n18682": 3,
                        "15n20353": 3,
                        "15n28777": 3,
                        "15n29526": 3,
                        "15n31070": 3,
                        "15n33167": 3,
                        "15n39609": 3,
                        "15n39756": 3,
                        "15n39792": 3,
                        "15n39829": 3,
                        "15n39838": 3,
                        "15n39866": 3,
                        "15n40188": 3,
                        "15n40203": 3,
                        "15n45334": 3,
                        "15n47776": 3,
                        "15n48100": 3,
                        "15n50732": 3,
                        "15n52424": 3,
                        "15n52723": 3,
                        "15n55025": 3,
                        "15n59277": 3,
                        "15n59777": 3,
                        "15n59987": 3,
                        "15n60899": 3,
                        "15n61859": 3,
                        "15n62066": 3,
                        "15n68367": 3,
                        "15n68469": 3,
                        "15n72570": 3,
                        "15n75241": 3,
                        "15n77155": 3,
                        "15n78784": 3,
                        "15n78786": 3,
                        "15n81011": 3,
                        "15n84209": 3,
                        "15n85291": 3,
                        "15n93105": 3,
                        "15n95263": 3,
                        "15n95294": 3,
                        "15n98814": 3,
                        "15n99593": 3,
                        "15n100351": 3,
                        "15n105142": 3,
                        "15n122147": 3,
                        "15n126255": 3,
                        "15n127744": 3,
                        "15n132539": 3,
                        "15n134183": 3,
                        "15n134435": 3,
                        "15n135170": 3,
                        "15n137023": 3,
                        "15n142082": 3,
                        "15n145384": 3,
                        "15n146140": 3,
                        "15n147033": 3,
                        "15n151780": 3,
                        "15n152852": 3,
                        "15n153976": 3,
                        "15n154660": 3,
                        "15n154766": 3,
                        "15n159959": 3,
                        "15n160415": 3,
                        "15n160533": 3,
                        "15n165706": 3,
                        "15n165708": 3,
                        "15n165735": 3,
                        "15n165748": 3,
                        "15n166307": 3,
                        "15n167633": 3,
                        "15n167645": 3,
                        "15n168004": 3,
                        "15n168014": 3,
                        "10_123": 5,
                        "14n7478": 5,
                        "15a40549": 5,
                        "15a53966": 5,
                        "15a64035": 5,
                        "15a69121": 5,
                        "15a76651": 5,
                        "15a84903": 5,
                        "15a85262": 5,
                        "15n35157": 5,
                        "15n113162": 5,
                        "15n142117": 5,
                        "14a19470": 7,
                        "15n162490": 7,
                        "15a74206": 9,
                        "15n154766": 9,
                        "15n160415": 9,
                        "15n165706": 9,
                        }
        set_to_check |= set(self.success_dict.keys())

        # knots that are known to be (not) periodic
        self.periods_dict = {
                            "3_1": [3],
                            "5_1": [5],
                            "7_1": [7],
                            "8_19": [3],
                            "9_1": [3, 9],
                            "9_35": [3],
                            "9_40": [3],
                            "9_41": [3],
                            "9_47": [3],
                            "9_49": [3],
                            "10_3": [3],
                            "10_123": [5],
                            "10_124": [3, 5],
                            "11a367": [11],
                            "12a503": [3],
                            "12a561": [3],
                            "12a615": [3],
                            "12a1019": [3],
                            "12a1022": [3],
                            "12a1202": [3],
                            "14a19470": [7],
                            "15a64035": [5],
                            "15a84903": [5],
                            "15a85262": [5],
                            "15a85263": [5],
                            "12a100": [-3],
                            "12a348": [-3],
                            "12a376": [-3],
                            "12a1206": [-3],
                            "13a2142": [-5],
                            "13a2907": [-5],
                            "13a3010": [-5],
                            "15a23599": [-5],
                            "15a23902": [-5],
                            "15a40549": [-5],
                            "15a53966": [-5]
                            }
        set_to_check |= set(self.periods_dict.keys())

        # knots that have Alexander polynomial  = 1
        self.alexander_1 = set(["11n34",
                                "11n42",
                                "12n313",
                                "12n430",
                                "13n65",
                                "13n71",
                                "13n866",
                                "13n1019",
                                "13n1496",
                                "13n1756",
                                "13n1757",
                                "13n3871",
                                "13n3872",
                                "13n3897",
                                "13n3934",
                                "13n3936",
                                "13n3938",
                                "13n4582",
                                "13n4591",
                                "14n3798",
                                "14n4425",
                                "14n5152",
                                "14n5486",
                                "14n6082",
                                "14n7469",
                                "14n7708",
                                "14n9023",
                                "14n9290",
                                "14n9684",
                                "14n9773",
                                "14n9882",
                                "14n10011",
                                "14n10119",
                                "14n10990",
                                "14n11063",
                                "14n11129",
                                "14n11515",
                                "14n11680",
                                "14n12763",
                                "14n14735",
                                "14n14833",
                                "14n15285",
                                "14n15581",
                                "14n18909",
                                "14n18911",
                                "14n21673",
                                "14n22185",
                                "14n22589",
                                "14n23325",
                                "14n23411",
                                "14n23417",
                                "14n23940",
                                "14n24036",
                                "14n24396",
                                "14n25281",
                                "15n2810",
                                "15n3240",
                                "15n4018",
                                "15n4646",
                                "15n11287",
                                "15n11296",
                                "15n11568",
                                "15n11570",
                                "15n15829",
                                "15n16056",
                                "15n17501",
                                "15n21288",
                                "15n21905",
                                "15n21939",
                                "15n21944",
                                "15n24436",
                                "15n25044",
                                "15n27582",
                                "15n27824",
                                "15n28998",
                                "15n29401",
                                "15n29559",
                                "15n29563",
                                "15n30723",
                                "15n31075",
                                "15n34773",
                                "15n36113",
                                "15n37062",
                                "15n38863",
                                "15n40132",
                                "15n40402",
                                "15n40938",
                                "15n42200",
                                "15n42279",
                                "15n42516",
                                "15n44873",
                                "15n45781",
                                "15n45782",
                                "15n46093",
                                "15n46536",
                                "15n48362",
                                "15n49081",
                                "15n49735",
                                "15n49992",
                                "15n50050",
                                "15n50051",
                                "15n50147",
                                "15n50819",
                                "15n51748",
                                "15n51847",
                                "15n52282",
                                "15n52651",
                                "15n54221",
                                "15n58433",
                                "15n58501",
                                "15n59917",
                                "15n59918",
                                "15n61482",
                                "15n62093",
                                "15n62150",
                                "15n63468",
                                "15n64468",
                                "15n65084",
                                "15n65980",
                                "15n71170",
                                "15n73226",
                                "15n74185",
                                "15n77245",
                                "15n77247",
                                "15n80534",
                                "15n82843",
                                "15n83995",
                                "15n85041",
                                "15n85314",
                                "15n87941",
                                "15n88033",
                                "15n89822",
                                "15n91092",
                                "15n91760",
                                "15n95983",
                                "15n95989",
                                "15n95995",
                                "15n96014",
                                "15n103703",
                                "15n108850",
                                "15n108966",
                                "15n110439",
                                "15n113775",
                                "15n115135",
                                "15n115375",
                                "15n117232",
                                "15n120055",
                                "15n120219",
                                "15n121343",
                                "15n121598",
                                "15n121834",
                                "15n121916",
                                "15n122603",
                                "15n123337",
                                "15n123414",
                                "15n123479",
                                "15n124496",
                                "15n124511",
                                "15n124640",
                                "15n124838",
                                "15n124849",
                                "15n125351",
                                "15n126042",
                                "15n126050",
                                "15n127500",
                                "15n128163",
                                "15n130096",
                                "15n130504",
                                "15n130528",
                                "15n131977",
                                "15n132396",
                                "15n132965",
                                "15n134216",
                                "15n135221",
                                "15n135706",
                                "15n138033",
                                "15n138051",
                                "15n139247",
                                "15n139256",
                                "15n139840",
                                "15n140327",
                                "15n140449",
                                "15n142843",
                                "15n143482",
                                "15n143825",
                                "15n143856",
                                "15n143985",
                                "15n144034",
                                "15n144436",
                                "15n144439",
                                "15n145339",
                                "15n145981",
                                "15n146982",
                                "15n151010",
                                "15n154389",
                                "15n155056",
                                "15n155464",
                                "15n156539",
                                "15n163337",
                                "15n165398",
                                ])
        set_to_check |= self.alexander_1

        set_to_check |= set(["10_123"])
        set_to_check |= set(["15n166130"])

        return set_to_check

class PeriodicityTester(knots_periodicity.PeriodicityTester):
    def save_results(self, f_out, f_homfly_out=None):
        for result in self.results:
            line_to_write = self.name + "," + ",".join(map(str, result))
            if settings.check_old_results and (result[0] in [3, 5, 7, 9, 11]):
                line = f_old_results.readline()
                # name, q, murasugi, naik_1, naik_2, borodzik, przytycki
                while line:
                    line = line.split(',')
                    if line[0] == self.name and line[1] == str(result[0]):
                        old_results = [int(x) for x in line[2:]]
                        # if old_results[:-1] != result[1:-1]:
                        if old_results[:] != result[1:]:
                            print("#" * 30 + " ERROR " + line[0] + " " +
                                   "#" * 30)
                            print("q = " + line[1])
                            print("result      " + str(result[1:]))
                            print("old_results " + str(old_results))
                        break
                    line = f_old_results.readline()
                if not line:
                    print("No data to compare.")
            f_out.writelines(line_to_write + "\n")

        if self.przytycki_tester is not None and f_homfly_out is not None:
            lm_polynomial = self.przytycki_tester.homflypt_polynomial
            line_to_write = self.name + "," + str(lm_polynomial) + "\n"
            f_homfly_out.writelines(line_to_write)

    def print_results(self):
        super().print_results()
        if self.name in settings.periods_dict:
            print("periods: " + str(settings.periods_dict[self.name]))

    def print_przytycki_result(self, q, result):
        if not result:
            print("\t\tPrzytycki: fail, q = " + str(q))
        elif result == -1:
            print("Przytycki: not applicable, q = " + str(q))
        else:
            print("Przytycki: pass, q = " + str(q))

    def print_data_for_murasugi(self, q):

        if self.murasugi:
            print("\n" + "#" * 30 + " Knot " + str(self.name) +
                   " passes Murasugi condition for q = " +
                   str(q) + " " + "#" * 30)
        else:
            print("\nKnot " + str(self.name) +
                   " fails Murasugi condition for q = " + str(q))

        quotient_delta = self.delta.change_ring(GF(q))
        quotient_delta = quotient_delta.polynomial_construction()[0]
        print("delta: \t" + str(self.delta))
        if self.delta == 1:
            print("delta factors: " + str([1]))
        else:
            print("delta factors: " + str(self.delta.factor()))
        print("delta mod q =  " + str(quotient_delta))
        delta_degree = quotient_delta.degree()
        self.print_murasugi_fulfilling(q)
        # self.print_candidates_that_fail_murasugi(q)

    def print_murasugi_fulfilling(self, q):
        quotient_delta = self.delta.change_ring(GF(q))
        quotient_delta = quotient_delta.polynomial_construction()[0]
        delta_degree = quotient_delta.degree()
        print("\nNumber of candidates that pass Murasugi = " +
               str(len(self.murasugi_fulfilling)))
        for i, (delta_prime, r) in enumerate(self.murasugi_fulfilling):
            print("\n" + str(i + 1) + ". delta_prime:\t" + str(delta_prime))
            t_polynomial = get_t_polynomial(q, r)
            print("polynomial^(q-1) = " + str(t_polynomial))
            right_side = t_polynomial * delta_prime^q
            print("*" * 50)
            print("delta == delta_prime^q * polynomial^(q-1) mod q")
            # print("right side:\t" + str(right_side.factor()))
            # print("left side: \t" + str(quotient_delta.factor()))

    def print_candidates_that_fail_murasugi(self, q):
        quotient_delta = self.delta.change_ring(GF(q))
        quotient_delta = quotient_delta.polynomial_construction()[0]
        delta_degree = quotient_delta.degree()
        for candidate in self.delta_factors:
            quotient_candidate = candidate.change_ring(GF(q))
            power_candidate = quotient_candidate^q
            shifted_candidate = power_candidate.polynomial_construction()[0]
            r = (delta_degree - shifted_candidate.degree()) / (q - 1) + 1
            if r > 0 and r.is_integer():
                t_polynomial = get_t_polynomial(q, r)
                right_side = t_polynomial * shifted_candidate
                if (quotient_delta == right_side or
                        (-quotient_delta) == right_side):
                    continue
            print("\nFor candidate =     " + str(candidate))
            print( "quotient_candidate = " + str(quotient_candidate))
            print("candidate^q =       " + str(power_candidate))
            print("shifted           = " + str(shifted_candidate))
            print("delta degree = " + str(delta_degree))
            print("candidate^q degree " + str(shifted_candidate.degree()))
            print("r = " + str(r))
            if r > 0 and r.is_integer():
                print("right_side = " + str(right_side))
                print("delta mod q = " + str(quotient_delta))

    def print_data_for_naik_1(self, q):
        if not self.murasugi:
            return None
        if not self.naik_1:
            print("\nKnot " + str(self.name) +
                   " fails Naik 1 condition for q = " + str(q))
        else:
            print("\n" + "#" * 30 + " Knot " + str(self.name) +
                   " passes Naik 1 condition for q = " + str(q) +
                   " " + "#" * 30)
        print("delta:       " + str(self.delta))
        print("delta at -1: " + str(self.delta(-1)))
        print("factors for evaluated: " + str(self.delta(-1).factor()))
        self.print_naik_1_fulfilling(q)

    def print_naik_1_fulfilling(self, q):
        print("\nNumber of candidates that pass Naik 1 = " +
               str(len(self.naik_1_fulfilling)))
        for delta_prime, p_list in self.naik_1_fulfilling:
            print("delta prime:        " + str(delta_prime))
            print("delta prime at -1:    " + str(delta_prime(-1)))
            t_delta = self.delta(-1)/delta_prime(-1)
            print("delta/delta_prime(-1):\t\t" + str(t_delta))
            print("delta/delta_prime(-1) factors:\t" + str(t_delta.factor()))
            if not p_list:
                print("List of factors was empty.")
            for p in p_list:
                g = abs(naik_number_dict[(q, p)])
                print("factor of del/del'(-1): " + str(p))
                print("Naik number: " + str(g))
                print("2 * Naik number:\t" + str(2 * g))
                test_naik_number = p^g % q
                print(str(p) + "^" + str(g) + " % " + str(q) + " = " +
                      str(test_naik_number) + " = " +
                      str(test_naik_number - q))
                t_delta_dict = {i[0]: i[1] for i in factor(t_delta)}
                print("The power of factor:\t" + str(t_delta_dict[p]))

    def print_data_for_naik_2(self, q):
        if not self.naik_1:
            return None
        if not self.naik_2:
            return None
        print("\n" + "#" * 30 + " Knot " + str(self.name) +
               " passes Naik 2 condition for q = " + str(q) + " " + "#" * 30)
        print("delta:\t\t\t" + str(self.delta))
        print("delta at -1:\t\t" + str(self.delta(-1)))
        print("factors for evaluated:\t" + str(self.delta(-1).factor()))
        if self.naik_2 == -1:
            self.print_naik_2_not_applicable(q)
            return None
        self.print_naik_2_fulfilling(q)

    def print_naik_2_not_applicable(self, q):
        for delta_prime, p_list in self.naik_1_fulfilling:
            delta_prime_factors = set([d[0] for d in factor(delta_prime(-1))])
            p_list = [p for p in p_list if p not in delta_prime_factors]
            if not p_list:
                print("\nChecking Naik 2 condition for candidate " +
                       str(delta_prime) + " and q = " + str(q) + ".")
                print("The list of factors was empty or all factors " +
                       "were dela'(-1) factors.")
                print("Naik 2 and Borodzik can not exclude periodicity.\n")

    def print_naik_2_fulfilling(self, q):
        for delta_prime, delta_prime_bases in self.naik_2_fulfilling:
            print("\ndelta prime:\t\t\t" + str(delta_prime))
            print("delta prime at -1:\t\t" + str(delta_prime(-1)))
            t_delta = self.delta(-1)/delta_prime(-1)
            print("delta/delta_prime(-1):         " + str(t_delta))
            print("delta/delta_prime(-1) factors: " + str(t_delta.factor()))

            for p, bases_for_p in delta_prime_bases:
                print("\nfactor p for delta prime:\t\t\t" + str(p))
                g = abs(naik_number_dict[(q, p)])
                print("Naik number:\t\t" + str(g))
                print("2 * Naik number:\t" + str(2 * g))
                test_naik_number = p^g % q
                print(str(p) + "^" + str(g) + " % " + str(q) + " = " +
                      str(test_naik_number) + " = " +
                      str(test_naik_number - q))
                t_delta_dict = {i[0]: i[1] for i in factor(t_delta)}
                print("The power of factor:\t" + str(t_delta_dict[p]))
                print("diagonal: " + str(self.diagonal))
                print("p^k basis")
                for k, b in enumerate(bases_for_p):
                    print("k = " + str(k + 1))
                    print("basis:\t" + str(b))

    def print_data_for_borodzik(self, q):

        if self.naik_2 != 1:
            return None
        if self.borodzik:
            print("\n" + "#" * 30 + " Knot " + str(self.name) +
                  " passes Borodzik condition for q = " +
                  str(q) + " " + "#" * 30)
        else:
            print("%" * 200)
            print("\nKnot " + str(self.name) +
                  " fails Borodzik condition for q = " + str(q))

        if settings.print_matrices:
            self.print_matrices_for_borodzik(q)

        for delta_prime, delta_prime_bases in self.naik_2_fulfilling:
            print("\nResults for candidate delta_prime = " + str(delta_prime))
            for p, bases_for_p in delta_prime_bases:
                print("Results for p = " + str(p))
                for k, p_k_basis in enumerate(bases_for_p):
                    self.print_borodzik_for_p_k_basis(p, k, p_k_basis, q)
            print("%" * 200 + "\n" * 3)

    def print_matrices_for_borodzik(self, q):
        print("\n\nSeifert matrix A:")
        print(str(self.seifert))
        print("\n\nA + A^T:")
        print(str(self.seifert + self.seifert.transpose()))
        print("\n\nC")
        print(str(self.matrix_C))
        # print("\nE^(-1)")
        # print(str(self.E_inverse))
        print("\n\nD - diagonal")
        print(str(self.diagonal))
        print("\n\nE")
        print(str(self.matrix_E_inverse.inverse()))
        print("\n\nC^T * E^{-1} * D^{-1}")
        print(self.get_C_tran_E_inv_D_inv())

    def print_borodzik_for_p_k_basis(self, p, k, p_k_basis, q):

        # X matrix
        X = np.diagflat(p_k_basis)
        zero_columns = np.nonzero(X.sum(axis=0) == 0)
        X = np.delete(X, zero_columns, axis=1)
        n = X.shape[1]
        X = matrix(X)

        # P deterinant and epsilon_1
        P = p^(k + 1) * X.transpose() * self.get_C_tran_E_inv_D_inv() * X
        P_det = P.determinant()
        if settings.print_matrices:
            print("\nsubmatrix:")
            print(self.C_tran_E_inv_D_inv[-n:, -n:])
            print("\nP\n" + str(P))
        print("\ndet(P) = " + str(P_det))
        if mod(P_det, p).is_square():
            print("det(P) % p = " + str(P_det % p) +
                  " is a square => epsilon_1 := 1")
            epsilon_1 = 1
        else:
            print("det(P) % p = " + str(P_det % p) +
                  " isn't a square => episilon_1 := -1")
            epsilon_1 = -1

        # p % 4 and n % 4, and epsilon_2
        print("\np % 4 = " + str(p) + " % 4 = " + str(p % 4))
        print("n % 4 = " + str(n) + " % 4 = " + str(n % 4))
        if p % 4 == 3 and n % 4 == 2:
            print("(p % 4 == 3 and n % 4 == 2) => episilon_2 := -1")
            epsilon_2 = -1
        else:
            print("(p % 4 != 3 or n % 4 != 2) => episilon_2 := 1")
            epsilon_2 = 1

        # epsilon and eta
        print("epsilon = epsilon_1 * epsilon_2 = " + str(epsilon_1 * epsilon_2))
        q_p = naik_number_dict[(q, p)]
        d = n / (2 * abs(q_p))
        print("\nnaik_sign = " + str(sign(q_p)))
        print("eta = naik_sign^d = " + str(sign(q_p)^d))
        if sign(q_p)^d == epsilon_1 * epsilon_2:
            print("eta == epsilon\n")
        else:
            print("eta != epsilon\n")


class PrzytyckiTester(knots_periodicity.PrzytyckiTester):

    def __init__(self, K, name, f_homfly_in=None):
        super().__init__(self, K, name, f_homfly_in)
        self.verbose = settings.debugging
        if self.verbose:
            print("\n" + "Knot " + name)
            print("HOMFLYPT = " + str(homflypt))
            print("HOMFLYPT(a, -z) - HOMFLYPT(a^-1, -z) = " +
                  str(homfly_difference))
            print()

    def check_congruence(self, q):
        for i in range(q + 1):
            z_coefficient = self.homfly_difference.coefficient(z^(i+1))
            ideal = (a + a^-1)^(q - i)  # for i == q will be 1
            coefficient_modulo_ideal = z_coefficient.quo_rem(ideal)[1]
            coefficient_modulo_q = coefficient_modulo_ideal.change_ring(GF(q))
            if self.verbose:
                print("\nv_" + str(i) + " = " + str(z_coefficient))
                print("v_" + str(i) + " mod (a + a^-1)^(q - i) = " +
                      str(coefficient_modulo_ideal))
                print("(v_" + str(i) + " mod (a + a^-1)^(q - i)) mod q = " +
                      str(coefficient_modulo_q))
            if coefficient_modulo_q != 0:
                return 0
        return 1


def check_criteria(name, pd_code, A=None, f_homfly_in=None):

    if settings.only_chosen and name not in settings.set_to_check:
        return None

    tester = PeriodicityTester(name, pd_code, None, f_homfly_in)

    for i, q in enumerate(settings.periods):

        if settings.only_periods:
            if tester.name not in settings.periods_dict:
                continue
            if (q not in settings.periods_dict[tester.name] and
                    (-q) not in settings.periods_dict[tester.name]):
                continue

        if settings.only_periods_where_borodzik:
            if tester.name not in settings.fails_dict:
                if tester.name not in settings.success_dict:
                    continue
                if q != settings.success_dict[tester.name]:
                    continue
            else:
                if q != settings.fails_dict[tester.name]:
                    continue

        tester.check_criteria_for_period(q)
        if settings.debugging:
            print("\n" + "#" * 30 + " Calculations for knot " + tester.name +
                   " and q = " + str(q) + " " + "#" * 30 + "\n")
            tester.print_data_for_murasugi(q)
            tester.print_data_for_naik_1(q)
            tester.print_data_for_naik_2(q)
            tester.print_data_for_borodzik(q)

        tester.results.append([q, tester.murasugi, tester.naik_1,
                               tester.naik_2, tester.borodzik,
                               tester.przytycki])
    if settings.print_results:
        tester.print_results()

    return tester


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
            tester.save_results(f_out) #, f_homfly_out)


def check_up_to_10(f_out, f_homfly_out=None, f_homfly_in=None):
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
            tester.save_results(f_out) #, f_homfly_out)


def test_all(f_out, f_homfly_out=None, f_homfly_in=None):
    check_up_to_10(f_out, f_homfly_out, f_homfly_in)
    if f_homfly_out is not None:
        f_homfly_out.flush()
    if f_out is not None:
        f_out.flush()
    check_11_to_15(f_out, f_homfly_out, f_homfly_in)


def main():
    if not os.path.isfile(settings.f_old_results) \
       or not settings.check_old_results:
        settings.check_old_results = False
        if settings.save_homfly and settings.input_file_with_homflypt:
            with open(settings.f_results_out, 'w') as f_out,\
                 open(settings.f_homfly_lm_out, 'w') as f_homfly_out,\
                 open(settings.f_homfly_lm_in, 'r') as f_homfly_in:
                test_all(f_out, f_homfly_out, f_homfly_in)
        elif settings.save_homfly:
            with open(settings.f_results_out, 'w') as f_out,\
                 open(settings.f_homfly_lm_out, 'w') as f_homfly_out:
                test_all(f_out, f_homfly_out)
        elif settings.input_file_with_homflypt:
            with open(settings.f_results_out, 'w') as f_out,\
                 open(settings.f_homfly_lm_in, 'r') as f_homfly_in:
                test_all(f_out, None, f_homfly_in)
        else:
            with open(settings.f_results_out, 'w') as f_out:
                test_all(f_out)


    with open(settings.f_old_results, 'r') as f_old_results:
        if settings.save_homfly and settings.input_file_with_homflypt:
            with open(settings.f_results_out, 'w') as f_out,\
                 open(settings.f_homfly_lm_out, 'w') as f_homfly_out,\
                 open(settings.f_homfly_lm_in, 'r') as f_homfly_in:
                test_all(f_out, f_homfly_out, f_homfly_in)
        elif settings.save_homfly:
            with open(settings.f_results_out, 'w') as f_out,\
                 open(settings.f_homfly_lm_out, 'w') as f_homfly_out:
                test_all(f_out, f_homfly_out)
        elif settings.input_file_with_homflypt:
            with open(settings.f_results_out, 'w') as f_out,\
                 open(settings.f_homfly_lm_in, 'r') as f_homfly_in:
                test_all(f_out, None, f_homfly_in)
        else:
            with open(settings.f_results_out, 'w') as f_out:
                test_all(f_out)


if __name__ == '__main__':

    settings = MySettings()
    S.<a, z> = LaurentPolynomialRing(ZZ)
    R.<t> = LaurentPolynomialRing(ZZ)
    prime_numbers = Primes()
    naik_number_dict =  knots_periodicity.naik_number_dict
    if '__file__' in globals():
        main()
