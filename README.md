# A Sage script to verify periodicity

This repository contains a SageMath script used for computations in [Finite orthogonal groups and periodicity of links](https://arxiv.org/abs/1810.03881v1).



-----


The following criteria for knots periodicity were implemented:

- Przytycki's criterion for HOMFLYPT polynomials,
- Murasugi's criterion,
- Naik's homological criterion,
- criterion introduced in Theorem 1.4. in above mentioned paper (called in the script and the following description *Borodzik's criterion*).



-----



For a given list of knots [encoded in PD-code](http://katlas.org/wiki/Planar_Diagrams) (an input file is required) the script generates output with results for each criterion.


-----


We tested all knots up to 15 crossings and for periods 3, 5, 7, 9, 11, 13.

## Running the script


### Prerequisites

- [SageMath 8.1](http://www.sagemath.org).

- For Przytycki's criterion the module [libhomfly](https://github.com/miguelmarco/libhomfly/blob/master/README.md) is needed to calculate HOMFLYPT polynomial. There is also possibility to read saved polynomials from a file (file `homflypt.input` in the repository contains HOMFLYPTY polynomials for some knots).



### Input data

Files with PD-codes for all knots up to 15 crossing are available in the repository (`knots_11_15.txt` and `knots_3_10.txt`) and should be saved in the same folder as the script or a path in the `class MySettings` should be adjusted. Similarly, if other input files will be used, settings have to be changed. A correction in parsing or reading functions may be also necessary (`check_11_to_15()`, `check_up_to_10()`, `parse_pd_code()`, `parse_knot_name()`).


### Running in terminal

The script should be run in terminal: `sage periodicity.sage`.

### Results

All results will be saved in a file `results.out`. A line, for example `12a100,3,1,1,1,0,1`, should be read as follow. A knot 12a100 in tests for a period equal 3 (second position) gives results saved in five numbers: `1,1,1,0,1`, which correspond to criteria: Murasugi's, Naik's - part 1, Naik's - part 2, Borodzik's, Przytycki's. Only `0` means that results excluded periodicity. As computations of first 4 numbers depend on previous results, ones zero occurs, the following criteria will be skipped. The Przytycki's criterion is calculated independently of the others. `-1` means that criterion doesn't exclude periodicity, but the criterion wasn't applicable or that the testing algorithm met an empty list of arguments to be verified.




## Authors

* **Maria Marchwicka**
* **Wojciech Politarczyk**

## License

This project is licensed under the GNU General Public License - see the [LICENSE.md](LICENSE.md) file for details.
