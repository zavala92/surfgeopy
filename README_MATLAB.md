# How to Call surfgeopy/minterpy/ minterpy_levelsets/ in MATLAB

## 1. Configure Your System to Use Python

- [Configure MATLAB to Use Python](https://www.mathworks.com/help/matlab/matlab_external/install-supported-python-implementation.html)

- To determine which version MATLAB is using, call `pyenv`. For example:

    ```matlab
    pe = pyenv
    ```

    Output:
    ```
    PythonEnvironment with properties:

        Version: "3.9"
        Executable: "/usr/bin/python3"
        Library: "libpython3.9.so.1.0"
        Home: "/usr"
        Status: NotLoaded
        ExecutionMode: InProcess
    ```

- Specify Python executable:
  
    ```matlab
     py_exe =  "/usr/bin/python3";
    ```
    
    In case you would like to install the `matplotlib` package:
  
    ```matlab
     system([py_exe, ' -m pip install matplotlib']);
    ```

## 2. 🚀 Installing from a Cloned Repository

If you prefer installing from a cloned repository, such as [minterpy_levelsets](https://github.com/minterpy-project/minterpy-levelsets), here's how you can do it:

```matlab
    
    % path_minterpy_levelsets = 'path where minterpy-levelsets is located';
    system([py_exe, ' -m pip install ',  path_minterpy_levelsets]);
    
```

In case you would like to uninstall minterpy_levelsets package:

 ```matlab
    system([py_exe, ' -m pip uninstall -y minterpy_levelsets']);
```

## 3. Calling Python Codes from MATLAB

After you have installed all the dependencies, you may start calling Python codes from MATLAB.

### Example

Assume we have the following file: `level_ellipsoid.py`

```python
import numpy as np
import sympy as sp
import minterpy as mp
import minterpy_levelsets as ls

a = 1
b = 2
c = 3

x, y, z = sp.symbols('x y z')
expr = (x**2/a**2) + (y**2/b**2) + (z**2/c**2) - 1
poly = sp.Poly(expr, x, y, z)
newt_poly_exact = ls.sympy_to_mp(poly, mp.NewtonPolynomial)
point_data = ls.sample_points(newt_poly_exact,  200,   bounds=4.0, tol=1e-15, random_seed=42) # random seed
ls.output_VTK(point_data)
ls.output_VTR(newt_poly_exact, bounds=1.0)
poly = ls.LevelsetPoly(point_data, method='BK', tol=1e-11, verbose=True)
distance_errors = poly(point_data) / np.linalg.norm(poly.compute_gradients_at(point_data), axis=1)
linf_error = np.max(np.abs(distance_errors))
print(f"L_inf error = {linf_error}")
validation_point_data = ls.sample_points(newt_poly_exact,  # Polynomial in Newton basis
                          100,        # Number of points to be sampled
                          bounds=4.0, # Boundary of the Cubic domain to be sampled
                          tol=1e-15,  # Tolerance in solution
                          random_seed=1729) # random seed
validation_distance_errors = poly(validation_point_data) / np.linalg.norm(poly.compute_gradients_at(validation_point_data),axis=1)
linf_validation_error = np.max(np.abs(validation_distance_errors))
print(f"L_inf validation error = {linf_validation_error:.3}")
vtk_data = ls.output_VTK(point_data);
linf_validation_error = linf_validation_error.item();

```
To run the level_ellipsoid.py in MATLAB command window, you need to write:


```matlab
    level_ellipsoid_module = py.importlib.import_module('level_ellipsoid');
    %py.importlib.reload(level_ellipsoid_module);
 ```
The output after running would look like below:

```matlab
    Levelset error (method = 'BK'), n = 2, lp = 2.0 : 0.011339074022902632
    Levelset error (method = 'BK'), n = 3, lp = 2.0 : 0.00015041454565755763
    Levelset error (method = 'BK'), n = 4, lp = 2.0 : 8.326672684688674e-17
    L_inf error = 6.478176983772814e-12
    L_inf validation error = 4.19e-12
```
