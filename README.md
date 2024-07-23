# python-linearpartition

Unofficial CPython binding to LinearPartition

### Installation

Use `pip` to install the module.

```bash
pip install linearpartition-unofficial
```

You may build from the source code for unsupported Python versions or platforms.

```bas
git clone --recursive https://github.com/ChangLabSNU/python-linearpartition
cd python-linearpartition
pip install .
```

### Usage

The module currently only has one function called `partition(seq)`.
The seq parameter should be an RNA sequence in *uppercase* letters,
and any `T` should be converted to `U` before passing it to the function.

```python
>>> import linearpartition as lp
>>> seq = 'UGUCGGGGUUGGCUGUCUGACA'
>>> pred = lp.partition(seq)
>>> pred['free_energy']
-7.216465644007023
>>> pred['structure']
'(((((((........)))))))'
>>> import pandas as pd
>>> pd.DataFrame(pred['bpp']).sort_values('prob', ascending=False).head()
    i   j      prob
19  3  18  0.999201
18  2  19  0.998801
17  1  20  0.997717
21  5  16  0.996692
22  4  17  0.996508
```

### Functions

#### linearpartition.partition()

The `linearpartition.partition` function is a Python C extension function that
calls [LinearPartition](https://github.com/LinearFold/LinearPartition) to
perform a linear partitioning operation and get the base pairing probability
matrix.

```python
linearpartition.partition(seq, mode='eterna', beamsize=100, dangles=2)
```

##### Parameters

- `seq` (required): A string containing the RNA sequence to be analyzed.
  The sequence must be in uppercase and only contain A, C, G, and U.
  This parameter is required.
- `mode` (optional): The name of free energy parameters to use. Use
  `'vienna'` for Vienna RNA parameters, or `'eterna'` for EternaFold
  parameters.
- `beamsize` (optional): An integer representing the beam size for the
  operation. Larger value requires more computational time and memory.
  The default value is 100.
- `dangles` (optional): An integer representing the number of dangles for
  the partitioning operation. The default value is 2.

##### Return Value

This function returns a dictionary containing the MEA structure,
base-pairing probability matrix and free energy of the ensemble
structure in kcal/mol from the result of the partitioning operation.

### Author

Hyeshik Chang &lt;hyeshik@snu.ac.kr&gt;

### License

This Python binding is licensed under [the MIT-style license](LICENSE).
However, the compiled binary includes code from the LinearPartition
package, which is licensed for non-commercial use.
