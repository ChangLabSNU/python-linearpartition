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

The module currently only has one function called `partition(seq)`, and
it doesn't have any customizable options other than the defaults.
The seq parameter should be an RNA sequence in *uppercase* letters,
and any `T` should be converted to `U` before passing it to the function.

```python
>>> import linearpartition as lp
>>> seq = 'UGUCGGGGUUGGCUGUCUGACA'
>>> bpmtx, fe = lp.partition(seq)
>>> fe
-7.216465644007023
>>> import pandas as pd
>>> pd.DataFrame(bpmtx).sort_values('prob', ascending=False).head()
    i   j      prob
19  3  18  0.999201
18  2  19  0.998801
17  1  20  0.997717
21  5  16  0.996692
22  4  17  0.996508
```

The function returns a tuple with two elements. The first element is the
predicted base pairing probability, and the second element is the free
energy of the ensemble structure in kcal/mol.

### Author

Hyeshik Chang &lt;hyeshik@snu.ac.kr&gt;

### License

This Python binding is licensed under [the MIT-style license](LICENSE).
However, the compiled binary includes code from the LinearPartition
package, which is licensed for non-commercial use.
