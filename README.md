# pyreflect
Python codes to work with George Randall's synthetic seismogram code

Some code copied from my old TclTk code [TkReflect](https://github.com/crotwell/TkReflect).

# Example

```
from pyreflect.earthmodel import EarthModel
mod = EarthModel.loadFromFile('testmgen.ger')
print(mod.asJSON())
print(mod.asGER())
gradmod = mod.gradient(1, .01, .01, 5)
print(gradmod.asGER())
eftmod = gradmod.eft(5)
print(eftmod.asGER())
```

# Build

See https://packaging.python.org/tutorials/packaging-projects/#creating-the-package-files

```
python3 -m pip install --upgrade build
python3 -m build
pip install dist/pyreflect-0.0.3-py3-none-any.whl --force-reinstall
```
