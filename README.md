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
