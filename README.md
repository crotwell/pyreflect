# pyreflect
Python codes to work with George Randall's synthetic seismogram code

# Example

```
from pyreflect.earthmodel import EarthModel
mod = EarthModel.loadFromFile('testmgen.ger')
print(mod.asJSON())
print(mod.asGER())
gradmod = mod.gradient(1, .05, .05, 5)
print(gradmod.asGER())
eftmod = gradmod.eft(5)
print(eftmod.asGER())
```
