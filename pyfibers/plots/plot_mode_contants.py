from pyfibers.fibers import LeKienFiber
from pyfibers.modes.lekien import LeKienRadMode
from functools import partial
import numpy as np
import matplotlib
matplotlib.use('TKAgg')
from matplotlib import pyplot as plt

attrs = set()
value_dicts = []

Vs = np.linspace(0.1, 5, 101)

exclude = ('discrete_modes', 'fiber')

for V in Vs:
    fiber = LeKienFiber(1.45, 1, 1.0, V)
    mode = LeKienRadMode(fiber, 1, V*1.1, 1, 1)
    value_dict = {}
    for attr_name in dir(mode):
        if attr_name.startswith('_') or attr_name in exclude:
            continue
        attr = getattr(mode, attr_name)
        #if hasattr(attr, 'im_self'):
        #    attr = partial(attr, attr.im_self)
        if not callable(attr):
            val = attr
        else:
            try:
                val = attr()
            except (TypeError, ValueError, AttributeError), e:
                try:
                    val = attr(1)
                except (TypeError, ValueError, AttributeError), e:
                    continue
        attrs.add(attr_name)
        value_dict[attr_name] = val
    value_dicts.append(value_dict)

print attrs
print value_dicts

for i, attr_name in enumerate(attrs):
    values = [d[attr_name] for d in value_dicts]
    plt.figure(i)
    plt.title(attr_name)
    plt.plot(Vs, values)
    plt.show()

