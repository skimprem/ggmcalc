import numpy as np
from mycode import square_array as sa # импорт модуля, созданного f2py

x = np.array([1.0, 2.0, 3.0], dtype=np.float64)
x = np.asfortranarray(x)  # для совместимости по памяти (не всегда обязательно, но безопасно)

a = 2.3

text = 'test'

conf = {
    'ii': 10,
    'rr': 33.4,
    'ss': 'text'
}

sa(x, a, text)
