
import numpy as np

a = np.array([5., np.nan])
print(np.where(~np.isnan(a))[0].size)
