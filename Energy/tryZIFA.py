import numpy as np
from ZIFA import ZIFA
from ZIFA import block_ZIFA
import pandas as pd

# This gives an example for how to read in a real data called input.table. 
# genes are columns, samples are rows, each number is separated by a space. 
# If you do not want to install pandas, you can also use np.loadtxt: https://docs.scipy.org/doc/numpy/reference/generated/numpy.loadtxt.html

file = pd.read_csv('/home/patrick/gdrive/CAIL/Energy/df3.csv', index_col = 0) 
table = np.array(file)
Z, model_params = ZIFA.fitModel(table, 2)

pd.write

np.savetxt('/home/patrick/gdrive/CAIL/Energy/ZIFAScores.tab', Z, fmt='%.2f')
np.savetxt('/home/patrick/gdrive/CAIL/Energy/ZIFALoad.tab', model_params['A'], fmt='%.2f')