from apexpy import Apex
from datetime import datetime

utcs = [[2, 10, 9, 51, 27],
       [2, 12, 10, 22, 11],
       [3, 4, 7, 30, 12],
       [3, 4, 10, 13, 49],
       [3, 14, 6, 49, 7],
       [3, 19, 8, 23, 30]]

locs = [[65.5, 216],
        [66.0, 212],
        [63.5, 216],
        [66.5, 208],
        [64.0, 217],
        [67.0, 215]]

for i in range(6):
    utc = utcs[i]
    loc = locs[i]
    dt = datetime(2023, utc[0], utc[1], utc[2], utc[3], utc[4])
    print(dt)
    A = Apex(dt)
    mlat, mlt = A.convert(loc[0], loc[1], 'geo', 'mlt', height=110, datetime=dt)
    print('{:.1f}'.format(mlt))
