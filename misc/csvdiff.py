from __future__ import print_function
import sys
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline
from ahkab import csvlib

if __name__ == '__main__':
    if '-i' in sys.argv:
        sys.argv.remove('-i')
        interpolate = True
    else:
        interpolate = False
    if len(sys.argv) < 3:
        print("CSV subtractor. Usage:\n\tcsvdiff [-i] file1.csv file2.csv [ouput.csv]\n" +
              "Saves to the optional output file, or to stdout, file1.csv - file2.csv.")
        sys.exit(1)
    a, h1, _, _ = csvlib.load_csv(sys.argv[1])
    b, h2, _, _ = csvlib.load_csv(sys.argv[2])

    assert len(h1) == len(h2)
    if h1 == h2:
        h = h1
    else:
        h = [a + "-" + b for a, b in zip(h1, h2)]

    outfile = sys.argv[3] if len(sys.argv) > 3 else sys.stdout

    if interpolate and a.shape[0] > 1:
        data = None
        for j in range(a.shape[0]):
            a_int = InterpolatedUnivariateSpline(a[0, :], a[j, :])
            b_int = InterpolatedUnivariateSpline(b[0, :], b[j, :])
            if data is None:
                data = (a_int(a[0, :]) - b_int(a[0, :])).reshape((1, -1))
            else:
                data = np.vstack((data, (a_int(a[0, :]) - b_int(a[0, :]).reshape((1, -1)))))
    else:
        shortest_len = min(a.shape[1], b.shape[1]) - 1
        a = a[:, :shortest_len]
        b = b[:, :shortest_len]
        data = a - b

    csvlib.write_csv(outfile, data, h)
