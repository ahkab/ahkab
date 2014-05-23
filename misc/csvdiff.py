from __future__ import print_function
import sys
from ahkab import csvlib

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("CSV subtractor. Usage:\n\tcsvdiff file1.csv file2.csv [ouput.csv]\n" +
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
    csvlib.write_csv(outfile, a-b, h)
