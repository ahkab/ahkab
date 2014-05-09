# Package testing can be done remotely, without display. This would make
# matplotlib fail (and, sometimes, the test as well).
# If we check for $DISPLAY, that makes us probably lose in portability,
# because does Windows have the DISPLAY env variable defined?
import os
import matplotlib
if not os.system('python -c "import matplotlib.pyplot as plt;plt.figure()"'):
    matplotlib.use('Agg')

from . import *
