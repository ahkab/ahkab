# -*- coding: iso-8859-1 -*-
# ticker.py
# Progress indicator class
# Copyright 2006 Giuseppe Venturini

# This file is part of the ahkab simulator.
#
# Ahkab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 2 of the License.
#
# Ahkab is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License v2
# along with ahkab.  If not, see <http://www.gnu.org/licenses/>.

"""
A progress indicator.
"""

from __future__ import (unicode_literals, absolute_import,
                        division, print_function)

import sys

class ticker:

    """This is a progress indicator class.

    If activated, you shouldn't print anything to screen before calling
    ticker.hide().

    If you wish to change the progress indicator, change self.progress to
    something else.
    """
    def __init__(self, increments_for_step=10):
        self.progress = ("-", "\\", "|", "/")
        self._index = 0
        self._step = 0
        self._display = False
        self.increments_for_step = increments_for_step

    def step(self):
        """After calling this function ticker.increments_for_step times
        the status is incremented."""
        if not self._display:
            return
        if (self._index + 1) % self.increments_for_step == 0:
            if (self._step + 1) % len(self.progress) == 0:
                self._step = 0
            else:
                self._step = self._step + 1
            self._index = 0
            if self._display:
                sys.stdout.write("\b" + self.progress[self._step])
                sys.stdout.flush()
        else:
            self._index = self._index + 1

    def hide(self, enable=None):
        """Before printing text to screen, call this to hide the progress
        indicator.
        """
        if enable == False:
            return
        sys.stdout.write("\b")
        sys.stdout.flush()
        self._display = False

    def display(self, enable=None):
        """Print to screen the progress indicator. Call hide to hide it
        again.
        """
        if enable == False:
            return
        sys.stdout.write(self.progress[self._step])
        sys.stdout.flush()
        self._display = True

    def reset(self):
        """Reset to initial status. Doesn't hide it."""
        self._step = 0
        self._index = 0

if __name__ == "__main__":
    # test
    import time
    tk = ticker()
    tk.display(True)
    for i in range(1000):
        time.sleep(0.005)
        tk.step()
    tk.hide()
