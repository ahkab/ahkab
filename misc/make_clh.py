import os
os.system("ahkab -h > command")

fin = open('command', 'r')
fp = open('Command-Line-Help.rst', 'w')

fp.write(
"""Command line help
=================

The ``ahkab`` simulator has a command line interface that
allows for quick simulation of netlist decks, without the need
to load the Python interpreter explicitely.

Several switches are available, to set the input and output files
and to override some built-in options.

Notice that options set on the command line always take precedence
on any netlist option or any value set in :module:`ahkab.options`.

""")

indent = False
for line in fin:
    line = line.strip()
    if not len(line):
        fp.write('\n')
        indent = False
        continue
    elements = line.split(' ')
    if len(elements) and elements[0][-1] == ':' and elements[0] != 'Default:':
        fp.write(elements[0] + '\n' + '-'*len(elements[0]) +'\n'*2)
        elements = elements[1:]
    ew = []
    suppress = False
    for i in range(len(elements)):
        if elements[i] == '' and ((len(elements) > i + 1
                                   and elements[i+1] == '')
                                  or len(elements) == i+1):
            suppress = True
        else:
            if suppress == True:
                ew.append('\n   ')
                indent = True
                suppress = False
            elif elements[i][0] == '\t' and i == 0:
                ew.append('::\n\n')
                ew.append('    ' + elements[i][1:])
            else:
                ew.append(elements[i])
    if len(ew):
        if ew[0][0] != '-' and indent:
            fp.write('    ')
        fp.write(' '.join(ew) + '\n')
        if ew[0][0] == '-':
            indent == True
fp.close()
fin.close()
os.remove('command')
