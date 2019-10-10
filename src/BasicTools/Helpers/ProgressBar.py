# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
                       
"""
you can filter the output to "reinterpret" the control chars using :

    col -b  < log  > log2

    where log is the original logfile

"""
def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '¦'):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end = '\r')
    import sys
    sys.stdout.flush()
    # Print New Line on Complete
    if iteration == total:
        print()


def CheckIntegrity():

  from time import sleep

  # A List of Items
  items = list(range(0, 57))
  l = len(items)

  # Initial call to print 0% progress
  printProgressBar(0, l, prefix = 'Progress:', suffix = 'Complete', length = 50)
  for i, item in enumerate(items):
    # Do stuff...
    sleep(0.1)
    # Update Progress Bar
    printProgressBar(i + 1, l, prefix = 'Progress:', suffix = 'Complete', length = 50)
  return 'ok'

if __name__ == '__main__':
    CheckIntegrity() # pragma: no cover

