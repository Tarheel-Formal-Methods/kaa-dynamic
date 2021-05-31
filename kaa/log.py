from enum import Enum, auto
from termcolor import cprint, colored

'''
Basic logger utility module for debugging purposes.
'''

#Just basic write to file feature that nicely formats the relevant variable states and outputs them into a readable format.
#Just assume that the debugging info is already known beforehand.
#Quick and dirty logging feature. Not extensive or efficient by any means in its current state.
#
spaces = lambda x: ''.join('  ' for _ in range(x))

class DebugCodes(Enum):
    STEP = auto()
    MINMAX = auto() #For min/max point for parallelotopes
    POLY = auto() #For min/max polynomals calculated for bernstein expansion
    LOCAL_BOUND = auto() #Local bounds for each parallelotope for offl, offu
    GLOBAL_BOUND = auto() #Final bound after considering all parallelotopes
    A_B = auto()
    PROJ_MINMAX = auto()

_Debug_Strings = {
DebugCodes.MINMAX: 'Min/Max points for Parall {0}: {1}    {2}',
DebugCodes.POLY:  'UPoly: {0}   LPoly: {1}',
DebugCodes.LOCAL_BOUND: 'MaxB: {0}   MinB: {1}    for P: {2}',
DebugCodes.GLOBAL_BOUND: 'New Offu: {0}    New Offl: {1}',
DebugCodes.STEP: 'STEP: {0} \n ----------------------------------------------',
DebugCodes.A_B: 'A:  {0}  \n\n  b: {1}',
DebugCodes.PROJ_MINMAX: 'STEP {0}: \n  MIN: {1}   MAX: {2}'
}

def write_log(*args):
    debug_cat = args[-1]
    with open('log.txt','a') as f:
        f.write('\n' + spaces(debug_cat.value) + _Debug_Strings[debug_cat].format(*args) + '\n')


"""
Logging class should have the following features
- Dumping Directons/Ptope matrix into readable format
- Dumping Strategy moves (Adding/Removing templates) into readable format
- Which points were chosen during PCA/LinApp Strategy.
- Which directions were computed as output from PCA/LinApp.
"""

class Output:

    @staticmethod
    def warning(output):
        cprint(''.join((output.upper(), '\n')), 'red', attrs=['blink'])

    @staticmethod
    def prominent(output):
        cprint(''.join((output.upper(), '\n')), 'white', attrs=['bold'])

    @staticmethod
    def write(output):
        with open('log3.txt','a') as f:
            f.write(output + '\n')

    @staticmethod
    def bold_write(output):
            with open('log2.txt','a') as f:
                f.write(colored(output, 'white', attrs=['bold']) + '\n')
