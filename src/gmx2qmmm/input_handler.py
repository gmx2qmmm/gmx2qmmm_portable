#   // INITIAL DESCRIPTION //
"""Short Module Description; Reference To Readme"""

#   // MEATDATA //
__author__ = 'Florian Anders'
__date__ = '2024-01-09'

#   // IMPORTS //

#   Imports Of Existing Libraries
import os
import sys

#   Imports From Existing Libraries
from collections import defaultdict
from typing import Optional

#   Imports Of Custom Libraries

#   Imports From Custom Libraries
from gmx2qmmm.logging import Logger
from gmx2qmmm.types import StrPath

#   // TODOS & NOTES //
#   TODO: This Is An Example To Do
#   NOTE: This Is An Example Note

#   // CLASS & METHOD DEFINITIONS //
class FileReader():

    '''
    This Class Reads Specified Files And Provides Methods To Extract Specific Information Written In These Files
    '''

    @staticmethod
    def read_file(file: StrPath, logfile: Optional[StrPath]) -> defaultdict:

        '''
        ------------------------------ \\
        EFFECT: \\
        --------------- \\
        NONE \\
        ------------------------------ \\
        INPUT: \\
        --------------- \\
        NONE \\
        ------------------------------ \\
        RETURN: \\
        --------------- \\
        NONE \\
        ------------------------------ \\
        '''

        try:

            with open(file) as fp:
                list_content_file: list[str] = fp.readlines()

        except FileNotFoundError:

            Logger.log(
                logfile,
                (
                    'An error occured when trying to read the parameter file - \n'
                    'Make sure that the file exists and try running gmx2qmmm again.'
                )
            )

            print(
                '\n\nAn error occured when trying to read the parameter file - \n'
                'Make sure that the file exists and try again, running the command\n\n'
                '\"python gmx2qmmm.py -p <ParameterFile>\"\n\n'
            )

            sys.exit(1)

        except ValueError:

            Logger.log(
                logfile,
                (
                    'An error occured when trying to read the parameter file - \n'
                    'We encountered a value mismatch. Please check the correctness of your input types!'
                )
            )
            sys.exit(1)


        #   Initialise An Empty Defaultdict;
        #   This Will Hold Key / Value Pairs Based On The Inputs
        defaultdict_parameters_input: defaultdict = defaultdict()

        for str_pair_parameter in list_content_file:

            #   Ignore Newline Characters ('\n');
            #   Ignore Empty Lines;
            #   Ignore Lines Beginning With A Hashtag Symbol (#) From The List Of Parameters;
            #   For Now, Ignore Lines Beginning With Exclamation Mark (!); This Is Used For Include Statements, Which We Handle Later
            if str_pair_parameter == '\n' or str_pair_parameter == '' or str_pair_parameter.startswith('#') or str_pair_parameter.startswith('!'):

                continue

            # XX Florian not adding empty parameters to the list
            elif str_pair_parameter.strip().split('=')[1] == '':

                continue

            else:

                #   Populate The Defaultdict For Each Line In The Parameter Input File With:
                #   str_pair_parameter.strip().split('=')[0] -> The First Part Of The Line Text, Before The Equality Sign ('=') - This Becomes The Key Of The Defaultdict
                #   str_pair_parameter.strip().split('=')[1] -> The Second Part Of The Line Text, After The Equality Sign ('=') - This Becomes The Value Of The Defaultdict
                #   See The Syntax Of Constructing A Defaultdict For Further Information On The Following Line Of Code
                defaultdict_parameters_input[str_pair_parameter.strip().split('=')[0]] = str_pair_parameter.strip().split('=')[1]

        #   Check For Usage Of 'Library' Files After We Have Populated The Initial Defaultdict;
        #   These Should be Included In The Parameter File With '!INCLUDE';
        #   Use The Specified Parameters In The Include File And Override The Previously Created Defaultdict Values
        for str_pair_parameter in list_content_file:

            if str_pair_parameter.startswith('!INCLUDE'):

                #   Check, If A Value Was Set
                if str_pair_parameter.strip().split('=')[1] != '':

                    #   If Value Was Set, It Has To Be The File Path + Name (Absolute Path) Of The File To Be Included
                    file_parameters_include = str_pair_parameter.strip().split('=')[1]

                    with open(os.path.abspath(file_parameters_include)) as fp:

                        list_content_file_include: list = fp.readlines()

                    #   Iterate Over The List Of Superior Parameters
                    for str_pair_parameter_superior in list_content_file_include:

                        #   XX AJ again checking for empty line etc., otherwise we're getting an index error later
                        #   Ignore Newline Characters ('\n');
                        #   Ignore Empty Lines;
                        #   Ignore Lines Beginning With A Hashtag Symbol (#) From The List Of Parameters;
                        #   For Now, Ignore Lines Beginning With Exclamation Mark (!); This Is Used For Include Statements, Which We Handle Later
                        if str_pair_parameter_superior == '\n' or str_pair_parameter_superior == '' or str_pair_parameter_superior.startswith('#'):

                            continue

                        elif str_pair_parameter_superior.strip().split('=')[1] == '':

                            continue

                        else:
                            #   Check, If For Each parameter A Value Has Been Set
                            #   XX AJ added strip to remove '\n', otherwise '\n' as value in dict
                            if str_pair_parameter_superior.strip().split('=')[1] != '':

                                #   Use Parameter Only If Not Defined Before Or With Priority Parameter

                                if int(defaultdict_parameters_input['includeprio']) or str_pair_parameter_superior.strip().split('=')[0] not in defaultdict_parameters_input:

                                    #   Overwrite The Respective Parameter In The Defaultdict With The Value To be Used
                                    defaultdict_parameters_input[str_pair_parameter_superior.strip().split('=')[0]] = str_pair_parameter_superior.strip().split('=')[1]

        # clean the dictionary
        for key, value in defaultdict_parameters_input.items():
            if value.startswith('[') and value.endswith(']'):
                # If the value is enclosed in square brackets, it's already a list
                # Remove the brackets and split the contents into a list
                value = value[1:-1]
                if value:
                    value = value.split(',')
                else:
                    defaultdict_parameters_input[key] = []
            else:
                # Try converting to integer
                try:
                    defaultdict_parameters_input[key] = int(value)
                except ValueError:
                    # If conversion fails, keep the value as is
                    pass

        #   Check, If The Inputs In The Defaultdict Are:
        #   1) Complete For All Required Inputs,
        #   2) All Required Inputs Have Values Assigned,
        #   3) All Numerical Values Are, In Fact, Numerical (int Should Be int, float Should Be float, etc.)
        #   TODO: Implement Assessment!

        Logger.log(logfile, 'Asserting basic correctness of inputs')

        try:

            Asserter.assert_defaultdict_input(defaultdict_parameters_input)

            Logger.log(logfile, 'Basic correctness check COMPLETE')

            return defaultdict_parameters_input

        except AssertionError:

            Logger.log(logfile,
                (
                    'Basic correctness check FAILED;\n'
                    'Please make sure that all necessary parameters are provided in your parameter file!'
                )
            )

            print('An error occured! Check the logfile (\'{0}\') for more information.'.format(logfile))

            sys.exit(1)


class Asserter:

    '''
    This Class Assesses The Correctness Of Inputs (Currently Only On A Basic Level)
    '''

    def __init__(self) -> None:

        '''
        ------------------------------ \\
        EFFECT: \\
        --------------- \\
        NONE \\
        ------------------------------ \\
        INPUT: \\
        --------------- \\
        NONE \\
        ------------------------------ \\
        RETURN: \\
        --------------- \\
        NONE \\
        ------------------------------ \\
        '''

        pass


    @staticmethod
    def assert_defaultdict_input(defaultdict_to_assert: defaultdict) -> None:

        '''
        ------------------------------ \\
        EFFECT: \\
        --------------- \\
        Asserts The Basic Correctness Of Inputs Stored in A Defaultdict \\
        ------------------------------ \\
        INPUT: \\
        --------------- \\
        defaultdict_to_assert: defaultdict -> Defaultdict Holding Key / Value Pairs Based On Inputs \\
        ------------------------------ \\
        RETURN: \\
        --------------- \\
        NONE \\
        ------------------------------ \\
        '''

        #   TODO: Expand Assertion Functionality!
        bool_parameters_required_existing = False

        list_parameters_input_required =\
            [
                'gaussianpath',
                'gromacspath',
            ]
        #gaussianexepath=
        #gromacsexepath=
        #gromacscoordinatespath=
        #gromacstopologypath=
        #jobtype=nma
        ###pccores=3
        ###pcmemory=
        #qmprogram=gaussian
        ###qmbasisset=
        ###qmmethod=
        #qmatomslist=[1, 2, 3, 4]
        ###systemcharge=
        ###systemmultiplicity=
        #activeatomslist=[1, 2, 3, 4, 5, 6, 7, 8, 9]
        #useinnerouter=False
        ##inneratomslist=[]
        ##outeratomslist=[]
        #   TODO: Assertion HERE
        bool_parameters_required_existing = True

        if bool_parameters_required_existing:

            pass

        else:

            raise AssertionError

if __name__ == '__main__':

    #   TODO: Unit Testing Should Follow Below!

    pass

