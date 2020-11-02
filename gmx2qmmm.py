def userInputs():
    from argparse import ArgumentParser
    from argparse import RawDescriptionHelpFormatter
    import textwrap
    parser = ArgumentParser(description="gmx2qmmm, a python interface for Quantum mechanics/Molecular mechanics (QM/MM) calculation.",
                            formatter_class=RawDescriptionHelpFormatter,
                            epilog=textwrap.dedent('''\
                                gmx2qmmm (GNU licenced) by Jan Philipp Goetze.
                                gmx2qmmm is a free QM/MM interface for Gromacs. Compatible with most versions, starting from 5.X.
                                gmx2qmmm uses an additive QM/MM scheme with a numerical charge shift procedure and correction of cut bonds via model potentials.\n
                                ''')
                            )
    
    parser.add_argument("-c", "--coord", help="coordinate(.g96 or .gro)", type=str, default='conf.g96')
    parser.add_argument("-p", "--top", help="Topology (.top)", type=str, default='topol.top')
    parser.add_argument("-n", "--qmatoms", help="QM atoms file(.ndx)", default='qmatoms.ndx')
    parser.add_argument("-qm", "--qmFile", help="QM parameters (.dat)", type=str, default='qm.dat')
    parser.add_argument("-mm", "--mmFile", help="MM parameters (.dat)", type=str, default="mm.dat")
    parser.add_argument("-qmmm", "--qmmmFile", help="QM/MM parameters (.dat)", type=str, default="qmmm.dat")
    parser.add_argument("-act", "--act", help="Active atoms (.ndx)", type=str, default="act.ndx")
    parser.add_argument("-path", "--pathFile", help="path file (.dat)", type=str, default="path.dat")
    parser.add_argument("-g", "--logfile", help=" logfile(.log)", type=str, default="logfile")
    args = parser.parse_args()
    return args

def gmx2qmmm():
    print("gmx2qmmm, a python interface for Quantum mechanics/Molecular mechanics (QM/MM) calculation")



if __name__ == "__main__":
    inputFile = userInputs()
    gmx2qmmm()

