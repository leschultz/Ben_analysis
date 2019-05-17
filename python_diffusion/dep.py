'''
Find pertinent information in LAMMPS input files.
'''


def info(item):
    '''
    Parse the input file to grab run information.
    '''

    # Open file and iterate per line
    runsteps = []
    with open(item) as file:
        for line in file:
            value = line.strip().split(' ')

            # Find timestep if word in line
            if 'timestep' in value[0]:
                for i in value:
                    try:
                        timestep = float(i)
                    except Exception:
                        pass

                # Find all defined run steps and save under runsteps
            if 'run' in value:
                for i in value:
                    try:
                        runsteps.append(int(i))
                    except Exception:
                        pass

            # Determine the elements in the system
            if 'pair_coeff' in value:
                types = [i for i in value if i != ''][4:]

                elements = []
                for element in types:
                    elements.append(element)

            if 'dump' in value:
                trajdumprate = int(value[-7])

        increment = sum(runsteps[-2:])  # The steps to the start of hold

        hold1 = sum(runsteps[:-2])  # The steps until temperature steps
        hold2 = runsteps[-2]  # The start of hold
        hold3 = runsteps[-1]  # The end of hold

        # Save paramters for each trajectory file
        parameters = {
                      'timestep': timestep,
                      'runsteps': runsteps,
                      'trajdumprate': trajdumprate,
                      'hold1': hold1,
                      'hold2': hold2,
                      'hold3': hold3,
                      'increment': increment,
                      'elements': elements,
                      }

    return parameters
