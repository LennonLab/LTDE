from __future__ import division

import numpy as np
import gillespy2 as gp
import matplotlib.pyplot as plt

from gillespy2.solvers.numpy.basic_ode_solver import BasicODESolver

class Simple1(gp.Model):
    """
    This is a simple example for mass-action degradation of species S.
    """

    def __init__(self, parameter_values=None,volume=1.0):

        # Initialize the model.
        gp.Model.__init__(self, name="survivorship")



        # Parameters
        # maximum uptake rate
        v = gp.Parameter(name='v', expression=220) # fg C per cell per day
        # half-saturation constant
        K = gp.Parameter(name='K', expression=1.4e8) # fg C L-1
        # carbon content of cell
        c = gp.Parameter(name='c', expression=100) # fg C cell-1
        # maintenance cost as fraction of cellular carbon
        m = gp.Parameter(name='m', expression=5)
        # fraction of dead cells that are recycled to substrate
        r = gp.Parameter(name='r', expression=0.005)
        # death rate of viable cells
        d = gp.Parameter(name='d', expression=0.03)


        self.add_parameter(v)
        self.add_parameter(K)
        self.add_parameter(c)
        self.add_parameter(m)
        self.add_parameter(r)
        self.add_parameter(d)

        # Species
        N = gp.Species(name='N', initial_value=int(1e9))
        self.add_species(N)
        #D = gp.Species(name='D', initial_value=int(0))
        #self.add_species(D)
        S = gp.Species(name='S', initial_value=int(0))
        self.add_species(S)


        #b_N = Reaction(name = "r1", reactants = {}, products = {N:1},
        #          propensity_function = "(N*v * (S/(S+K))) / (c*m)")

        #b_S = Reaction(name = "r1", reactants = {}, products = {S:1},
        #          propensity_function = "(N*v * (S/(S+K))) / (c*m)")


        # Reactions
        # N + S -> N
        N_S_to_N = gp.Reaction(
                name = 'b_N',
                reactants = {N:1, S:1},
                products = {N:1},
                propensity_function = '(N*v * (S/(S+K))) / (c*m)' )

            

        N_to_S = gp.Reaction(
                name = 'd_N',
                reactants = {N:1},
                products = {S:1},
                propensity_function = 'r*d*N' )

        self.add_reaction(N_S_to_N)
        self.add_reaction(N_to_S)

        #dDdt = gp.Reaction(
        #        name = 'dDdt',
        #        reactants = {N:1},
        #        products = {D:1},
        #        propensity_function = 'd*N - r*D' )

        #self.add_reaction(dDdt)

        #dSdt = gp.Reaction(
        #        name = 'dSdt',
        #        reactants = {S:1,D:1},
        #        products = {S:1},
        #        propensity_function = 'c*r*N - N*v*(S/(S+K))' )

        #self.add_reaction(dSdt)

        self.timespan(np.linspace(0, 1000, 200))



if __name__ == '__main__':

    num_trajectories = 1
    model = Simple1()

    simple_1trajectories = model.run(number_of_trajectories = num_trajectories, show_labels=False)

    time = np.array(simple_1trajectories[0][:,0])

    #print(time)

    print(simple_1trajectories)

    #print(simple_1trajectories[0][:,1])
    #print(simple_1trajectories[0][:,2])
    #print(time)

    #d_results = model.run(solver = BasicODESolver, show_labels = False)


    #plt.plot(time, d_results[0][:,1], '-r', label='U')
    #plt.plot(d_results[0][:,0], d_results[0][:,2], '-b', label='V')
    #plt.plot([0], [11])
