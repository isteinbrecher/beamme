# -*- coding: utf-8 -*-
# -----------------------------------------------------------------------------
# MeshPy: A beam finite element input generator
#
# MIT License
#
# Copyright (c) 2021 Ivo Steinbrecher
#                    Institute for Mathematics and Computer-Based Simulation
#                    Universitaet der Bundeswehr Muenchen
#                    https://www.unibw.de/imcs-en
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
# -----------------------------------------------------------------------------
"""
This script is used to test the functionality of the simulation manager.
"""

# Python imports.
import os
import unittest
import shutil
import numpy as np
import subprocess

# User imports.
from meshpy import (mpy, InputFile, MaterialReissner, BoundaryCondition,
    Beam3rHerm2Line3, Function, set_header_static, set_runtime_output)
from meshpy.mesh_creation_functions import create_beam_mesh_line
from meshpy.simulation_manager import (Simulation, SimulationManager)
from tests.testing_utility import (testing_temp, get_baci_path,
    testing_input, compare_strings)


def create_cantilever(convergence_base_dir, subpath, n_el):
    """
    Create a cantilever beam for a convergence analysis.
    """

    input_file = InputFile()
    set_header_static(input_file, time_step=0.25, n_steps=4)
    set_runtime_output(input_file, output_energy=True)
    mat = MaterialReissner(radius=0.1, youngs_modulus=10000.0)
    beam_set = create_beam_mesh_line(input_file, Beam3rHerm2Line3, mat,
        [0, 0, 0], [1, 0, 0], n_el=n_el)
    input_file.add(BoundaryCondition(beam_set['start'],
        ('NUMDOF 9 ONOFF 1 1 1 1 1 1 0 0 0 VAL 0 0 0 0 0 0 0 0 0 ' +
        'FUNCT 0 0 0 0 0 0 0 0 0'),
        bc_type=mpy.bc.dirichlet
        ))
    fun = Function('COMPONENT 0 FUNCTION t')
    input_file.add(fun, BoundaryCondition(beam_set['end'],
        ('NUMDOF 9 ONOFF 0 0 1 0 0 0 0 0 0 VAL 0 0 -{} 0 0 0 0 0 0 ' +
        'FUNCT 0 0 {} 0 0 0 0 0 0'),
        format_replacement=[0.5, fun],
        bc_type=mpy.bc.dirichlet
        ))
    base_dir = os.path.join(convergence_base_dir, subpath)
    if os.path.isdir(base_dir):
        shutil.rmtree(base_dir)
    os.makedirs(base_dir)
    input_file_path = os.path.join(base_dir, 'cantilever.dat')
    input_file.write_input_file(input_file_path, add_script_to_header=False,
        header=False)
    return input_file_path


class TestSimulationManager(unittest.TestCase):
    """Test various stuff of the simulation manager."""

    def setUp(self):
        """
        This method is called before each test and sets the default meshpy
        values for each test. The values can be changed in the individual
        tests.
        """

        # Set default values for global parameters.
        mpy.set_default_values()

    def create_simulations(self, convergence_base_dir):
        """
        Create the convergence study and return the SimulationManager
        """

        os.makedirs(convergence_base_dir, exist_ok=True)
        manager = SimulationManager(convergence_base_dir)
        manager.add(Simulation(
            create_cantilever(convergence_base_dir, 'ref', 40), n_proc=4))
        simulations = [Simulation(
            create_cantilever(os.path.join(convergence_base_dir, 'sim'),
                str(n_el), n_el)) for n_el
            in range(1, 7, 2)]
        manager.add(simulations)
        return manager

    def check_convergence_results(self, convergence_base_dir):
        """
        Check the results from the convergence study.
        """

        sim_dir = os.path.join(convergence_base_dir, 'sim')
        sub_dirs = [os.path.join(sim_dir, path) for path in os.listdir(sim_dir)
            if os.path.isdir(os.path.join(sim_dir, path))]
        sub_dirs.append(os.path.join(convergence_base_dir, 'ref'))
        results = {}
        for sub_dir in sub_dirs:
            my_data = np.genfromtxt(sub_dir + '/xxx_energy.csv',
                delimiter=',')
            key = sub_dir.split('/')[-1]
            results[key] = my_data[-1, 2]

        results_ref = {'5': 0.335081498526998, '3': 0.335055487040675,
            '1': 0.33453718896204, 'ref': 0.335085590674607}
        for key in results_ref.keys():
            self.assertTrue(abs(results[key] - results_ref[key]) < 1e-12)

    def xtest_simulation_manager(self, command):
        """
        Create a convergence study and check the results.

        Args
        ----
        command: str
            Name of the command that is used to run the study.
        """

        # Create the simulations.
        convergence_base_dir = os.path.join(testing_temp, 'simulation_manager')
        manager = self.create_simulations(convergence_base_dir)

        if shutil.which(command) is not None:
            if command == 'mpirun':
                manager.run_simulations_and_wait_for_finish(
                    baci_build_dir=os.path.dirname(get_baci_path()),
                    status=False)
            else:
                manager.submit_batch_files_and_wait_for_finish(
                    baci_build_dir=os.path.dirname(get_baci_path()),
                    check_interval=1)
        else:
            self.skipTest('{} was not found'.format(command))

        # Check the results.
        self.check_convergence_results(convergence_base_dir)

    def test_simulation_manager_mpirun(self):
        """
        Create a convergence study and check the results. The simulations are
        run with mpirun.
        """
        self.xtest_simulation_manager('mpirun')

    def test_simulation_manager_cluster(self):
        """
        Create a convergence study on the cluster and check the results.
        """
        self.xtest_simulation_manager('sbatch')

    def test_batch_file(self):
        """
        Test that the created batch file is correct.
        """

        sim = Simulation(os.path.join(testing_temp, 'dummy.dat'),
            n_proc=37,
            n_nodes=69,
            exclusive=True,
            output_prefix='xxxXXXxxx',
            wall_time='66:66:66',
            restart_step=17,
            restart_dir='../old_sim',
            restart_from_prefix='xxx_old',
            job_name='awsome_job',
            feature='skylake')

        sim.create_batch_file(testing_temp, 'batch_name')

        compare_strings(self,
            'test_batch_file',
            os.path.join(testing_temp, 'batch_name'),
            os.path.join(testing_input, 'batch_name_reference'))


if __name__ == '__main__':
    # Execution part of script.
    unittest.main()
