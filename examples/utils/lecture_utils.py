# The MIT License (MIT)
#
# Copyright (c) 2018-2025 BeamMe Authors
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
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
"""This file contains functionality for the lecture lab."""

import subprocess
from pathlib import Path

from beamme.four_c.header_functions import set_header_static, set_runtime_output
from beamme.four_c.input_file import InputFile
from beamme.four_c.run_four_c import clean_simulation_directory


def run_four_c(
    *,
    mesh=None,
    simulation_name=None,
    total_time=1.0,
    n_steps=1,
    tol_residuum=1e-10,
):
    """Run a 4C simulation with given parameters."""

    # Setup the input file with mesh and parameters.
    input_file = InputFile()
    input_file.add(mesh)
    set_header_static(
        input_file,
        total_time=total_time,
        n_steps=n_steps,
        max_iter=20,
        tol_residuum=tol_residuum,
        tol_increment=1,
        create_nox_file=False,
        predictor="ConstDis",
    )
    set_runtime_output(
        input_file,
        output_solid=False,
        output_stress_strain=False,
        btsvmt_output=False,
        btss_output=False,
        output_triad=True,
        every_iteration=False,
        absolute_beam_positions=True,
        element_owner=True,
        element_gid=True,
        element_mat_id=True,
        output_energy=False,
        output_strains=True,
    )
    input_file["IO/RUNTIME VTK OUTPUT/BEAMS"]["MATERIAL_FORCES_GAUSSPOINT"] = True

    # Dump the file to disc.
    simulation_directory = Path.cwd() / simulation_name
    input_file_path = simulation_directory / f"{simulation_name}.4C.yaml"
    clean_simulation_directory(simulation_directory)
    input_file.dump(input_file_path)

    # Run the simulation and process the results line by line.
    with open(simulation_directory / f"{simulation_name}.log", "w") as logfile:
        # Command to run 4C
        four_c_exe = "/data/a11bivst/dev/4C/release/4C"
        command = [four_c_exe, input_file_path.absolute(), simulation_name]

        # Start simulation
        print("Start simulation")
        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,  # merge stderr into stdout (optional)
            text=True,  # get str instead of bytes
            bufsize=1,  # line-buffered
            cwd=simulation_directory,
        )

        # Process the output line by line
        nonlinear_solver_step_count = 0
        is_error = False
        finished = False
        for line in process.stdout:
            line = line.rstrip("\n")

            # Write line to logfile
            logfile.write(line + "\n")

            # Flush file so log is always up to date
            logfile.flush()

            # Process the line however you want
            if "Nonlinear Solver Step" in line:
                nonlinear_solver_step_count = int(line.split(" ")[4])
            elif "||F||" in line:
                if not nonlinear_solver_step_count == 0:
                    residuum = float(line.split(" ")[2])
                    print(
                        f"  Nonlinear Solver Step {nonlinear_solver_step_count}: Residuum = {residuum:.3e}"
                    )
            elif "Finalised step" in line:
                split = line.split(" ")
                step = int(split[2])
                time = float(split[7])
                print(f"Finished time step {step} for time {time:.3e}")
            elif "OK (0)" in line:
                finished = True
            elif (
                "========================================================================="
                in line
                and not finished
            ):
                if is_error:
                    print(line)
                is_error = not is_error

            if is_error:
                print(line)

        _return_code = process.wait()
