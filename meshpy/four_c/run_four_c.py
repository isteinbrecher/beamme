# -*- coding: utf-8 -*-
# -----------------------------------------------------------------------------
# MeshPy: A beam finite element input generator
#
# MIT License
#
# Copyright (c) 2018-2024
#     Ivo Steinbrecher
#     Institute for Mathematics and Computer-Based Simulation
#     Universitaet der Bundeswehr Muenchen
#     https://www.unibw.de/imcs-en
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
"""Provide a function that allows to run 4C"""


import os
import sys
import subprocess


def run_four_c(
    input_file,
    output_dir,
    *,
    four_c_exe=None,
    mpi_command=None,
    n_proc=None,
    output_name="xxx",
    restart_step=None,
    restart_from=None,
    log_to_console=False,
):
    """Run a 4C simulation and return the exit code of the run

    This function looks into the environment variables for some parameters:
        "MESHPY_FOUR_C_EXE"
        "MESHPY_MPI_COMMAND"
        "MESHPY_MPI_NUM_PROC"
    If the corresponding keyword arguments are set, they overwrite the environment
    variable.

    If an environment variable "MESHPY_RUN_4C_WITH_QUEENS" exists, then the QUEENS
    "run_subprocess_with_logging" functionality is used to run 4C. This additionally
    requires the environment variable "QUEENS_JOB_ID".

    Args
    ----
    input_file: str
        Path to the dat file on the filesystem
    output_dir: str
        Directory where the simulation should be performed (will be created if
        it does not exist)
    four_c_exe: str
        Optionally explicitly specify path to the 4C executable
    mpi_command: str
        Command to launch MPI, defaults to "mpirun"
    n_proc: int
        Number of process used with MPI, defaults to 1
    output_name: str
        Base name of the output files
    restart_step: int
        Time step to restart from
    restart_from: str
        Path to initial simulation (relative to output_dir)
    log_to_console: bool
        If the 4C simulation output should be shown in the console. This has no
        effect if 4C is run with QUEENS.

    Return
    ----
    return_code: int
        Return code of 4C run
    """

    def get_env_variable(name, *, fail_if_not_set=False, default=None):
        """Return the value of an environment variable"""
        if name in os.environ.keys():
            return os.environ[name]
        elif fail_if_not_set:
            raise ValueError(f"Environment variable {name} is not set")
        return default

    # Fist get all needed parameters
    if four_c_exe is None:
        four_c_exe = get_env_variable("MESHPY_FOUR_C_EXE")
    if mpi_command is None:
        mpi_command = get_env_variable("MESHPY_MPI_COMMAND", default="mpirun")
    if n_proc is None:
        n_proc = get_env_variable("MESHPY_MPI_NUM_PROC", default="1")

    # Setup paths and actual command to run
    os.makedirs(output_dir, exist_ok=True)
    log_file = os.path.join(output_dir, output_name + ".log")
    error_file = os.path.join(output_dir, output_name + ".err")
    command = mpi_command.split(" ") + [
        "-np",
        str(n_proc),
        four_c_exe,
        input_file,
        output_name,
    ]
    if restart_step is None and restart_from is None:
        pass
    elif restart_step is not None and restart_from is not None:
        command.extend([f"restart={restart_step}", f"restartfrom={restart_from}"])
    else:
        raise ValueError(
            "Provide either both or no argument of [restart_step, restart_from]"
        )

    # Actually run the command, depending on the method provided
    run_with_queens = (
        get_env_variable("MESHPY_RUN_4C_WITH_QUEENS", default=None) is not None
    )
    if not run_with_queens:
        with open(log_file, "w") as stdout_file, open(error_file, "w") as stderr_file:
            process = subprocess.Popen(
                command,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                cwd=output_dir,
                text=True,
            )

            for stdout_line in process.stdout:
                if log_to_console:
                    sys.stdout.write(stdout_line)
                stdout_file.write(stdout_line)

            for stderr_line in process.stderr:
                if log_to_console:
                    sys.stderr.write(stderr_line)
                stderr_file.write(stderr_line)

            process.stdout.close()
            process.stderr.close()
            return_code = process.wait()
        return return_code

    else:
        from queens.utils.run_subprocess import run_subprocess_with_logging

        job_id = get_env_variable("QUEENS_JOB_ID")
        command_string = " ".join(command)
        command_string = f"cd {output_dir} && " + command_string
        return_code, process_id, stdout, stderr = run_subprocess_with_logging(
            command_string,
            terminate_expression="PROC.*ERROR",
            logger_name=__name__ + f"_{job_id}",
            log_file=log_file,
            error_file=error_file,
            streaming=False,
            raise_error_on_subprocess_failure=False,
        )
        return return_code
