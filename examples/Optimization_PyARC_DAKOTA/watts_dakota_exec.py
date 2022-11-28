# SPDX-FileCopyrightText: 2022 UChicago Argonne, LLC
# SPDX-License-Identifier: MIT

"""
This example demonstrates how to use DAKOTA to perform 
optimization with WATTS. This example uses PyARC as the
coupled code for Dakota. There are a number of files in 
this example. The primary files are 'watts_dakota_exec.py',
'watts_pyarc_exec.py', 'dakota_watts_opt.in', 'pyarc_template',
and 'dakota_driver.py', with several other secondary files
needed to run PyARC. Note that all of the files mentioned above 
can be templated with the 'extra_template_inputs' option in 
the Dakota plugin (refer below for details).

The execution of Dakota with WATTS is essentially a two-step process. 
When the user runs the 'watts_dakota_exec.py' script, it generates 
an input file for Dakota based on the 'dakota_watts_opt.in' template. 
It then runs Dakota which in turn drives the execution of the coupled 
code (PyARC) through 'dakota_driver.py' by running 'watts_pyarc_exec.py'.
The Dakota driver is responsible for facilitating the communication
between Dakota and the coupled code. Note that this is done through 
Dakota's `interfacing` library. The user needs to ensure that this 
library is available prior to running Dakota with WATTS.

The user first needs to identify the parameters that will be optimized 
by Dakota (These are the parameters that will be varied in PyARC).
In this example, the parameters are the 'assembly_pitch' and 
'assembly_length' and are represented as 'AP' and 'AL' in the 
'dakota_watts_opt.in' file. The input values of 'assembly_pitch' 
and 'assembly_length' are set to the values generated by Dakota in
'watts_pyarc_exec.py'.

The user then needs to identify the response parameters that Dakota
will use to judge the optimization. In this example, 
'KeffOpt', 'CoreWeight', and 'KeffCrit' are chosen as the response 
parameters, where they are represented as 'KO', 'CW', and 'KC', respectively.
Then, in the 'watts_pyarc_exec.py', the user has to set an object known as 
'dakota_descriptors' to 'params'. The 'dakota_descriptors' is a dictionary
that links 'KeffOpt' to 'KO', 'CoreWeight' to 'CW', and so on. The order of 
the descriptors MUST match the order of the response descriptors in 
'dakota_watts_opt.in'. 

After that, the user needs to save the results from PyARC as objects 
to 'params'. Note that the names of these objects must match the names 
of the response descriptors described above ('KeffOpt', 'CoreWeight', etc).
Lastly, 'params', MUST be saved as a pickle file with the name 'opt_res.out'
which will be loaded and returned to Dakota.
"""

import watts

# watts.Database.set_default_path('/default/directory') # Set default save directory if necessary

params = watts.Parameters()

# Provide optional parameters to generate the dakota_driver.py file
params['dakota_path'] = '/software/Workbench/dakota/share/dakota/Python/dakota' # Specify path to Dakota's 'interface' library
params['coupled_code_exec'] = 'watts_pyarc_exec.py' # Specify the script of the coupled code
params['dakota_driver_name'] = 'dakota_driver.py'   # Specify the file name of Dakota driver.
                                                    # If user wishes to template the Dakota driver file, 
                                                    # make sure to set the key as 'dakota_driver_name'
                                                    # because WATTS relies on this key name to do 
                                                    # additional internal checks to ensure Dakota
                                                    # has the right permission to run this file.
                                                    # An error will be raised if a different key is used here.
                                                    # Note: If users do not wish to template the Dakota driver file,
                                                    # make sure that the mode of path of the file is correct.
                                                    # If not, users can use os.chmod(<dakota_driver_file_name>, 0o755)
                                                    # to change it to the right mode.

# Dakota parameters
params['real'] = 2
params['temp'] = 26.85

# Response descriptors for Dakota (Make sure the descriptors match those in the watts_pyarc_exec.py)
params['dakota_descriptor_1'] = 'KO'
params['dakota_descriptor_2'] = 'CW'
params['dakota_descriptor_3'] = 'KC'

# File name of the final output from Dakota.
# The default name of the output file is 'dakota_opt.dat'
# If the user wishes to change the name of the output file, make sure to use the 'dakota_out_file' key.
params['dakota_out_file'] = 'dakota_opt.dat'

params.show_summary(show_metadata=False, sort_by='key')

# Dakota Workflow
# Use 'extra_template_inputs' to template optional extra files.
# Use 'extra_inputs' to copy files to the temporary directory for WATTS operation.
# Set 'auto_link_files = <string_name_for_files>' if users wish to automatically add all files 
# in 'extra_template_inputs' and 'extra_inputs' to the 'link_files' option in the Dakota 
# input file. In the Dakota input file, users must set "link_files = {{ <string_name_for_files> }}".
dakota_plugin = watts.PluginDakota(
    template_file='dakota_watts_opt.in',
    extra_template_inputs=['watts_pyarc_exec.py', 'dakota_driver.py'],
    extra_inputs=['pyarc_input.isotxs', 'pyarc_template', 'lumped.son', 'watts_dakota_exec.py'],
    auto_link_files = 'auto_link_file_string_name',
    show_stdout=True) # show all the output

dakota_result = dakota_plugin(params)

for key in dakota_result.output_data:
    print(key, dakota_result.output_data[key])