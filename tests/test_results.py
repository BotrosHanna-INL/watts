# SPDX-FileCopyrightText: 2022 UChicago Argonne, LLC
# SPDX-License-Identifier: MIT

from datetime import datetime
from pathlib import Path

import watts
import numpy as np
from astropy import units as u


def test_results_openmc(run_in_tmpdir):
    params = watts.Parameters(city='Chicago', population=2.7e6)
    now = datetime.now()

    # Create some fake input files
    geom_xml = Path('geometry.xml')
    geom_xml.touch()
    mats_xml = Path('material.xml')
    mats_xml.touch()
    inputs = [geom_xml, mats_xml]

    # Create fake output files
    sp = Path('statepoint.50.h5')
    sp.touch()
    log_file = Path('OpenMC_log.txt')
    log_file.write_text("this is output\n")
    outputs = [sp, log_file]


    results = watts.ResultsOpenMC(params, now, inputs, outputs)

    # Sanity checks
    assert results.plugin == 'OpenMC'
    assert results.parameters == params
    assert results.time == now
    assert results.inputs == inputs
    assert results.outputs == outputs
    assert results.stdout == "this is output\n"

    # Other attributes
    assert len(results.statepoints) == 1

    # Saving
    p = Path('myresults.h5')
    results.save(p)
    assert p.is_file()

    # Ensure results read from file match
    new_results = watts.Results.from_hdf5(p)
    assert isinstance(new_results, watts.ResultsOpenMC)
    assert new_results.parameters == results.parameters
    assert new_results.time == results.time
    assert new_results.inputs == results.inputs
    assert new_results.outputs == results.outputs
    assert new_results.stdout == results.stdout


def test_results_moose(run_in_tmpdir):
    params = watts.Parameters(city='Chicago', population=2.7e6)
    now = datetime.now()

    # Create some fake input files
    moose_inp = Path('MOOSE.i')
    moose_inp.touch()
    inputs = [moose_inp]

    # Create fake output files
    csv = Path('MOOSE_csv.csv')
    csv.write_text("""\
prop1,prop2
3.5,1
4.0,2
5.0,3
    """)
    stdout = Path('MOOSE_log.txt')
    stdout.write_text('MOOSE standard out\n')
    outputs = [csv, stdout]

    results = watts.ResultsMOOSE(params, now, inputs, outputs)

    # Sanity checks
    assert results.plugin == 'MOOSE'
    assert results.parameters == params
    assert results.time == now
    assert results.inputs == inputs
    assert results.outputs == outputs

    # Other attributes
    assert results.stdout == 'MOOSE standard out\n'
    np.testing.assert_equal(results.csv_data['prop1'], [3.5, 4.0, 5.0])
    np.testing.assert_equal(results.csv_data['prop2'], [1, 2, 3])

    # Saving
    p = Path('myresults.h5')
    results.save(p)
    assert p.is_file()

    # Ensure results read from file match
    new_results = watts.Results.from_hdf5(p)
    assert isinstance(new_results, watts.ResultsMOOSE)
    assert new_results.parameters == results.parameters
    assert new_results.time == results.time
    assert new_results.inputs == results.inputs
    assert new_results.outputs == results.outputs
    assert new_results.stdout == results.stdout


def test_results_unit_conversion(run_in_tmpdir):
    # Uses Astropy for unit conversion
    u.imperial.enable()    # Enable imperial units
    Quantity = u.Quantity

    # Create a fake template file
    file = open("sam_template", "w")
    file.close()

    params = watts.Parameters()

    # Test with various unit conversion formats
    params['He_inlet_temp'] = Quantity(600, "Celsius")  # 873.15 K
    params['He_cp'] = Quantity(4.9184126, "BTU/(kg*K)") # 5189.2 J/kg-K
    params['He_Pressure'] = Quantity(7.0, "MPa") # 7e6 Pa
    params['Height_FC'] = Quantity(2000, "mm") # 2 m

    # Check that unit conversion in the MOOSE plugin is correct
    moose_plugin = watts.PluginMOOSE('sam_template', show_stderr=True)
    params_moose = moose_plugin.convert_unit(params, unit_system='si', unit_temperature='K')

    assert params_moose["He_inlet_temp"] == 873.15
    assert round(params_moose["He_cp"], 1) == 5189.2
    assert params_moose["He_Pressure"] == 7000000.0
    assert params_moose["Height_FC"] == 2.0    # Check that this parameter is in unit of 'm'

    # Check that unit conversion in the openmc plugin is correct
    openmc_plugin = watts.PluginOpenMC("foo")
    params_openmc = openmc_plugin.convert_unit(params, unit_system='cgs', unit_temperature='K')

    assert params_moose["He_inlet_temp"] == 873.15
    assert round(params_moose["He_cp"], 1) == 5189.2
    assert params_moose["He_Pressure"] == 7000000.0
    assert params_openmc["Height_FC"] == 200.0  # Check that this parameter is in unit of 'cm'