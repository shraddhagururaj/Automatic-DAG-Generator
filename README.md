# Automatic-DAG-Generator - DIRECTED RESEARCH

Automatic DAG Generator is customized for Jupiter Orchestrator (available here: [https://github.com/ANRGUSC/Jupiter]).


`
### dummy_app_gen.py

This code uses the output DAG of ``rand_task_gen.py``(by Diyi), which is used to generate a random DAG, to generate the corresponding dummy application for Jupiter.This code has been updated to generate the tasks in C script. To run this script, you should set `--conf task_config.yml` as the parameter. It will generate the ``dummy_app`` folder and all required content which can be used as a sample application for [Jupiter](https://github.com/ANRGUSC/Jupiter). The ``dummy_app`` has been tested and work with **Jupiter Version 3 & 4**.

## User Instructions
Run the dummy_app_gen code to generate the C scripts:<br />
**python3 dummy_app_gen.py --task_config.yml**<br />
To execute the C scripts: <br />
**gcc -o taskx.c** <br />
**./taskx** <br />
 To run the python wrapper:<br />
**python3 taskx_wrapper.py**
 

## Acknowledgment
This material is based upon work supported by Defense Advanced Research Projects Agency (DARPA) under Contract No. HR001117C0053. Any views, opinions, and/or findings expressed are those of the author(s) and should not be interpreted as representing the official views or policies of the Department of Defense or the U.S. Government.

