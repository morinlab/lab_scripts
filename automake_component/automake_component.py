#!/usr/bin/env python

"""
automake_component.py
=====================
This script  generates a Pipeline Factory component
for a given Python script that uses argparse.

Inputs:
- Python script using argparse as interface

Outputs:
- A automatically generated Pipeline Factory component

Requirements
------------
- autopep8 for code formatting (optional)
- The Python script must use argparse as its interface.
- The argument parsing "code block" should be together.
  In other words, the argparse.ArgumentParser call and
  the argparse.parse_args call should delimit the argu-
  ment parsing code.

Known Issues
------------
- It's not trivial to determine if an argument is an input
  file/directory, an output file/directory or a parameter.
  Hence, the script prompts the user to classify each
  command-line argument (as input, output or param).
  Optionally, there's an automatic mode (--automatic) that
  uses simple heuristics to classify the arguments. It is
  recommended that the user verifies the automatic
  classification.
"""

import argparse
import os.path
import time
import shutil


def main():

    # Argument parsing
    parser = argparse.ArgumentParser(description="Generate a Pipeline Factory component "
                                     "from a Python script that uses argparse.")
    parser.add_argument("input_script", help="Input Python script")
    parser.add_argument("--output_dir", default=".", help="Output directory")
    parser.add_argument("--automatic", action="store_true",  help="Automatically determine "
                        "classification of command-line arguments (input, output or param)")
    args = parser.parse_args()

    # Initialize variables
    script_path = args.input_script
    script_filename = os.path.basename(script_path)
    if script_filename.endswith(".py"):
        script_name = script_filename[:-3]
    else:
        script_name = script_filename
    properties = {"description": "",
                  "input_files_dict": {},
                  "input_files_list": "",
                  "opt_args": {},
                  "output_files_dict": {},
                  "output_files_list": "",
                  "params_dict": {},
                  "params_list": "",
                  "parser_block": "",
                  "pos_args": [],
                  "script_filename": script_filename,
                  "script_name": script_name,
                  "script_version": "v1.0.0",
                  "today": time.strftime("%b %d %Y")}

    # Script parsing
    in_parser = False
    parser_block = ""
    parser_var = ""
    with open(script_path) as infile:
        for line in infile:
            # Look for argparse import
            if "import argparse" in line or "from argparse" in line:
                parser_block += line
            # Look for version definition in script
            elif "__version__" in line and "=" in line:
                var, value = line.split("=")
                if type(eval(value)) is str:
                    properties["script_version"] = eval(value)
            # Look for beginning of argument parsing code block
            elif "ArgumentParser" in line:
                in_parser = True
                parser_block += line.lstrip()
                parser_var = line.strip().split("=")[0].strip()
            # Stop once the parse_args call is reached
            elif parser_var + ".parse_args(" in line:
                break
            # Append every line of argument parsing code block
            elif in_parser:
                parser_block += line.lstrip()
    # Add parser block for component_ui
    properties["parser_block"] = parser_block
    properties["parser_block"] += "\nargs, unknown = {}.parse_known_args()\n".format(parser_var)

    # Analyze parser
    exec(parser_block)
    parser = eval(parser_var)
    # Set description
    properties["description"] = parser.description
    # Iterate over parameters
    for action in parser._actions:
        # Skip help option
        if action.dest == "help":
            continue
        is_optional = False
        if len(action.option_strings) >= 1:
            is_optional = True
        if is_optional:
            properties["opt_args"][action.dest] = action.option_strings[0]
        else:
            properties["pos_args"].append(action.dest)
        if args.automatic:
            # Determine if input file, output file or param using heuristics
            is_file = False
            is_output = False
            # If file or dir is in the dest or if the type is argparse.FileType
            if ("file" in action.dest.lower() or "dir" in action.dest.lower() or
                    "input" in action.dest.lower() or "output" in action.dest.lower() or
                    isinstance(action.type, argparse.FileType)):
                is_file = True
            # If out is in dest or help, then mark as output
            if is_file and ("out" in action.dest.lower() or "out" in action.help.lower()):
                is_output = True
            # If nargs is 0, then the argument isn't a file
            if action.nargs == 0:
                is_file = False
            # If there is only one positional argument, then assume it's the input
            if (sum([int(len(a.option_strings) == 0) for a in parser._actions]) == 1 and
                    not is_optional):
                is_file = True
                is_output = False
        else:
            is_answered = False
            while not is_answered:
                question_template = ("\nIs the following command-line argument an input "
                                     "file/directory, an output file/directory or a parameter?\n"
                                     "    Name:  {}\n"
                                     "    Desc:  {}\n"
                                     "Select an option ([i]nput, [o]utput or [p]arameter): ")
                answer = raw_input(question_template.format(action.dest, action.help))
                if answer.lower().startswith("i"):
                    is_file = True
                    is_output = False
                    is_answered = True
                elif answer.lower().startswith("o"):
                    is_file = True
                    is_output = True
                    is_answered = True
                elif answer.lower().startswith("p"):
                    is_file = False
                    is_output = False
                    is_answered = True
                else:
                    print
                    print "*** Error: Invalid input. Try again.***"
        if is_file and not is_output:  # Input files
            if is_optional:
                properties["input_files_dict"][action.dest] = "__OPTIONAL__"
            else:
                properties["input_files_dict"][action.dest] = "__REQUIRED__"
        if is_file and is_output:  # Output files
            if is_optional:
                properties["output_files_dict"][action.dest] = "__OPTIONAL__"
            else:
                properties["output_files_dict"][action.dest] = "__REQUIRED__"
        if not is_file:  # Params
            if action.default is None and is_optional:
                properties["params_dict"][action.dest] = "__OPTIONAL__"
            elif action.default is None and not is_optional:
                properties["params_dict"][action.dest] = "__REQUIRED__"
            else:
                properties["params_dict"][action.dest] = action.default

    # Post processing
    properties["input_files_list"] = ", ".join(properties["input_files_dict"].keys())
    properties["output_files_list"] = ", ".join(properties["output_files_dict"].keys())
    properties["params_list"] = ", ".join(properties["params_dict"].keys())

    # Generating component
    template_dir = os.path.join(os.path.dirname(__file__), "component_template")
    component_dir = os.path.join(args.output_dir, script_name)
    shutil.copytree(template_dir, component_dir)
    for component_file in os.listdir(component_dir):
        if component_file.startswith("."):
            continue
        component_path = os.path.join(component_dir, component_file)
        incomp = open(component_path)
        incomp_contents = incomp.read()
        outcomp = open(component_path, "w")
        outcomp_contents = incomp_contents.format(**properties)
        try:
            import autopep8
            if component_file.endswith(".py"):
                outcomp_contents = autopep8.fix_code(outcomp_contents)
        except ImportError:
            pass
        outcomp.write(outcomp_contents)


if __name__ == '__main__':
    main()
