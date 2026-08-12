import yaml
import sys

def simplify_dependencies(yml_file, output_file):
    with open(yml_file, 'r') as file:
        yml_dict = yaml.safe_load(file)

    dependencies = yml_dict['dependencies']
    simplified_dependencies = []

    for dep in dependencies:
        # remove version information
        simplified_dependency = dep.split('=')[0]
        simplified_dependencies.append(simplified_dependency)

    yml_dict['dependencies'] = simplified_dependencies

    with open(output_file, 'w') as file:
        yaml.safe_dump(yml_dict, file)

if __name__ == "__main__":
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    simplify_dependencies(input_file, output_file)
