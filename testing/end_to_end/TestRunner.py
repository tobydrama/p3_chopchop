import argparse
import filecmp
import json
import os
import subprocess


def main():
    with open('test_data/index.json') as config:
        data = json.load(config)

    for name, info in data.items():
        print(f"Testing '{name}' (chopchop.py {info['args']})... ", end='')
        try:
            subprocess.run(
                ['python', './chopchop.py', '-o', 'testing/end_to_end/temp', '-BED', '-J', '--logLevel', 'DEBUG']
                + info['args'].split(' '), capture_output=True, cwd='../../', check=True)
            # print(result.stdout.decode())
            if filecmp.cmp('temp/results.bed', "test_data/" + info['results']):
                print(f"\r'{name}'\tSUCCESS: Output matches.")
            else:
                print(f"\r'{name}'\tERROR: Output does not match.")
        except subprocess.CalledProcessError as e:
            print(f"\r'{name}'\tERROR:\t{e.stderr.decode()}")
            pass

    # Empty temp folder run.

    for filename in os.listdir('temp/'):
        path = os.path.join('temp/', filename)
        if os.path.isfile(path):
            os.unlink(path)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--showDiff', action='store_false', dest='show_diffs',
                        help="Show which files did not match on test failure.")
    return parser.parse_args()


def compare_results(loc: str, verbose: bool = False) -> any:
    diff = []
    for filename in os.listdir('temp'):
        if not (os.path.isfile(loc + filename) and filecmp.cmp('temp/' + filename, loc + filename)):
            diff.append(filename)




def run_test(name: str, info: dict, show_diffs: bool = False) -> bool:
    print(f"Testing '{name}' (chopchop.py {info['args']})... ", end='')

    command = ['python', './chopchop.py', '-o', 'testing/end_to_end/temp', '-BED', '-J', '--logLevel', 'DEBUG'] \
              + info['args'].split(' ')

    try:
        subprocess.run(command, capture_output=True, cwd='../../', check=True)

        if compare_results(info['results']):
            print(f"\r'{name}' SUCCESS: Output matches")
            return True
        else:
            if show_diffs:
                print(f"\r'{name}' ERROR: Files do not match: {get_diffs(info['results'])}.")
            else:
                print(f"\r'{name}' ERROR: Output does not match.")

    except subprocess.CalledProcessError as e:
        print(f"\r'{name}' ERROR:\nSTDOUT:\n{e.stdout.decode()}\n\nSTDERR:\n{e.stderr.decode()}\n\n")

    return False


def new_testrunner(args):
    with open('test_data/index.json') as config:
        data = json.load(config)

    successes = 0
    for name, info in data.items():
        if run_test(name, info, args.show_diffs):
            successes += 1

    print(f"{successes} / {len(data)} tests passed.")


if __name__ == '__main__':
    main()
