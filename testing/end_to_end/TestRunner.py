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


if __name__ == '__main__':
    main()
