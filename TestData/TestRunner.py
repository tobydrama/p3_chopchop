import json
import filecmp
import os
import subprocess


def main():
    with open('config.json') as config:
        data = json.load(config)

    for name, info in data.items():
        print(f"Testing '{name}' (chopchop.py {info['args']})... ", end='')
        try:
            subprocess.run(['python', './chopchop.py', '-o', 'TestData/temp', '-BED', '-J'] + info['args'].split(' '),
                           capture_output=True, cwd='../', check=True)
            if filecmp.cmp('temp/results.bed', info['results']):
                print(f"\r{name}\tSUCCESS: Output matches.")
            else:
                print(f"\r{name}\tERROR: Output does not match.")
        except subprocess.CalledProcessError as e:
            print(f"{name} failed:")
            print(e.stderr.decode())
            pass
    
    # Empty temp folder run.
    for filename in os.listdir('temp/'):
        path = os.path.join('temp/', filename)
        if os.path.isfile(path):
            os.unlink(path)


if __name__ == '__main__':
    main()
