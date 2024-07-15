#! /usr/bin/env python
# To remove relion jobs' intermediate files, only map end with .mrc(s) will be removed.
import os
import sys

# Use argparse to handle command line arguments
import argparse

parser = argparse.ArgumentParser(description="Clean up intermediate files from relion jobs")
parser.add_argument('--dry', action='store_true', help="Enable dry run (won't actually delete files)")
args = parser.parse_args()

def get_all_jobfolders(path):
    # Recursively find all job folders under the given path
    for root, dirs, files in os.walk(path):
        if 'job' in root:
            yield root
        for dir in dirs:
            yield from get_all_jobfolders(os.path.join(root, dir))

def remove_intermediates(path, dry_run=False):
    # Remove intermediate files from a given directory
    print(f'  in {path}')
    files = os.listdir(path)
    iterations = [int(f.split('_it')[1].split('_')[0]) for f in files if '_it' in f]
    if not iterations:
        return
    max_iteration = max(iterations)
    print(f'    detected {max_iteration} iterations')
    for filename in files:
        if filename.endswith(('.mrc', '.mrcs')) and '_it' in filename and f'it{max_iteration:03}' not in filename:
            print(f'      removing {os.path.join(path, filename)}')
            if not dry_run:
                os.remove(os.path.join(path, filename))
    print('    remove completed')

def main():
    if not args.dry:
        command = input('Not in dry mode, print y to continue: ')
        if command.lower() != 'y':
            sys.exit()
    if os.getcwd().split('/')[-1].startswith('job'):
        remove_intermediates('.', args.dry)
    else:
        job_folders = list(get_all_jobfolders('.'))
        print(f'  detected {len(job_folders)} job folders')
        for job_folder in job_folders:
            remove_intermediates(job_folder, args.dry)

if __name__ == '__main__':
    main()

