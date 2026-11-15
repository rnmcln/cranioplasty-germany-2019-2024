"""Run the full pipeline from the repository root: analysis, validation, figures."""
import subprocess
import sys

for step in ('code/01_analysis.py', 'code/03_validate.py', 'code/02_figures.py'):
    print(f'\n=== {step}')
    subprocess.run([sys.executable, step], check=True)
