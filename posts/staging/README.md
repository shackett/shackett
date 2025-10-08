## Preparing a post

### Setting up a Python environment

```python
uv venv --python 3.11
source .venv/bin/activate

# Core dependencies
# uv pip install "napistu[scverse]==0.5.4"
# Personal utilities package with genomics analysis functions
# uv pip install "git+https://github.com/shackett/shackett-utils.git@v0.1.2[all]" 
# Additional dependencies
uv pip install openpyxl ipykernel nbformat nbclient
python -m ipykernel install --name=blog-staging
```