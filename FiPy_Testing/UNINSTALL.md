To completely delete a Conda environment along with all its installed dependencies, you can use the following steps:

---

### 1. **List All Conda Environments**
First, list all available Conda environments to confirm the name of the one you want to delete:

```bash
conda env list
```
This will display a list of all environments, e.g.:
```
# conda environments:
#
base                  /home/user/miniconda3
fipy-env              /home/user/miniconda3/envs/fipy-env
```

---

### 2. **Deactivate the Environment**
Make sure the environment you want to delete is not active. If it is active, deactivate it by running:

```bash
conda deactivate
```

---

### 3. **Delete the Environment**
To delete the environment (e.g., `fipy-env`), run:

```bash
conda env remove -n fipy-env
```

Here:
- `-n fipy-env` specifies the name of the environment you want to delete.

Alternatively, if the environment is activated, you can use:

```bash
conda env remove --all
```

---

### 4. **Confirm Deletion**
You can verify that the environment is deleted by listing environments again:

```bash
conda env list
```

---

### 5. **Remove Leftover Files (Optional)**
Sometimes, deleting an environment leaves behind cached files or configurations. To completely clean up:
- Navigate to the environments folder:
  ```bash
  cd ~/miniconda3/envs/
  ```
- Delete any folders related to the removed environment:
  ```bash
  rm -rf fipy-env
  ```

This step ensures that all residual files are gone.