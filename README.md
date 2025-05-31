# bands_NMREGIS

**bands_NMREGIS** is a Python package for analyzing and visualizing electronic band structures, designed to work with data from DFT calculations using VASP.

---

## 🚀 Features

- Load and process electronic band structure data
- Generate plots with Matplotlib
- Command-line interface support

---

## 📦 Installation

You can install the package using **Git**:

### Clone and install locally:

```bash
git clone https://github.com/yourusername/bands_NMREGIS.git
cd bands_NMREGIS
pip install .
```

### Or install directly with pip from GitHub:

```bash
pip install git+https://github.com/yourusername/bands_NMREGIS.git
```

For development mode (editable install):

```bash
git clone https://github.com/yourusername/bands_NMREGIS.git
cd bands_NMREGIS
pip install -e .
```

---

## 🛠️ Usage

### As a script:

```bash
python bands_NMREGIS.py
```

### As a command-line tool (after installation):

```bash
bands_NMREGIS
```

---

## 📄 Example

```python
import bands_NMREGIS

# Example function call if exposed via API
bands_NMREGIS.run_analysis("path/to/data/file.dat")
```

---

## 🧩 Requirements

- Python 3.8+
- numpy
- matplotlib

(Dependencies will be installed automatically via `pip`.)

---

## 📁 Project Structure

```
bands_NMREGIS/
├── bands_NMREGIS.py
├── setup.py
├── README.md
├── LICENCE
```

---

## 📃 License

GPL-v3 License. See `LICENSE` file for details.

---

## 👤 Author

- Natan M. Regis (natan.moreira.regis12@gmail.com)

Feel free to contribute or open issues!