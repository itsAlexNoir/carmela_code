# 🌟 Carmela Code

> ✨ A Python library for solving the one-electron one-dimensional Time-Dependent Schrödinger Equation (TDSE)

Welcome to **Carmela Code**! 👋 This friendly library lets you explore quantum dynamics in a fun, interactive way. Whether you're simulating laser-matter interactions or just curious about how electrons behave in time-varying fields, Carmela has got you covered! 🚀

## 🎯 What is this?

Carmela Code solves the **time-dependent Schrödinger equation** for a one-electron system in one dimension. It's a toy model designed for rapid prototyping and experimentation—perfect for researchers, students, or anyone who wants to play with quantum mechanics without getting bogged down in complexity! 🎓⚛️

## ✨ Features

- 🔬 **Full TDSE solver** with customizable Hamiltonians
- 💡 **Laser pulse generation** with flexible envelope shapes
- 📊 **Built-in visualization** for wavefunctions and observables
- 🎨 **Beautiful plots** powered by matplotlib and palettable
- 📐 **Momentum analysis** via Fourier transforms and t-SURFF
- 🛡️ **Absorbing boundaries** to prevent reflection artifacts
- ⚡ **Fast and flexible** numerical propagation

## 🚀 Quick Start

### Installation

Clone this repo and install dependencies:

```bash
git clone <your-repo-url>
cd carmela_code
pip install -e .
```

### Usage

Import the library in your Python scripts:

```python
from carmela import *
```

Check out the example in `calculations/example_XUV/` to see Carmela in action! 🎬

## 📦 Modules

Carmela is organized into clean, focused modules:

- 🧮 **`absorber`** – Absorbing boundary conditions
- 📏 **`axes`** – Spatial and momentum grids
- 🔢 **`constants`** – Physical constants in atomic units
- 📐 **`fdrule`** – Finite-difference rules
- 🌊 **`flux`** – Flux calculations for momentum distributions
- 📊 **`graph`** – Visualization utilities
- ⚛️ **`hamiltonian`** – Time-dependent Hamiltonians
- 📄 **`input_reader`** – Parse YAML/INP configuration files
- 💡 **`laser`** – Laser field generation
- 🚀 **`momentum`** – Momentum and kinetic energy analysis
- 🌀 **`wavefunction`** – Wavefunction manipulation and observables

## 🛠️ Requirements

- Python 3.13+
- numpy
- matplotlib
- palettable

All dependencies are managed automatically when you install the package! 📦

## 🎓 Example

Run the XUV ionization example:

```bash
cd calculations/example_XUV
python example_XUV.py
```

This will propagate a hydrogen ground state under an XUV pulse and calculate ionization yields. The configuration is in `input.yaml` for easy tweaking! ⚙️

## 🤝 Contributing

Got ideas? Found a bug? Want to add features? Contributions are welcome! Feel free to open issues or submit pull requests. Let's make quantum simulations more accessible together! 💪

## 📝 License

This project is open source. Check the repository for license details.

## 💬 Questions?

If you have questions or need help getting started, feel free to reach out or open an issue. Happy simulating! 🎉

---

Made with ❤️ and ☕ for the quantum physics community
