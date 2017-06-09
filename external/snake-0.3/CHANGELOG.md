# Snake Change Log

---

## Current development

---

## 0.3

---

### Added
* Experimental and computational results for flapping wing at Re=75.
* Tests for class `CuIBMSimulation`.

### Changed
* Use python module `unittest` to rewrite tests.

### Fixed
* PetIBM simulation: write 2D grid with correct header.
* Convergence: return x-stations when getting the field values along an horizontal gridline.
* CuIBM: reshape correctly the pressure field when reading from file.

### Removed
* Pressure and vorticity images obtained by Anush with cuIBM.
* Snake forces obtained by Anush with cuIBM.

---

## 0.2

---

### Added

* Compatibility with Python-3.5. (Note that PETSc Python scripts are not compatible with Python-3.5.)
* A change log.
* Compatibility with VisIt-2.12.1 for plotting the fields from an IBAMR solutions.

### Changed

* Review examples.

### Fixed


### Removed

