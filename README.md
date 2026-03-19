# RDSCM

`RDSCM` implements a Region-Dependent Nonstationary Spatial Correlation Model for modeling and predicting spatial correlations of earthquake ground-motion residuals across California seismic regions.

The package supports a range of intensity measures, including:

- Peak Ground Acceleration (`PGA`)
- Peak Ground Velocity (`PGV`)
- Pseudo Spectral Accelerations (`PSAs`)

For additional methodological details, please refer to the [GIRS report](https://doi.org/10.34948/N38C7X).

---

## Installation

Before installing `RDSCM`, make sure that the dependency package `sf` is installed.

### 1. Install `sf`

```r
install.packages("sf")
```

### 2. Install `devtools`

```r
install.packages("devtools")
```

### 3. Install `RDSCM` from GitHub

```r
library(devtools)
install_github("wltcwpf/RDSCM")
```

---

## Usage

The main function in this package is `CA_spa_cor()`.

```r
CA_spa_cor(location_1, location_2, period, alpha = 0.05)
```

### Arguments

#### `location_1`
A numeric vector of length 2 specifying the first site location:

- the first element is longitude
- the second element is latitude

Coordinates must be given in decimal degrees using the WGS84 coordinate reference system (`EPSG:4326`).

> Note: This function currently supports only locations within California, USA.

#### `location_2`
A numeric vector of length 2 specifying the second site location, using the same format and coordinate reference system as `location_1`.

#### `period`
A numeric value indicating the oscillator period.

Accepted values are:

- `-1` for Peak Ground Velocity (`PGV`)
- `0` for Peak Ground Acceleration (`PGA`)
- any value in `[0, 10]` for Pseudo Spectral Acceleration (`PSA`) periods

#### `alpha`
A numeric value specifying the significance level for the confidence interval.

Default:

```r
alpha = 0.05
```

---

## Output

The function returns a data frame with two columns:

- `estimate`: the estimated spatial correlation between the two locations
- `interval`: a character string representing the confidence interval, formatted as `"[lower, upper]"`

---

## Example

```r
library(RDSCM)

location_1 <- c(-122.4194, 37.7749)  # San Francisco
location_2 <- c(-121.8863, 37.3382)  # San Jose

spa_cor <- CA_spa_cor(location_1, location_2, period = -1)

print(spa_cor)
```

---

## Reporting Issues

If you encounter any bugs, errors, or unexpected results, please report them on the [Issues](https://github.com/wltcwpf/RDSCM/issues) page.

Feedback and suggestions are also very welcome.

