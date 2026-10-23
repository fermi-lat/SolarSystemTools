---
package: SolarSystemTools
summary: Livetime and exposure cubes binned in angle from a moving solar-system body (Sun or Moon), and the Sun/Moon intensity templates built from them for likelihood analysis (gtltcubesun, gtltsumsun, gtexpcubesun, gtexphpsun, gtsuntemp).
layer: analysis
depends_on: [Likelihood, map_tools, healpix, irfs, dataSubselector, astro, tip, st_app, st_stream, st_facilities, facilities, CLHEP, cfitsio]
used_by: []
last_reviewed_commit: baaf257
---

# SolarSystemTools Architecture

**SolarSystemTools** models emission from a source that moves across the
sky: the Sun (disk plus inverse-Compton halo) or the Moon. A normal
livetime cube (`gtltcube`) bins livetime by sky pixel and instrument
inclination θ. This package adds a second angle, **θ_sun**, the angle
between the sky pixel and the moving body at the time of each FT2
interval. From that cube it computes exposure binned in θ_sun. It then
folds the exposure with a radial intensity profile of the source into a
time-averaged intensity map, which `Likelihood` can use like a diffuse
map cube.

The package has a STATIC library (about 3,000 lines, of which roughly 400
are dead; see *Invariants and gotchas*) and five `st_app` tools (about
1,100 lines). It was written by G. Jóhannesson in 2012 and has had only
build-system changes since 2015. Recent history is the `auto_ptr`
removal, CMake and CI.

There is no Doxygen mainpage. `doc/release.notes` ends in 2015 and its
header says "package dataSubselector" (a copy-paste error). User help
for four of the tools is in `fermitools-fhelp` (`gtltcubesun.txt`,
`gtltsumsun.txt`, `gtexphpsun.txt`, `gtsuntemp.txt`). There is none for
`gtexpcubesun`.

## Purpose and scope

**Does:**
- Accumulate FT2 livetime into a HEALPix (equatorial, NESTED) array of
  `CosineBinner2D` cells indexed by (cos θ, θ_sun[, φ]) (`gtltcubesun`).
  For each FT2 interval, the body's direction is computed with
  `astro::SolarSystem` (JPL DE405) from the spacecraft position, so
  parallax is included (it matters for the Moon). Livetime is weighted by
  1/distance² for the solar disk, or for the whole sky for the Moon, and
  the time-weighted mean distance is recorded.
- Add such cubes together (`gtltsumsun`).
- Integrate the effective area over the cube to get exposure per sky
  pixel, energy and θ_sun bin. `gtexphpsun` does this on the cube's
  HEALPix grid, and `gtexpcubesun` on a WCS grid.
- Combine the θ_sun exposure, an ordinary `gtexpcube2` exposure and a
  source profile I(θ, E) into an intensity template (`gtsuntemp`).

**Does not:**
- Provide the source profile. `sunprof` is an external FITS file
  (`ANGLES`, `ENERGIES`, `SST_PROFILE`), for example the
  `solar_profile_v2r0.fits` and `lunar_profile_v2r0.fits` files named in
  the help.
- Fit anything or contain Python bindings.
- Track the Sun's own variability. The profile is fixed.

## Position in the system

| Direction | Package | Why |
|-----------|---------|-----|
| Depends on | `Likelihood` | `AppHelpers` (IRFs from `irfs`/`evtype`/CALDB, output-file checks, time-cut checks), `Observation` (wrapped by `SolarSystemTools::Observation`), `RoiCuts` (reads FT1 GTIs and DSS cuts, writes the GTI extension and DSS time keywords), `CountsMap` (`cmap` geometry), `BinnedExposure` (`avgexp` in `gtsuntemp`). For Likelihood internals see `src/Likelihood/doc/architecture/`. `gtltcubesun` is a modified copy of Likelihood's `gtltcube` (`makeExposureCube`). |
| Depends on | `map_tools` | Header-only use of the `BasicExposure<Sky, Cos>` template in `map_tools/Exposure.h`: `ExposureSun` derives from `BasicExposure<HealpixArray<CosineBinner2D>, CosineBinner2D>`. No compiled `map_tools` symbol is used. |
| Depends on | `healpix` | `HealpixArray`, `Healpix`, `HealpixArrayIO` for the cube's sky grid. `CosineBinner2D` replaces `healpix::CosineBinner`. |
| Depends on | `irfs` | `irfInterface::IEfficiencyFactor` (livetime-efficiency correction) and `Irfs`/`aeff()` through the Likelihood `Observation`. CMake links `irfLoader`. |
| Depends on | `dataSubselector` | `Gti` (copy GTIs to outputs) and `Cuts` (IRF version DSS keyword). It is **not** in this package's CMake; the headers and library arrive transitively through `Likelihood`'s `PUBLIC` links. |
| Depends on | `astro` | `SolarSystem` (body direction and distance), `JulianDate::missionStart`, `SkyDir`, `SkyProj`, `PointingTransform`. `ExposureSun.cxx` includes `EarthCoordinate.h` but does not use it. |
| Depends on | `tip`, `st_app`, `st_stream`, `st_facilities`, `facilities` | FITS I/O, tool framework, output, `Util::resolve_fits_files`/`file_ok`/`pixel2SkyDir`/`writeDateKeywords`, `Util::expandEnvVar`/`getDataPath`. |
| Depends on (external) | CLHEP, cfitsio | `Hep3Vector`/`HepRotation`. Direct `fits_*` calls in `ExposureCubeSun.cxx` (variable-length column format) and in the dead `LikeExposureSun.cxx`. Neither is declared in `CMakeLists.txt`; `gtexpcubesun` alone links `CLHEP::VectorS`. |
| Link only | `evtbin` | Linked `PUBLIC`, but nothing includes `evtbin/`. |
| Used by | — | No package includes `SolarSystemTools/` headers. `sane`'s `gt_apps` does not wrap these tools. |

## Key abstractions

```
CosineBinner2D : vector<double>          one sky pixel: bins in cosθ × θ_sun (× φ), sparse in θ_sun
        ^ cell of
healpix::HealpixArray<CosineBinner2D>    = SkyBinner2D
        ^ data of
map_tools BasicExposure<SkyBinner2D, CosineBinner2D>
        ^
ExposureSun        fill(dirz, dirsun, dirzenith, Δt, 1/d², distCosCut); FITS read/write
        ^ owns two (unweighted, livetime-fraction-weighted)
ExposureCubeSun    FT2 loop + astro::SolarSystem; value(dir, cosθsun, aeff[, E]); +=; writeFile
        ^ referenced by
Observation (SolarSystemTools)  = Likelihood::Observation + ExposureCubeSun&
        ^ used by
BinnedExposureSun (WCS, gtexpcubesun)     HealpixExposureSun (HEALPix, gtexphpsun; read by gtsuntemp)
                                                  ^
SolarProfile <- FitsSolarProfile          SolarTemplate (gtsuntemp)
```

| Class | Header | Role |
|-------|--------|------|
| `CosineBinner2D` | `CosineBinner2D.h` | Livetime histogram for one pixel. It has `nbins` bins in cos θ (instrument) and `nbins2` bins in θ_sun. Optionally it also has φ bins folded into 0–45°. θ_sun bins are stored **sparsely**: only θ_sun rows that were filled are allocated (`m_icostheta2toi`). **The binning is static (class-wide)**: `setBinning`, `setPhiBins`, `setPower2`. The θ_sun index is `floor(((1−cosθsun)/(1−cosmin2))^(1/power2) · nbins2)`. Integration helpers take functors of cos θ (and φ). |
| `ExposureSun` | `ExposureSun.h` | The HEALPix sky of `CosineBinner2D`s. It caches pixel directions. Its `Filler` applies the zenith cut and the 1/d² weight (only where cos θ_sun > `distCosCut`). It reads and writes the `EXPOSURESUN`-style table (`Index`/`Values` variable-length columns, plus `NBRBINS`, `NBRBINS2`, `COSMIN2`, `POWER2`, `PHIBINS`, `THETABIN`, `ZENMAX` keywords). |
| `ExposureCubeSun` | `ExposureCubeSun.h` | The livetime cube. `load(FT2 table)` walks the FT2 rows, weights each one by its GTI overlap fraction, and fills both the plain and the livetime-fraction-weighted `ExposureSun`. It writes `EXPOSURESUN`, `WEIGHTED_EXPOSURESUN`, `CTHETABOUNDS`, `THETASUNBOUNDS` and `GTI`, with keywords `SSBODY`, `DISTTCUT`, `TIMEDIST`, `TIME`. `value(dir, cosθsun, aeff, E)` applies the IRF efficiency factors, as `Likelihood::ExposureCube` does. `stringToBody` accepts any `astro::SolarSystem::Body` name, although the par file mentions only `SUN` and `MOON`. |
| `Observation` | `Observation.h` | Bundles a `Likelihood::Observation` with an `ExposureCubeSun` so that `Aeff` functors can reach both. |
| `BinnedExposureSun` | `BinnedExposureSun.h` | Exposure on a WCS grid × energy × θ_sun (4-D primary image, `ENERGIES`, `COSTHETASUN`, `GTI`). Output of `gtexpcubesun`. Nothing reads it back. |
| `HealpixExposureSun` | `HealpixExposureSun.h` | Exposure on the cube's HEALPix grid × energy × θ_sun (`EXPOSURE` table, `ENERGIES`, `THETASUN`, `GTI`, `AVGDIST`/`DISTTCUT` keywords). Output of `gtexphpsun`, input of `gtsuntemp`. `integrate(E, ra, dec, f)` sums exposure × f(cos θ_sun) over θ_sun bins. |
| `SolarProfile`, `FitsSolarProfile` | `SolarProfile.h`, `FitsSolarProfile.h` | Intensity I(cos θ, E): linear interpolation in cos θ, power-law in E (linear if either value ≤ 0). `average(cmin, cmax, E)` uses the trapezoid rule. It throws outside the tabulated range. `FitsSolarProfile` reads `ANGLES.Angle` (degrees), `ENERGIES.Energy` (MeV), `SST_PROFILE.Profile` and the `AVGDIST` keyword. |
| `SolarTemplate` | `SolarTemplate.h` | Output of `gtsuntemp`. Per pixel and energy it computes `Σ_θsun exposure(θ_sun) · ⟨I⟩(θ_sun, E) / avgexp(E, pixel)`. The disk bin(s) inside `distCosCut` are scaled by `(profile AVGDIST / cube AVGDIST)²`. The result is an intensity map cube (primary image, `ENERGIES`, `GTI`). |
| `LikeExposureSun` | `LikeExposureSun.h` | **Dead.** An older livetime-cube class. It is compiled into the library but nothing uses it. Commit `e49fb96` says it was merged into `ExposureCubeSun`. It is the only user of `data/LivetimeCubeSunTemplate`. |

## Data flow and control flow

Typical solar-template workflow (see the `fermitools-fhelp` texts):

```
FT1 (GTIs, DSS) + FT2 ──gtltcubesun──▶ ltcube_sun.fits   (EXPOSURESUN, WEIGHTED_EXPOSURESUN, CTHETABOUNDS, THETASUNBOUNDS, GTI)
                         [gtltsumsun: add several such cubes]
ltcube_sun + irfs ──gtexphpsun──▶ expsun.fits            (HEALPix × E × θ_sun exposure)
ordinary ltcube + cmap ──gtexpcube2 (Likelihood)──▶ avgexp.fits
expsun + avgexp + sunprof (+ cmap or geometry) ──gtsuntemp──▶ template.fits (intensity map cube for Likelihood)
```

`gtexpcubesun` (a WCS-grid version of `gtexphpsun`) is not part of this
chain; see the gotchas.

### `gtltcubesun` (`src/makeSolarExposureCube/makeSolarExposureCube.cxx`)

```
Prompt all; Save
readRoiCuts: Likelihood::RoiCuts::readCuts(FT1 files, evtable) → GTIs/time cuts; tmin/tmax pars override (if non-zero)
clobber=no and outfile exists → print error, std::exit(1)
createDataCube:
  choose θ_sun binning: first bin fixed at 0.25°; nbins2 = floor((1−cos θmax)/(1−cos 0.25°))^(1/powerbinsun); power2 re-derived
  ExposureCubeSun(binsz, dcostheta, 0.25, thetasunmax, power2, timeCuts, gtis, body, zmax)
  for each FT2 file: readTable(SC_DATA, filter START/STOP within [tmin−60, tmax+60]) → ExposureCubeSun::load
     per row in GTI: body direction at mid-interval from sc_position; d = |body − sc|;
     fill(scz, sun_dir, zenith, livetime·fraction, 1/d², distCosCut)  [or fill_zenith with φ if phibins>0]
writeFile(...) then writeDateKeywords on 6 HDUs; CREATOR "gtltcubesun <ver>", VERSION = file_version
```

### `gtexphpsun` / `gtexpcubesun`

Each tool builds `Likelihood::AppHelpers` (IRFs), turns off φ dependence
if `ignorephi=yes`, reads the sun cube (this resets the static
`CosineBinner2D` binning), and generates energies:
- from the `ENERGYBINS` table of `ebinfile`, converting keV to MeV, or
- logarithmically from `emin`/`emax`/`enumbins`. With `bincalc=CENTER`,
  geometric bin centres are used.

For each pixel × energy × θ_sun bin it then evaluates
`ExposureCubeSun::value`. `gtexphpsun` skips θ_sun bins that are empty
for that pixel. Finally it writes the output and copies the GTI from the
cube. `gtexpcubesun` with `cmap` set takes its geometry from the counts
map and does **not** copy the GTI.

### `gtsuntemp`

`gtsuntemp` reads `expsun` (`HealpixExposureSun`), `avgexp`
(`Likelihood::BinnedExposure`) and `sunprof`. It builds the map
geometry from `cmap` or from the parameters, then calls
`SolarTemplate::computeMap` and `writeOutput`. The energy grids of
`expsun`, `avgexp` and the output must match; otherwise the tool fails
with a "some energy bins are not found" error (per the help text).

## Public interfaces

### Applications

| Tool | Par file | Key parameters |
|------|----------|----------------|
| `gtltcubesun` | `pfiles/gtltcubesun.par` | `evfile`, `scfile`, `outfile`, `body` (`SUN`/`MOON`; any `astro` body name works), `dcostheta=0.025`, `thetasunmax=180`, `binsz=0.5`; hidden `phibins=0`, `tmin`/`tmax` (MET override if non-zero), `powerbinsun=2.7`, `zmax=180`, `file_version=1`. |
| `gtltsumsun` | `pfiles/gtltsumsun.par` | `infile1` (file or `@list`), `infile2` (default `none`), `outfile`. |
| `gtexpcubesun` | `pfiles/gtexpcubesun.par` | `infile` (sun cube), `cmap` (or `none` plus `nxpix`…`proj`), `outfile`, `irfs` (**default `P7SOURCE_V6`**; there is no `evtype` parameter), energy binning; hidden `bincalc`, `ignorephi`, `thmin`/`thmax`, `table`. |
| `gtexphpsun` | `pfiles/gtexphpsun.par` | `infile`, `outfile`, `irfs=CALDB` (see gotchas), hidden `evtype=INDEF`, energy binning, `bincalc`, `ignorephi`, `thmin`/`thmax`. `binsz` is hidden and unused: the grid is the cube's. |
| `gtsuntemp` | `pfiles/gtsuntemp.par` | `expsun`, `avgexp`, `sunprof`, `cmap` (or geometry), `outfile`, energy binning; hidden `irfs`/`evtype` ("Required by apphelpers, not used in the calculation"). |

All five tools follow the `st_app` exit-status contract
(`src/st_app/doc/architecture/README.md`): uncaught exceptions exit 1.
On Linux, `st_app` enables FPE trapping (`TRAP_FPE`), so a floating-point
fault becomes **SIGFPE (exit 136)**. That is how the `gtexpcubesun` bug
below shows up.

### C++ library

The headers in `SolarSystemTools/` are installed to
`include/SolarSystemTools`, but no other package uses them.

### Python

None.

### Environment variables

The package calls no `getenv` directly. It relies on:
- `TIMING_DIR`, for `JPLEPH.405`, read by `astro::SolarSystem` on every
  FT2 row in `gtltcubesun`;
- `CALDB` and related variables (IRFs), `PFILES`, and `FERMI_INST_DIR`.

See "Runtime environment" in ScienceTools' `docs/architecture/README.md`.
`readExposureCubeSun` expands `$VAR` in the file name
(`facilities::Util::expandEnvVar`).

### File formats

- **Sun livetime cube** (`gtltcubesun`): `EXPOSURESUN` and
  `WEIGHTED_EXPOSURESUN` hold one row per HEALPix pixel. The
  variable-length `Index` (`1PJ`) and `Values` (`1PE`) columns store only
  the filled bins. Header keywords: `NSIDE`, `ORDERING=NESTED`,
  `COORDSYS=EQU`, `NBRBINS`, `COSMIN`, `NBRBINS2`, `COSMIN2`, `POWER2`,
  `PHIBINS`, `ZENMAX`, `SSBODY`, `DISTTCUT`, `TIMEDIST`, `TIME`, DSS time
  keywords. Further HDUs: `CTHETABOUNDS`, `THETASUNBOUNDS`, `GTI`. This is
  **not** compatible with `gtltcube` or `healpix::CosineBinner` files.
- **`gtexphpsun` output:** `EXPOSURE` (HEALPix rows, NENERGIES ×
  NTHBINS values per pixel), `ENERGIES`, `THETASUN` (`Theta_Min`/
  `Theta_Max` in degrees), `GTI`.
- **`gtexpcubesun` output:** 4-D primary image (x, y, energy, θ_sun),
  `ENERGIES`, `COSTHETASUN`, `GTI`.
- **`gtsuntemp` output:** primary image (x, y, energy), `ENERGIES`, and
  `GTI` (copied only when `cmap=none`).
- **Profile input:** `ANGLES(Angle)`, `ENERGIES(Energy)`,
  `SST_PROFILE(Profile)` vector column, with keyword `AVGDIST`. Note the
  row/vector orientation question under *Invariants and gotchas*.

## Build and test

- **CMake targets:** `SolarSystemTools` (STATIC; links everything
  `PUBLIC`), `gtltcubesun` (built from `src/makeSolarExposureCube/`; that
  directory name is historical, and there is no separate
  `makeSolarExposureCube` executable), `gtltsumsun`, `gtexpcubesun`,
  `gtexphpsun`, `gtsuntemp`. `install(TARGETS …)` has `EXPORT
  fermiTargets` commented out.
- **Install:**
  - headers go to `include/SolarSystemTools`
  - `pfiles/` goes to `syspfiles`
  - **`data/` goes to `${FERMI_INSTALL_DATADIR}` itself**
    (`share/fermitools/data/LivetimeCubeSunTemplate`). That is neither
    `data/SolarSystemTools/` nor the `refdata/fermi/SolarSystemTools`
    directory that `getDataPath("SolarSystemTools")` returns. Only the
    dead `LikeExposureSun::writeFile` looks for the template, so nothing
    breaks today.
- **Tests:** none. There is no test program, `recipe/tests/ST-unit-test`
  has no entry for this package, and there is no `outref`. A smoke test
  with `src/evtbin/data/ft1tiny.fits` and `ft2tiny.fits` (`binsz=10
  dcostheta=0.1`) runs `gtltcubesun` in under a second. With the default
  `binsz=0.5` it takes about a minute even on these tiny files.
- **Legacy:** `SConscript`, `SolarSystemToolsLib.py`,
  `azure-pipelines.yml`.
- **Unmerged branches:**
  - `origin/addition_creator_fits` adds a `CREATOR` banner via
    `GenericUtils::creator_banner` and changes par files.
  - `origin/f26-cpp17-compliance` is mechanical.
  - `origin/Issue17` makes the libraries optionally static.

## Invariants and gotchas

**`gtexpcubesun` is broken for the default θ_sun range**
- `ExposureCubeSun::thetaBinsSun` returns θ_sun bin edges in **radians**
  (binning moved from cos θ_sun to θ_sun in commit `e49fb96`).
  `BinnedExposureSun::computeMap` passes the bin-centre **θ** to
  `ExposureCubeSun::value(dir, costhetasun, …)` as if it were cos θ_sun,
  and `writeOutput` writes the θ values into `CosTheta_Min`/`CosTheta_Max`.
  `HealpixExposureSun` does this correctly (`cos(...)`).
- When the cube's `thetasunmax` exceeds 1 rad (57.3°), which includes the
  default 180°, `costheta2_index` computes `pow(negative, 1/power2)`.
  That gives a NaN, which traps as **SIGFPE, exit 136** (confirmed at
  runtime with Fermitools 2.5.1 on a `thetasunmax=180` cube).
- With `thetasunmax=50` the tool runs, but the exposures are sampled at
  the wrong θ_sun bins, and the `COSTHETASUN` table holds radians
  (confirmed at runtime: first row `(0, 0.00436)`).
- Use `gtexphpsun` instead.

**IRF selection**
- **`gtexphpsun` with the default `irfs=CALDB` fails** with `File not
  found: ` (an empty name), exit 1 (confirmed at runtime). It constructs
  `AppHelpers(&pars, "none")`, so no input file is named from which to
  read DSS event-class keywords. Give an explicit IRF name and `evtype`,
  for example `irfs=P8R3_SOURCE_V3 evtype=3`, which works.
  `gtexpcubesun` with `cmap=none` has the same problem with `CALDB`.
  Its own default, `P7SOURCE_V6`, is a Pass 7 IRF.
- `gtexpcubesun` has no `evtype` parameter, so with a Pass 8 IRF name no
  event-type cut is applied (Likelihood `AppHelpers` behaviour; from
  reading the code).

**File handling**
- **`gtltsumsun` ignores `clobber`**: it always `std::remove`s the output
  first (confirmed at runtime: `clobber=no` still overwrote an existing
  file, exit 0). `gtexpcubesun`, `gtexphpsun` and `gtsuntemp` check
  through `AppHelpers::checkOutputFile`. `gtltcubesun` checks itself.
- `gtltsumsun` with a single `infile1` and the default `infile2=none`
  tries to open a file named `none` and exits 1. Pass an `@list` or two
  files.
- `gtltsumsun` checks that the `CosineBinner2D` binning matches across
  inputs, and `+=` checks that the bodies, the number of HEALPix pixels
  and the zenith cut match.

**Static binning state**
- `CosineBinner2D`'s binning is class-static, and every
  `ExposureSun(filename)` load **overwrites** it. Two cubes with
  different binning cannot coexist in one process. Comparing the static
  values before and after a load, as `gtltsumsun` does, is the only
  safeguard.

**Physics conventions**
- The first θ_sun bin is fixed at 0.25° (the solar and lunar disk).
  `powerbinsun` is adjusted so that an integer number of bins fits; the
  value actually used is written as `POWER2`.
- The distance weighting stores livetime × (1/d²) with d in
  light-seconds, inside `DISTTCUT`. That is 0.25° for the Sun, and
  everywhere (cut = −1) for the Moon. `gtsuntemp` rescales by
  `(AVGDIST_profile / AVGDIST_cube)²` only in θ_sun bins whose centre lies
  inside `distCosCut` (the code comment says "Assume that distCosCut
  aligns with one of the bins").
- The body position is computed at the mid-point of each FT2 interval,
  but from the spacecraft position at the interval **start**.
  `MJD = MET/86400 + missionStart` ignores leap seconds (about 5 s since
  2001). Both effects are negligible for 30 s FT2 bins (inferred).
- `zmax < 180` prints a warning that this is not equivalent to a
  `gtselect` zenith cut.

**Ephemeris range**
- `astro::SolarSystem` reads `$TIMING_DIR/JPLEPH.405`. The installed copy
  covers JD 2433264.5–2469808.5, that is **1949-12-14 to 2050-01-02**
  (MET ≈ 1.546 × 10⁹). `gtltcubesun` cannot process FT2 data after that
  date.

**Profile orientation**
- `FitsSolarProfile` treats **each `SST_PROFILE` row as one energy**, with
  the vector running over angles (`m_profile[row·nθ + i]`). The
  `fermitools-fhelp` text says the opposite: "each row gives the spectrum
  of the emission for a given distance". The code does not check the
  vector length, so a file with the other orientation would be read out
  of bounds (from reading the code; not checked with a real profile
  file). The comment in `SolarProfile.h` names the indices the other way
  round from the code.

**Robustness**
- The progress-dot code in `BinnedExposureSun`, `HealpixExposureSun` and
  `SolarTemplate` computes `iter % (npix/20)`. For maps with fewer than 20
  pixels this is integer division by zero (SIGFPE). This is inherited from
  the Likelihood code it was copied from (from reading the code).
- `gtltcubesun` crashed once with SIGSEGV (exit 139) after loading
  `ft2tiny.fits` with `binsz=10 dcostheta=0.1`. Six identical reruns
  succeeded. Possibly an uninitialized read; not diagnosed (observed at
  runtime, not reproduced).
- `gtexpcubesun::copyHeaderKeywords` is never called, so `gtexpcubesun`
  output lacks `TSTART`/`TSTOP`/`DATE-OBS` and the IRF DSS keywords.

## Extension points

- **Another body:** `gtltcubesun body=` already accepts any
  `astro::SolarSystem::Body` name (`stringToBody`). The 0.25° disk bin and
  the `distCosCut` logic are tuned to the Sun and Moon. Review
  `createDataCube` and the `ExposureCubeSun` constructor before using it
  for planets.
- **Another profile format:** subclass `SolarProfile` (fill
  `m_costheta` descending, `m_energies` ascending, `m_profile`,
  `m_avgDist`) and construct it in `gtsuntemp`.
- **Fixing `gtexpcubesun`:** convert θ to cos θ in
  `BinnedExposureSun::computeMap` and write θ (or cos θ) consistently in
  `writeOutput`, following `HealpixExposureSun`.
- Any change to the cube format must stay readable by `ExposureSun::load`
  (keywords listed above) and should bump `file_version`.

## Further reading

- `src/fermitools-fhelp/fhelp_files/gtltcubesun.txt`, `gtltsumsun.txt`,
  `gtexphpsun.txt`, `gtsuntemp.txt`: the user-level descriptions.
- `src/Likelihood/doc/architecture/`: `AppHelpers`, `Observation`,
  `ExposureCube`/`gtltcube` (the template for `gtltcubesun`), and
  `BinnedExposure`.
- `src/healpix/doc/architecture/` (`CosineBinner`, `HealpixArray`) and
  `src/map_tools/doc/architecture/` (`BasicExposure`).
- `src/astro/doc/architecture/`: `SolarSystem` and `TIMING_DIR`.
- E. Orlando & A. Strong 2008, A&A 480, 847, and Abdo et al. 2011, ApJ
  734, 116 (Fermi-LAT quiet Sun; the motivation for the IC halo template).
