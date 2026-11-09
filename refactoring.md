# Estru3D refactoring proposal

Task 7 of `revamp.txt`. Nothing here has been done yet; it is a proposal, in the order I would do it.

The goal, as stated: bring a 28-year-old program up to date and connect it to sumi-designer. Estru3D stays fully
parametric and does no visual design. The solver's results must not change.

## Where the code stands

| | |
|---|---|
| Own code | about 31,000 lines of Gambas, 34 forms |
| Three modules hold most of it | `grf.module` 5,485 lines, `e3d.module` 4,158, `clc.module` 3,743; `FMain.class` 2,467 |
| Global state | `e3d` 203 public symbols, `grf` 122, `FMain` 209 |
| Solver tied to the window | `clc.module` has 54 lines calling `FMain`, `Message` or `DoEvents`; `e3d.module` has 70 |
| Error handling | 115 `Try`, 3 `Catch`; 10 leftover `Stop` statements |
| Legacy OpenGL still in use | `grf.module` 69 lines, `Scene3d.class` 54 |
| Links to sumi-e | 118: 5 in `.src/Library`, `glCamera`, 111 in `.src/Sumi`, and the shaders folder |
| Tests | none in the repo |

The compiler does not check members, so calls to routines that no longer exist only fail when that line runs.

## 1. Stop the silent breakage (small, first)

**Calls to routines that do not exist.** Each of these crashes when reached:

| Call | Where | Effect |
|---|---|---|
| `e3d.addLosa`, `e3d.ReLoadElement` | `formPlano` 122, 261, 263 | the form cannot create a slab |
| `e3d.EliminarBarra` | `FMain` 997 | deleting a bar |
| `e3d.PurgarApoyos`, `e3d.purgar_cargas` | `formGenerarPortico` 96, 120 | the frame wizard |
| `clc.IniciarCalculosCompletos` | `FMain` 1894 | the "complete calculation" menu |
| `clc.IniciarCalculosEnvolventes` | `formEnvolventes` 33 | envelopes |
| `e3d.swap_nudo`, `clc.GSLMatRigEst` | `clc` 2528, 2639; `matrix8` 355 on | solver paths |

Also worth checking: `e3d.trenCargas`, `trenDistancia`, `Camino`, `envolventes`, `dezplazamientos`, `elastica`,
`elasticaXcargas`, `EsfExtBar`. They look like globals that moved into `e3d.model` in 2025. My scan is by name, so
each needs a look before it is called a bug.

For each one: restore it, or remove the menu item that reaches it. The 2025 code in git history may have some.

**Remove the junk from git.** `.hidden/imagenes` duplicates `imagenes` (44 files). `.info`, `.list`, `.settings`,
`.startup` and `.icon.png` are in `.gitignore` but still tracked, so every compile dirties the tree.

**Add a test before anything else moves.** Estru3D has no equivalent of sumi-e's `--smoke-test`. Open the examples and
an export, check counts and loads, pick a bar. The scratch harness used for this revamp is most of it already.

## 2. Separate the solver from the window (highest value)

`clc` should run with no `FMain`, no `Message` and no `DoEvents`. Progress and errors go through a small interface and
`AppLog`.

That gives three things at once:

- **Golden tests.** Calculate the benchmark examples and compare with stored results, so every later refactor can
  prove the numbers did not move. `ejemplos-benchmark.ce5` and `test-APP93.ce5` exist for this, but they are format
  5.0 and nothing opens those any more (`cargar5` has no callers). Reconnect it or convert them.
- **A `--calc file` mode**, which is what a real connection to sumi-designer needs.
- A solver that can later be replaced or moved to a thread without touching the interface.

No change to the formulas in this step. One thing to check first, by someone who knows the solver: in `GenMatPorEsp`
(`clc.module` 1233 to 1335) the test `(restriccion And 1) = 1` comes before `restriccion = 3`, so a bar pinned at both
ends takes the "pinned at start only" branch and the `= 3` branch is never reached. The plane-frame version compares
with `= 1`, `= 2`, `= 3` and does not have this problem.

## 3. Split `e3d.module`

Into modules that each own one thing:

| New | Takes from `e3d` |
|---|---|
| `FileCe5` | `cargar5`, `cargar6`, `grabar`: the text format, import only once JSON exists |
| `FileJson` | new, see below |
| `Materials` | `mats_load2`, `GetMaterialBy...`, units and factors |
| `Loads` | `ReLoadLosa`, `ReLoadViga`, `GenerarPesoPropio`, `GenerarSobrecargas` |
| `Checks` | `chequeoestructura`, `barras_coincidentes` |
| `e3d` | only what is really global: the current model and the flags |

`GetMaterialByName` and `GetMaterialByHandle` still return the first material when nothing matches. The file reader no
longer uses them that way, but `FormBarras` and `formGenerarPlantaTipo` do. Make them return `Null`.

The two slab-load routines should become one. `GenerarSobrecargas` loses the edge strip over beams (about 10% on a
6 x 4 m slab), and `ReLoadLosa` keeps the total but piles it on the interior bars. Use sumi-designer's rule: tributary
width, edge bars at half, total exactly load x area.

## 4. JSON file format

As you decided. I would copy sumi-e's scheme: `"format"`, `"version"`, a step-by-step `Upgrade`, written to a `.part`
file and renamed so a reader never gets half a file.

- Materials, springs and sections live in the file, always. The library only fills the forms.
- Nodes, bars, slabs, beams, columns, loads by state, combinations. Handles as strings.
- Member rotation as a real field, which ends the 0-or-90 limit once the solver honours it.
- `.ce5` stays as import, so old projects keep opening.
- sumi-designer writes the same JSON, replacing `Estru3DExport`.

## 5. Finish moving the drawing to `glx`

The main view now uses sumi-e's camera, matrices and shader. Still on the old fixed pipeline: node and bar numbers,
axes, the colour scale, the filled "render" mode and the envelopes.

Split `grf.module` by job: `GrfLimits`, `GrfStructure` (one builder per VBO), `GrfResults`, `GrfOverlay`. Then drop
`MatricesFijas`.

**Remove the Render tab, `Scene3d`, `GlCameraLegacy`, the STL writing and `VirtualReality/` (4.7 MB).** They are visual
output, which task 4 rules out, and the STL save is already off, so the tab shows nothing. I kept all of it for now
because deleting a feature should be your call.

## 6. Depend on less of sumi-e

`glx` and `VBO` pull in `Texture`, `Entity` and from there nearly all of sumi-e: 108 classes, forms included, to draw
lines. A change in sumi-e can break Estru3D, and often only at runtime.

Proposal: move what is truly shared (`glx`, `VBO`, `EBO`, `ShaderMatrices`, `glCamera`, `puntos`, a `PickBuffer` not
tied to `Entity`, `AppLog`) into one shared folder or Gambas library used by sumi-e, sumi-designer and Estru3D. The PBR
and texture parts stay in sumi-e behind a small interface. That is a change to sumi-e, so it needs your go-ahead.

`.src/Library/PickE3D.class` could then be retired in favour of the shared picker.

## 7. The interface

- `FMain` has 209 public symbols. Move the menu dispatcher and the view handling out of it.
- Fix what shows up in use: `<center>` printed literally in `formApoyos`, and support pictures clipped on their buttons.
- Convert the remaining bitmaps: the section and load diagrams in `imagenes` and the unused icons in `iconos`. The 40
  icons actually on screen are done. Their PNGs can go once you are happy with them.
- Move the old `Settings` file and `Config` to one JSON config, as sumi-e does.

## 8. Connecting to sumi-designer

- Reload when the exported file changes, the way sumi-e follows its `.glb`.
- Bug in `Estru3DExport.module:202`: a bar released at both ends is written as `4` (tensor) and should be `3`.
- Show sumi-designer's export warnings when opening the file.

## Suggested order

1. Section 1: broken calls, junk, smoke test.
2. Section 2: solver apart, golden tests.
3. Sections 3 and 4 together: the file code moves anyway.
4. Section 8: the live connection on top of JSON.
5. Sections 5, 6 and 7: drawing, sumi-e dependency, interface.

Each step leaves a working program, and from step 2 on the golden tests guard the numbers.

## Decisions for you

- Remove the Render tab, `Scene3d`, STL and `VirtualReality/`? (section 5)
- Create the shared library, which touches sumi-e? (section 6)
- Reconnect `cargar5`, or convert the old examples once and drop it? (section 2)
- Is the `restriccion` branch order a bug? (section 2)
