# Contributing to OpenPIC

## Code style

All C++ source in `src/` is formatted with clang-format. Apply it before committing:

```powershell
$cf = "C:/Program Files/Microsoft Visual Studio/18/Insiders/VC/Tools/Llvm/x64/bin/clang-format.exe"
Get-ChildItem src -Include *.cpp,*.h -Recurse | ForEach-Object { & $cf -i $_.FullName }
```

The configuration is in [`.clang-format`](.clang-format) at the repo root (LLVM base, 4-space indent, 120-column limit, Allman braces for functions and structs).

Key conventions:

- Left-aligned pointers and references: `Grid& grid`, `const Particle* p`
- No `auto` for non-trivial types where the type aids readability
- Comments explain **why**, not what — well-named identifiers carry the what
- No trailing summaries or task references in comments

## Units

All physical quantities use CGS-Gaussian units. Normalizations are set in `main.lua` per task. Do not mix SI into new code without an explicit conversion layer.

## Commit messages

Follow the existing log style:

```
<type>: <short imperative summary>

Body explaining the motivation if not obvious.
```

Types: `feat`, `fix`, `perf`, `refactor`, `style`, `docs`, `test`.

Keep the subject line under 72 characters. Reference the physics motivation in the body when the change affects numerical results.

## MPI builds

Test both the standard and MPI builds before submitting changes that touch `src/`:

```powershell
cmake -B build    -DOPENPIC_ENABLE_MPI=OFF && cmake --build build    --config Release
cmake -B build-mpi -DOPENPIC_ENABLE_MPI=ON  && cmake --build build-mpi --config Release
```

A single-rank MPI run (`mpiexec -n 1 open-pic.exe -mpi main.lua`) must produce output identical to the non-MPI run for the same task.

## Third-party libraries

All dependencies are vendored under `3rdparty/`. Do not add new external dependencies without discussion. The vendored versions are:

| Library | Version |
|---------|---------|
| Lua | 5.4.8 |
| sol2 | 3.3.0 |
| fmt | 10.2.1 |
| pocketfft | — |
