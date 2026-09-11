# CLAUDE.md — FDA_FEP_public session rules

- Work in narrow explicit slices only.
- No agents/advisor unless explicitly authorized.
- No Python.
- Use only the pinned R 4.4.2 / Pmetrics 3.0.9 environment.
- Never refit or overwrite Pmetrics runs unless explicitly authorized.
- Never rerun simulation unless explicitly authorized.
- Unexpected command failure ends the slice: STOP, no retries or changed invocation unless authorized.
- No temp validation scripts.
- Inspect exact returned objects before interpreting them.
- One diagnostic authorizes one diagnostic only.
- Do not modify manuscript-facing analysis while auditing.
- Generated/private artifacts remain untracked.
- STOP means STOP.
