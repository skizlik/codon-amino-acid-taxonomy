"""Check committed notebooks: size, errors, stderr output and execution order.

    python scripts/check_notebooks.py notebooks/*.ipynb

Fails if a notebook is over 500 KB, contains an error or any stderr output
(warnings, log lines), or was not executed top to bottom in one fresh kernel
(code cells must be numbered 1, 2, 3, ... in order). A notebook whose outputs
are all cleared passes the order check.
"""
import json
import sys
from pathlib import Path

MAX_KB = 500


def problems(path):
    out = []
    kb = path.stat().st_size / 1024
    if kb > MAX_KB:
        out.append(f"{kb:.0f} KB > {MAX_KB} KB")
    cells = json.loads(path.read_text(encoding="utf-8"))["cells"]
    code = [c for c in cells if c["cell_type"] == "code"]
    counts = [c.get("execution_count") for c in code]
    if any(n is not None for n in counts) and counts != list(range(1, len(code) + 1)):
        out.append(f"execution counts are not 1..{len(code)} in order: {counts}")
    for i, c in enumerate(cells):
        for o in c.get("outputs", []):
            if o["output_type"] == "error":
                out.append(f"cell {i}: error {o['ename']}")
            if o["output_type"] == "stream" and o["name"] == "stderr":
                out.append(f"cell {i}: stderr output: {''.join(o['text'])[:80]!r}")
        if c["cell_type"] == "code" and not "".join(c["source"]).strip():
            out.append(f"cell {i}: empty code cell")
    return out


def main(paths):
    bad = 0
    for p in map(Path, paths):
        found = problems(p)
        status = "ok" if not found else "FAIL"
        print(f"{status:4s} {p} ({p.stat().st_size / 1024:.0f} KB)")
        for msg in found:
            print(f"     {msg}")
        bad += bool(found)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
