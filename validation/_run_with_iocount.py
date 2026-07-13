"""Run a command and report its wall time + OS I/O counters (M0 probe for the read-once work).

Usage:  python validation/_run_with_iocount.py LABEL -- cmd arg1 arg2 ...
Prints one summary line:  IOCOUNT label=<label> rc=<rc> elapsed_s=<s> read_gb=<GB> write_gb=<GB>
(The child is polled twice a second; the final sample can miss the last <0.5 s of I/O.)
"""
import subprocess
import sys
import time

import psutil


def main():
    sep = sys.argv.index("--")
    label, cmd = sys.argv[1], sys.argv[sep + 1:]
    t0 = time.monotonic()
    child = subprocess.Popen(cmd)
    proc = psutil.Process(child.pid)
    last = None
    while child.poll() is None:
        try:
            last = proc.io_counters()
        except psutil.Error:
            break
        time.sleep(0.5)
    elapsed = time.monotonic() - t0
    rb = last.read_bytes / 1e9 if last else float("nan")
    wb = last.write_bytes / 1e9 if last else float("nan")
    print(f"IOCOUNT label={label} rc={child.returncode} elapsed_s={elapsed:.1f} "
          f"read_gb={rb:.3f} write_gb={wb:.3f}", flush=True)
    sys.exit(child.returncode or 0)


if __name__ == "__main__":
    main()
