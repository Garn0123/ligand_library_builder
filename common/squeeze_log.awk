# squeeze_log.awk -- collapse repeated RDKit warnings in db2_converter output.
#
#   build_ligand ... 2>&1 | awk -f common/squeeze_log.awk
#   ... | llb squeeze
#
# db2_converter's chemistry check prints "[hh:mm:ss] WARNING: Proton(s)
# added/removed" (plus a blank line) once per conformer -- thousands per
# molecule. Each unbroken run becomes one line with a count; every other line
# passes through unchanged, in order. Add patterns to NOISE to squeeze more.

BEGIN { NOISE = "WARNING: Proton\\(s\\) added/removed" }

function flush() {
    if (n > 0) {
        printf "[squeezed] %s x %d\n", last, n
        n = 0
    }
}

$0 ~ NOISE {
    sub(/^\[[0-9:]+\] /, "")      # drop the per-line timestamp
    last = $0
    n++
    next
}
/^[[:space:]]*$/ && n > 0 { next } # blank lines inside a run
{ flush(); print; fflush() }
END { flush() }
