proc simplePlot(
    a: [?d] real,
    minY: real = 0, maxY: real = 1, numY: int = 8,
    minX: real = 0, maxX: real = 1
) where d.rank == 1 {
    const deltaY = (maxY - minY) / numY;
    const steps = [i in -1..(numY+1)] i * deltaY + minY;

    var yLevels : [d] int;
    for (val, vi) in zip(a, d) do
        for i in 0..numY do
            if val >= steps[i-1] && val <= steps[i] then yLevels[vi] = i;

    for si in steps.domain by -1 {
        writef("%{##.###}-|", steps[si]);
        for l in yLevels {
            if l == si { write("."); } else { write(" "); }
        }
        writeln();
    }

    writef("       |");
    for l in 0..#d.dim(0).size do write("-");
    writeln();

    writef(" x: %{##.###}", minX);
    for l in 0..#d.dim(0).size-5 do write(" ");
    writef("%{##.###}", maxX);
    writeln("\n");
}
