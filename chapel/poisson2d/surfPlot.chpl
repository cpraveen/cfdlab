use IO;
import Subprocess.spawn;

proc plot(
    const ref a : [?d] real,
    lens: 2*real,
    title: string,
    varName: string
) throws {
    // create a data file
    const fileName = title + ".dat";
    var fw = openWriter(fileName, locking=false);

    // print array data
    fw.write(a);
    fw.close();

    // call the python plotting script
    spawn(["python3", "surfPlot.py", fileName, title, varName, lens[0]:string, lens[1]:string]);
}
