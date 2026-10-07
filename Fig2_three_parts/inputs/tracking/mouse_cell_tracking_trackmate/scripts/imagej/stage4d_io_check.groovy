import ij.IJ
import ij.ImagePlus
import java.io.File
import java.io.PrintWriter

String inputDirPath = '/data/mouse_cell_tracking/'
String logPath = '/app/data/projects/mouse_cell_tracking_trackmate/logs/IO_Check.txt'

File inputDir = new File(inputDirPath)
File logFile = new File(logPath)
if (!logFile.getParentFile().exists()) {
    logFile.getParentFile().mkdirs()
}

PrintWriter log = null
try {
    log = new PrintWriter(logFile)
    log.println('IO CHECK REPORT')
    log.println('Input directory: ' + inputDirPath)
    log.println('')

    if (!inputDir.exists() || !inputDir.isDirectory()) {
        String msg = 'FAIL: Input directory missing or not a directory.'
        log.println(msg)
        println(msg)
        return
    }

    File[] tifFiles = inputDir.listFiles({ File d, String name ->
        String n = name.toLowerCase()
        return n.endsWith('.tif') || n.endsWith('.tiff')
    } as FilenameFilter)

    if (tifFiles == null || tifFiles.length == 0) {
        String msg = 'FAIL: No .tif/.tiff files found.'
        log.println(msg)
        println(msg)
        return
    }

    Arrays.sort(tifFiles, { a, b -> a.getName() <=> b.getName() } as Comparator<File>)

    int total = tifFiles.length
    int openedOk = 0
    int failed = 0

    String firstName = tifFiles[0].getName()
    String lastName = tifFiles[total - 1].getName()

    int firstW = -1
    int firstH = -1
    int firstBit = -1

    int lastW = -1
    int lastH = -1
    int lastBit = -1

    // Verify read access for all TIFF files.
    IJ.runMacro('setBatchMode(true);')
    for (int i = 0; i < tifFiles.length; i++) {
        File f = tifFiles[i]
        try {
            ImagePlus imp = IJ.openImage(f.getAbsolutePath())
            if (imp == null) {
                failed++
                log.println('OPEN_FAIL: ' + f.getName())
                continue
            }

            openedOk++

            if (i == 0) {
                firstW = imp.getWidth()
                firstH = imp.getHeight()
                firstBit = imp.getBitDepth()
            }
            if (i == tifFiles.length - 1) {
                lastW = imp.getWidth()
                lastH = imp.getHeight()
                lastBit = imp.getBitDepth()
            }

            imp.close()
        } catch (Exception e) {
            failed++
            log.println('OPEN_EXCEPTION: ' + f.getName() + ' :: ' + e.getMessage())
        }
    }
    IJ.run('Close All')
    IJ.runMacro('setBatchMode(false);')

    log.println('Total TIFF files: ' + total)
    log.println('Opened successfully: ' + openedOk)
    log.println('Failed to open: ' + failed)
    log.println('')
    log.println('Example first file: ' + firstName)
    log.println('First dimensions: ' + firstW + ' x ' + firstH + ', bit depth: ' + firstBit)
    log.println('Example last file: ' + lastName)
    log.println('Last dimensions: ' + lastW + ' x ' + lastH + ', bit depth: ' + lastBit)

    if (failed == 0) {
        println('SUCCESS: IO check complete. All TIFF files opened.')
        log.println('STATUS: SUCCESS')
    } else {
        println('FAIL: IO check complete with open failures. See log for details.')
        log.println('STATUS: FAIL')
    }

} finally {
    if (log != null) {
        log.flush()
        log.close()
    }
}
