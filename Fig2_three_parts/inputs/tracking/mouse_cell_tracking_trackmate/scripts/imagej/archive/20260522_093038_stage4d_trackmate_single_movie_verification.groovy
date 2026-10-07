import ij.IJ
import ij.ImagePlus
import ij.ImageStack
import ij.process.ImageProcessor
import ij.process.ByteProcessor

import java.io.File
import java.io.FileWriter
import java.io.PrintWriter
import java.util.regex.Matcher
import java.util.regex.Pattern

import fiji.plugin.trackmate.Model
import fiji.plugin.trackmate.Settings
import fiji.plugin.trackmate.TrackMate
import fiji.plugin.trackmate.Logger
import fiji.plugin.trackmate.SelectionModel
import fiji.plugin.trackmate.detection.ThresholdDetectorFactory
import fiji.plugin.trackmate.tracking.jaqaman.SparseLAPTrackerFactory
import fiji.plugin.trackmate.gui.displaysettings.DisplaySettings
import fiji.plugin.trackmate.visualization.table.TrackTableView
import fiji.plugin.trackmate.visualization.hyperstack.HyperStackDisplayer

final String inputDirPath = '/data/mouse_cell_tracking/'
final String dataOutDirPath = '/app/data/projects/mouse_cell_tracking_trackmate/data/'
final String figOutDirPath = '/app/data/projects/mouse_cell_tracking_trackmate/figures/'

final int maxFramesToLoad = 50
final int fixedThreshold = 150
final int backgroundGray = 123
final int spotWhite = 255

final double linkingMaxDistance = 60.0d
final double gapClosingMaxDistance = 100.0d
final int maxFrameGap = 4
final double alternativeLinkingCostFactor = 1.05d
final boolean allowGapClosing = true
final boolean allowTrackSplitting = false
final boolean allowTrackMerging = false

File inputDir = new File(inputDirPath)
File dataOutDir = new File(dataOutDirPath)
File figOutDir = new File(figOutDirPath)
if (!dataOutDir.exists()) dataOutDir.mkdirs()
if (!figOutDir.exists()) figOutDir.mkdirs()

if (!inputDir.exists() || !inputDir.isDirectory()) {
    println('FAILURE: Input images directory missing or invalid: ' + inputDirPath)
    IJ.log('FAILURE: Input images directory missing or invalid: ' + inputDirPath)
    return
}

File[] tifFiles = inputDir.listFiles({ File d, String name ->
    String n = name.toLowerCase()
    return n.endsWith('.tif') || n.endsWith('.tiff')
} as FilenameFilter)

if (tifFiles == null || tifFiles.length == 0) {
    println('FAILURE: No TIFF frames found in ' + inputDirPath)
    IJ.log('FAILURE: No TIFF frames found in ' + inputDirPath)
    return
}

Pattern numPattern = Pattern.compile('(\\d+)')
List<File> sortedFiles = tifFiles.toList().sort { File a, File b ->
    Matcher ma = numPattern.matcher(a.getName())
    Matcher mb = numPattern.matcher(b.getName())
    long na = ma.find() ? Long.parseLong(ma.group(1)) : Long.MAX_VALUE
    long nb = mb.find() ? Long.parseLong(mb.group(1)) : Long.MAX_VALUE
    if (na != nb) return na <=> nb
    return a.getName() <=> b.getName()
}

int framesToLoad = Math.min(maxFramesToLoad, sortedFiles.size())
if (framesToLoad < 1) {
    println('FAILURE: No frames available after sorting.')
    IJ.log('FAILURE: No frames available after sorting.')
    return
}

ImagePlus firstImp = IJ.openImage(sortedFiles.get(0).getAbsolutePath())
if (firstImp == null) {
    println('FAILURE: Could not open first frame: ' + sortedFiles.get(0).getAbsolutePath())
    IJ.log('FAILURE: Could not open first frame: ' + sortedFiles.get(0).getAbsolutePath())
    return
}

int w = firstImp.getWidth()
int h = firstImp.getHeight()
ImageStack processedStack = new ImageStack(w, h)
firstImp.close()

for (int i = 0; i < framesToLoad; i++) {
    File f = sortedFiles.get(i)
    ImagePlus frameImp = IJ.openImage(f.getAbsolutePath())
    if (frameImp == null) {
        println('WARNING: Skipping unreadable frame: ' + f.getAbsolutePath())
        IJ.log('WARNING: Skipping unreadable frame: ' + f.getAbsolutePath())
        continue
    }

    ImageProcessor ip = frameImp.getProcessor().duplicate()
    ByteProcessor bp = ip.convertToByteProcessor(true)

    byte[] px = (byte[]) bp.getPixels()
    for (int p = 0; p < px.length; p++) {
        int v = px[p] & 0xff
        int bin = (v >= fixedThreshold) ? 0 : 255
        int inv = 255 - bin
        int out = (inv == 0) ? backgroundGray : spotWhite
        px[p] = (byte) (out & 0xff)
    }

    processedStack.addSlice(f.getName(), bp)
    frameImp.close()
}

if (processedStack.getSize() < 2) {
    println('FAILURE: Need at least 2 frames for tracking; got ' + processedStack.getSize())
    IJ.log('FAILURE: Need at least 2 frames for tracking; got ' + processedStack.getSize())
    return
}

ImagePlus imp = new ImagePlus('verification_sequence_preprocessed', processedStack)
def imp2 = imp.duplicate()
imp.close()

imp2.setDimensions(1, 1, imp2.getStackSize())
imp2.setOpenAsHyperStack(true)

IJ.log('Loaded frames: ' + imp2.getNFrames())
IJ.log('Preprocessing: 16-bit->8-bit, threshold=' + fixedThreshold + ', binary->invert, background=' + backgroundGray + ', spots=' + spotWhite)

Model model = new Model()
model.setLogger(Logger.IJ_LOGGER)

Settings settings = new Settings(imp2)
settings.detectorFactory = new ThresholdDetectorFactory()
settings.detectorSettings = [
    'TARGET_CHANNEL': (Integer) 1,
    'INTENSITY_THRESHOLD': 150.0d,
    'SIMPLIFY_CONTOURS': true
]
settings.initialSpotFilterValue = -1.0d
settings.addAllAnalyzers()

settings.trackerFactory = new SparseLAPTrackerFactory()
settings.trackerSettings = settings.trackerFactory.getDefaultSettings()
settings.trackerSettings['LINKING_MAX_DISTANCE'] = linkingMaxDistance
settings.trackerSettings['GAP_CLOSING_MAX_DISTANCE'] = gapClosingMaxDistance
settings.trackerSettings['MAX_FRAME_GAP'] = (Integer) maxFrameGap
settings.trackerSettings['ALTERNATIVE_LINKING_COST_FACTOR'] = alternativeLinkingCostFactor
settings.trackerSettings['ALLOW_GAP_CLOSING'] = allowGapClosing
settings.trackerSettings['ALLOW_TRACK_SPLITTING'] = allowTrackSplitting
settings.trackerSettings['ALLOW_TRACK_MERGING'] = allowTrackMerging

TrackMate trackmate = new TrackMate(model, settings)
if (!trackmate.checkInput()) {
    String err = trackmate.getErrorMessage()
    println('FAILURE: TrackMate input validation failed: ' + err)
    IJ.log('FAILURE: TrackMate input validation failed: ' + err)
    imp2.show()
    return
}
if (!trackmate.process()) {
    String err = trackmate.getErrorMessage()
    println('FAILURE: TrackMate processing failed: ' + err)
    IJ.log('FAILURE: TrackMate processing failed: ' + err)
    imp2.show()
    return
}

int trackCount = model.getTrackModel().nTracks(true)
int totalDetections = 0
for (def s in model.getSpots().iterable(false)) {
    totalDetections++
}

File summaryCsv = new File(dataOutDir, 'Verification_TrackMate_Summary.csv')
PrintWriter pw = new PrintWriter(new FileWriter(summaryCsv))
pw.println('metric,value')
pw.println('track_count,' + trackCount)
pw.println('total_detections,' + totalDetections)
pw.println('loaded_frames,' + imp2.getNFrames())
pw.println('detector_threshold,' + fixedThreshold)
pw.println('background_gray,' + backgroundGray)
pw.println('spot_gray,' + spotWhite)
pw.println('linking_max_distance_px,' + linkingMaxDistance)
pw.println('gap_closing_max_distance_px,' + gapClosingMaxDistance)
pw.println('max_frame_gap,' + maxFrameGap)
pw.println('alternative_linking_cost_factor,' + alternativeLinkingCostFactor)
pw.println('allow_gap_closing,' + allowGapClosing)
pw.println('allow_track_splitting,' + allowTrackSplitting)
pw.println('allow_track_merging,' + allowTrackMerging)
pw.close()

boolean tableExported = false
try {
    SelectionModel sm = new SelectionModel(model)
    def ds = DisplaySettings.defaultStyle()
    TrackTableView tableView = new TrackTableView(model, sm, ds)
    File spotsCsv = new File(dataOutDir, 'Verification_TrackMate_Spots.csv')
    File tracksCsv = new File(dataOutDir, 'Verification_TrackMate_Tracks.csv')
    tableView.getSpotTable().exportToCsv(spotsCsv)
    tableView.getTrackTable().exportToCsv(tracksCsv)
    tableExported = true
} catch (Exception e) {
    IJ.log('WARNING: TrackTableView CSV export failed: ' + e.getMessage())
}

try {
    SelectionModel sm2 = new SelectionModel(model)
    def ds2 = DisplaySettings.defaultStyle()
    HyperStackDisplayer displayer = new HyperStackDisplayer(model, sm2, imp2, ds2)
    displayer.render()
    displayer.refresh()
    ImagePlus qcFlat = imp2.flatten()
    File qcPng = new File(figOutDir, 'Verification_TrackMate_Overlay_QC.png')
    IJ.saveAs(qcFlat, 'PNG', qcPng.getAbsolutePath())
    qcFlat.close()
} catch (Exception e) {
    IJ.log('WARNING: QC overlay export failed: ' + e.getMessage())
}

println('SUCCESS: TrackMate verification completed.')
println('Tracks: ' + trackCount)
println('Total detections: ' + totalDetections)
println('Summary CSV: ' + summaryCsv.getAbsolutePath())
println('Tables exported via TrackTableView: ' + tableExported)
IJ.log('SUCCESS: TrackMate verification completed. Tracks=' + trackCount + ', Total detections=' + totalDetections)

imp2.show()