import ij.IJ
import ij.ImagePlus
import ij.ImageStack
import ij.process.ImageProcessor
import ij.process.ByteProcessor

import java.io.File
import java.io.FileWriter
import java.io.PrintWriter
import java.io.FilenameFilter
import java.text.SimpleDateFormat
import java.util.Date
import java.util.regex.Matcher
import java.util.regex.Pattern

import fiji.plugin.trackmate.Model
import fiji.plugin.trackmate.Settings
import fiji.plugin.trackmate.TrackMate
import fiji.plugin.trackmate.Logger
import fiji.plugin.trackmate.Spot
import fiji.plugin.trackmate.detection.ThresholdDetectorFactory
import fiji.plugin.trackmate.tracking.jaqaman.SparseLAPTrackerFactory

final String inputDirPath = '/data/mouse_cell_tracking/'
final String dataOutDirPath = '/app/data/projects/mouse_cell_tracking_trackmate/data/'
final String logsOutDirPath = '/app/data/projects/mouse_cell_tracking_trackmate/logs/'

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
File logsOutDir = new File(logsOutDirPath)
if (!dataOutDir.exists()) dataOutDir.mkdirs()
if (!logsOutDir.exists()) logsOutDir.mkdirs()

String timestamp = new SimpleDateFormat('yyyyMMdd_HHmmss').format(new Date())
File logFile = new File(logsOutDir, 'stage4b_trackmate_full_movie_batch_' + timestamp + '.log')
PrintWriter logWriter = new PrintWriter(new FileWriter(logFile, true))

def logMsg = { String msg ->
    println(msg)
    IJ.log(msg)
    logWriter.println(new SimpleDateFormat('yyyy-MM-dd HH:mm:ss').format(new Date()) + ' | ' + msg)
    logWriter.flush()
}

ImagePlus impOriginal = null
ImagePlus imp2 = null

IJ.runMacro('setBatchMode(true);')

try {
    logMsg('START: Full-movie TrackMate batch processing.')

    if (!inputDir.exists() || !inputDir.isDirectory()) {
        throw new RuntimeException('Input images directory missing or invalid: ' + inputDirPath)
    }

    File[] tifFiles = inputDir.listFiles({ File d, String name ->
        String n = name.toLowerCase()
        return n.endsWith('.tif') || n.endsWith('.tiff')
    } as FilenameFilter)

    if (tifFiles == null || tifFiles.length == 0) {
        throw new RuntimeException('No TIFF files found in input directory: ' + inputDirPath)
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

    logMsg('Found TIFF files: ' + sortedFiles.size())

    ImagePlus firstImp = IJ.openImage(sortedFiles.get(0).getAbsolutePath())
    if (firstImp == null) {
        throw new RuntimeException('Could not open first frame: ' + sortedFiles.get(0).getAbsolutePath())
    }

    int w = firstImp.getWidth()
    int h = firstImp.getHeight()
    firstImp.close()

    ImageStack processedStack = new ImageStack(w, h)

    int loadedCount = 0
    int skippedCount = 0

    for (int i = 0; i < sortedFiles.size(); i++) {
        File f = sortedFiles.get(i)
        try {
            ImagePlus frameImp = IJ.openImage(f.getAbsolutePath())
            if (frameImp == null) {
                skippedCount++
                logMsg('WARNING: Unreadable frame skipped: ' + f.getAbsolutePath())
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
            loadedCount++
            frameImp.close()
        } catch (Exception fe) {
            skippedCount++
            logMsg('WARNING: Failed processing frame ' + f.getAbsolutePath() + ' | ' + fe.getMessage())
        }
    }

    if (processedStack.getSize() < 2) {
        throw new RuntimeException('Need at least 2 processed frames for tracking; got ' + processedStack.getSize())
    }

    impOriginal = new ImagePlus('full_movie_preprocessed_stack', processedStack)
    imp2 = impOriginal.duplicate()
    impOriginal.close()

    imp2.setDimensions(1, 1, imp2.getStackSize())
    imp2.setOpenAsHyperStack(true)

    logMsg('Preprocessing complete: loaded=' + loadedCount + ', skipped=' + skippedCount + ', frames in stack=' + imp2.getNFrames())
    logMsg('Preprocessing params: 16-bit->8-bit, threshold=' + fixedThreshold + ', invert=true, background=' + backgroundGray + ', spots=' + spotWhite)

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
        throw new RuntimeException('TrackMate input check failed: ' + trackmate.getErrorMessage())
    }
    if (!trackmate.process()) {
        throw new RuntimeException('TrackMate processing failed: ' + trackmate.getErrorMessage())
    }

    int trackCount = model.getTrackModel().nTracks(true)
    int totalDetections = 0
    for (Spot s : model.getSpots().iterable(false)) {
        totalDetections++
    }

    File summaryCsv = new File(dataOutDir, 'TrackMate_Summary_FullMovie.csv')
    PrintWriter sw = new PrintWriter(new FileWriter(summaryCsv))
    sw.println('metric,value')
    sw.println('track_count,' + trackCount)
    sw.println('total_detections,' + totalDetections)
    sw.println('input_files_found,' + sortedFiles.size())
    sw.println('frames_loaded,' + loadedCount)
    sw.println('frames_skipped,' + skippedCount)
    sw.println('detector_threshold,' + fixedThreshold)
    sw.println('background_gray,' + backgroundGray)
    sw.println('spot_gray,' + spotWhite)
    sw.println('linking_max_distance_px,' + linkingMaxDistance)
    sw.println('gap_closing_max_distance_px,' + gapClosingMaxDistance)
    sw.println('max_frame_gap,' + maxFrameGap)
    sw.println('alternative_linking_cost_factor,' + alternativeLinkingCostFactor)
    sw.println('allow_gap_closing,' + allowGapClosing)
    sw.println('allow_track_splitting,' + allowTrackSplitting)
    sw.println('allow_track_merging,' + allowTrackMerging)
    sw.close()

    File spotsCsv = new File(dataOutDir, 'TrackMate_Spots_FullMovie.csv')
    PrintWriter spw = new PrintWriter(new FileWriter(spotsCsv))
    spw.println('spot_id,track_id,frame,x,y,quality')
    for (Spot s : model.getSpots().iterable(false)) {
        Integer tid = model.getTrackModel().trackIDOf(s)
        Double frame = s.getFeature('FRAME')
        Double x = s.getFeature('POSITION_X')
        Double y = s.getFeature('POSITION_Y')
        Double q = s.getFeature('QUALITY')
        spw.println(s.ID() + ',' + (tid == null ? '' : tid) + ',' + (frame == null ? '' : frame.intValue()) + ',' + (x == null ? '' : x) + ',' + (y == null ? '' : y) + ',' + (q == null ? '' : q))
    }
    spw.close()

    File tracksCsv = new File(dataOutDir, 'TrackMate_Tracks_FullMovie.csv')
    PrintWriter tw = new PrintWriter(new FileWriter(tracksCsv))
    tw.println('track_id,n_spots')
    for (Integer tid : model.getTrackModel().trackIDs(true)) {
        int nSpotsInTrack = model.getTrackModel().trackSpots(tid).size()
        tw.println(tid + ',' + nSpotsInTrack)
    }
    tw.close()

    File edgesCsv = new File(dataOutDir, 'TrackMate_Edges_FullMovie.csv')
    PrintWriter ew = new PrintWriter(new FileWriter(edgesCsv))
    ew.println('source_spot_id,target_spot_id,track_id')
    for (def edge : model.getTrackModel().edgeSet()) {
        Spot source = model.getTrackModel().getEdgeSource(edge)
        Spot target = model.getTrackModel().getEdgeTarget(edge)
        Integer tid = model.getTrackModel().trackIDOf(edge)
        ew.println(source.ID() + ',' + target.ID() + ',' + (tid == null ? '' : tid))
    }
    ew.close()

    logMsg('SUCCESS: TrackMate full-movie batch completed.')
    logMsg('Outputs:')
    logMsg(' - ' + summaryCsv.getAbsolutePath())
    logMsg(' - ' + spotsCsv.getAbsolutePath())
    logMsg(' - ' + tracksCsv.getAbsolutePath())
    logMsg(' - ' + edgesCsv.getAbsolutePath())
    logMsg('Tracks=' + trackCount + ' | TotalDetections=' + totalDetections)

} catch (Exception e) {
    logMsg('FAILURE: ' + e.getMessage())
    e.printStackTrace()
} finally {
    try {
        if (imp2 != null) imp2.close()
    } catch (Exception ignored) {}
    try {
        if (impOriginal != null) impOriginal.close()
    } catch (Exception ignored) {}
    try {
        IJ.run('Close All')
    } catch (Exception ignored) {}
    IJ.runMacro('setBatchMode(false);')
    logMsg('END: Batch mode off; images closed.')
    logWriter.close()
}
