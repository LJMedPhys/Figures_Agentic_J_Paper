import ij.IJ
import ij.ImagePlus
import ij.ImageStack
import ij.process.ByteProcessor

import java.io.File
import java.io.FileWriter
import java.io.PrintWriter
import java.io.FilenameFilter
import java.text.SimpleDateFormat
import java.util.Date
import java.util.regex.Matcher
import java.util.regex.Pattern

import org.jgrapht.graph.DefaultWeightedEdge

import fiji.plugin.trackmate.Model
import fiji.plugin.trackmate.Settings
import fiji.plugin.trackmate.TrackMate
import fiji.plugin.trackmate.Logger
import fiji.plugin.trackmate.SelectionModel
import fiji.plugin.trackmate.Spot
import fiji.plugin.trackmate.detection.ThresholdDetectorFactory
import fiji.plugin.trackmate.tracking.jaqaman.SparseLAPTrackerFactory
import fiji.plugin.trackmate.gui.displaysettings.DisplaySettings
import fiji.plugin.trackmate.visualization.hyperstack.HyperStackDisplayer
import fiji.plugin.trackmate.io.CSVExporter

final String inputDirPath = '/data/mouse_cell_tracking/'
final String dataOutDirPath = '/app/data/projects/mouse_cell_tracking_trackmate/data/'
final String figOutDirPath = '/app/data/projects/mouse_cell_tracking_trackmate/figures/'
final String logsOutDirPath = '/app/data/projects/mouse_cell_tracking_trackmate/logs/'

final String summaryCsvName = 'TrackMate_Summary_FullMovie_dimspots.csv'
final String spotsCsvName = 'TrackMate_Spots_FullMovie_dimspots.csv'
final String tracksCsvName = 'TrackMate_Tracks_FullMovie_dimspots.csv'
final String edgesCsvName = 'TrackMate_Edges_FullMovie_dimspots.csv'
final String qcPngName = 'TrackMate_FullMovie_dimspots_Overlay_QC.png'

final int trackmateThreshold = 150
final int preThreshold = 80
final int backgroundGray = 123
final int spotWhite = 255

final boolean enableRollingBall = true
final double rollingBallRadiusPx = 40.0d
final boolean enableGaussianBlur = true
final double gaussianSigma = 1.5d
final boolean enhanceContrastNormalize = true

final double linkingMaxDistance = 60.0d
final double gapClosingMaxDistance = 100.0d
final int maxFrameGap = 4
final double alternativeLinkingCostFactor = 1.05d
final boolean allowGapClosing = true
final boolean allowTrackSplitting = false
final boolean allowTrackMerging = false

String timestamp = new SimpleDateFormat('yyyyMMdd_HHmmss').format(new Date())
File logDir = new File(logsOutDirPath)
if (!logDir.exists()) logDir.mkdirs()
File logFile = new File(logDir, 'stage4b_trackmate_full_movie_batch_dimspots_' + timestamp + '.log')
PrintWriter logWriter = new PrintWriter(new FileWriter(logFile, true))

def logMsg = { String msg ->
    String t = new SimpleDateFormat('yyyy-MM-dd HH:mm:ss').format(new Date())
    String line = '[' + t + '] ' + msg
    println(line)
    IJ.log(line)
    logWriter.println(line)
    logWriter.flush()
}

def safeDouble = { Double v ->
    return (v == null || v.isNaN() || v.isInfinite()) ? '' : v.toString()
}

ImagePlus impRaw = null
ImagePlus imp2 = null

IJ.runMacro('setBatchMode(true);')

try {
    File inputDir = new File(inputDirPath)
    File dataOutDir = new File(dataOutDirPath)
    File figOutDir = new File(figOutDirPath)
    if (!dataOutDir.exists()) dataOutDir.mkdirs()
    if (!figOutDir.exists()) figOutDir.mkdirs()

    if (!inputDir.exists() || !inputDir.isDirectory()) {
        logMsg('FAILURE: Input images directory missing or invalid: ' + inputDirPath)
        return
    }

    File[] tifFiles = inputDir.listFiles({ File d, String name ->
        String n = name.toLowerCase()
        return n.endsWith('.tif') || n.endsWith('.tiff')
    } as FilenameFilter)

    if (tifFiles == null || tifFiles.length == 0) {
        logMsg('FAILURE: No TIFF frames found in ' + inputDirPath)
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

    if (sortedFiles.size() < 2) {
        logMsg('FAILURE: Need at least 2 frames for tracking; found ' + sortedFiles.size())
        return
    }

    ImagePlus firstImp = IJ.openImage(sortedFiles.get(0).getAbsolutePath())
    if (firstImp == null) {
        logMsg('FAILURE: Could not open first frame: ' + sortedFiles.get(0).getAbsolutePath())
        return
    }
    int w = firstImp.getWidth()
    int h = firstImp.getHeight()
    firstImp.close()

    ImageStack processedStack = new ImageStack(w, h)
    int loadedCount = 0

    for (int i = 0; i < sortedFiles.size(); i++) {
        File f = sortedFiles.get(i)
        try {
            ImagePlus frameImp = IJ.openImage(f.getAbsolutePath())
            if (frameImp == null) {
                logMsg('WARNING: Skipping unreadable frame: ' + f.getAbsolutePath())
                continue
            }

            def workImp = frameImp.duplicate()
            frameImp.close()

            if (enableRollingBall) {
                IJ.run(workImp, 'Subtract Background...', 'rolling=' + rollingBallRadiusPx)
            }

            if (enableGaussianBlur) {
                IJ.run(workImp, 'Gaussian Blur...', 'sigma=' + gaussianSigma)
            }

            IJ.run(workImp, '8-bit', '')

            if (enhanceContrastNormalize) {
                IJ.run(workImp, 'Enhance Contrast...', 'saturated=0.35 normalize')
            }

            ByteProcessor bp = (ByteProcessor) workImp.getProcessor().convertToByteProcessor(true)
            byte[] px = (byte[]) bp.getPixels()
            for (int p = 0; p < px.length; p++) {
                int v = px[p] & 0xff
                int out = (v >= preThreshold) ? spotWhite : backgroundGray
                px[p] = (byte) (out & 0xff)
            }

            processedStack.addSlice(f.getName(), bp)
            loadedCount++
            workImp.close()
        } catch (Exception frameEx) {
            logMsg('WARNING: Failed processing frame ' + f.getAbsolutePath() + ' :: ' + frameEx.getMessage())
        }
    }

    if (processedStack.getSize() < 2) {
        logMsg('FAILURE: Preprocessed stack has fewer than 2 frames: ' + processedStack.getSize())
        return
    }

    impRaw = new ImagePlus('full_movie_dimspots_preprocessed', processedStack)
    imp2 = impRaw.duplicate()
    impRaw.close()

    imp2.setDimensions(1, 1, imp2.getStackSize())
    imp2.setOpenAsHyperStack(true)

    int nFrames = imp2.getNFrames()
    logMsg('INFO: Frames detected from folder: ' + sortedFiles.size())
    logMsg('INFO: Frames loaded into stack: ' + loadedCount)
    logMsg('INFO: Frames in tracking image: ' + nFrames)

    Model model = new Model()
    model.setLogger(Logger.IJ_LOGGER)

    Settings settings = new Settings(imp2)
    settings.detectorFactory = new ThresholdDetectorFactory()
    settings.detectorSettings = [
        'TARGET_CHANNEL': (Integer) 1,
        'INTENSITY_THRESHOLD': (double) trackmateThreshold,
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
        logMsg('FAILURE: TrackMate input validation failed: ' + trackmate.getErrorMessage())
        return
    }
    if (!trackmate.process()) {
        logMsg('FAILURE: TrackMate processing failed: ' + trackmate.getErrorMessage())
        return
    }

    int totalDetections = 0
    for (def s in model.getSpots().iterable(false)) totalDetections++
    int nTracks = model.getTrackModel().nTracks(true)

    File summaryCsv = new File(dataOutDir, summaryCsvName)
    PrintWriter sw = new PrintWriter(new FileWriter(summaryCsv))
    sw.println('metric,value')
    sw.println('n_frames,' + nFrames)
    sw.println('total_detections,' + totalDetections)
    sw.println('n_tracks,' + nTracks)
    sw.println('detector,ThresholdDetectorFactory')
    sw.println('detector_threshold,' + trackmateThreshold)
    sw.println('prethreshold_binary_mask,' + preThreshold)
    sw.println('background_gray,' + backgroundGray)
    sw.println('spot_gray,' + spotWhite)
    sw.println('enable_rolling_ball,' + enableRollingBall)
    sw.println('rolling_ball_radius_px,' + rollingBallRadiusPx)
    sw.println('enable_gaussian_blur,' + enableGaussianBlur)
    sw.println('gaussian_sigma,' + gaussianSigma)
    sw.println('enhance_contrast_normalize,' + enhanceContrastNormalize)
    sw.println('tracker,SparseLAPTrackerFactory')
    sw.println('linking_max_distance_px,' + linkingMaxDistance)
    sw.println('gap_closing_max_distance_px,' + gapClosingMaxDistance)
    sw.println('max_frame_gap,' + maxFrameGap)
    sw.println('alternative_linking_cost_factor,' + alternativeLinkingCostFactor)
    sw.println('allow_gap_closing,' + allowGapClosing)
    sw.println('allow_track_splitting,' + allowTrackSplitting)
    sw.println('allow_track_merging,' + allowTrackMerging)
    sw.close()

    File spotsCsv = new File(dataOutDir, spotsCsvName)
    CSVExporter.exportSpots(spotsCsv.getAbsolutePath(), model, true)

    File tracksCsv = new File(dataOutDir, tracksCsvName)
    PrintWriter tw = new PrintWriter(new FileWriter(tracksCsv))
    tw.println('track_id,n_spots,n_edges,track_duration,track_start,track_stop,track_displacement,track_mean_speed')
    def fm = model.getFeatureModel()
    for (Integer tid : model.getTrackModel().trackIDs(true)) {
        int ns = model.getTrackModel().trackSpots(tid).size()
        int ne = model.getTrackModel().trackEdges(tid).size()
        Double dur = fm.getTrackFeature(tid, 'TRACK_DURATION')
        Double sta = fm.getTrackFeature(tid, 'TRACK_START')
        Double sto = fm.getTrackFeature(tid, 'TRACK_STOP')
        Double disp = fm.getTrackFeature(tid, 'TRACK_DISPLACEMENT')
        Double ms = fm.getTrackFeature(tid, 'TRACK_MEAN_SPEED')
        tw.println(tid + ',' + ns + ',' + ne + ',' + safeDouble(dur) + ',' + safeDouble(sta) + ',' + safeDouble(sto) + ',' + safeDouble(disp) + ',' + safeDouble(ms))
    }
    tw.close()

    File edgesCsv = new File(dataOutDir, edgesCsvName)
    PrintWriter ew = new PrintWriter(new FileWriter(edgesCsv))
    ew.println('track_id,source_spot_id,target_spot_id,source_frame,target_frame,source_x,source_y,target_x,target_y,edge_weight')
    def tm = model.getTrackModel()
    for (Integer tid : tm.trackIDs(true)) {
        Set<DefaultWeightedEdge> edges = tm.trackEdges(tid)
        for (DefaultWeightedEdge e : edges) {
            Spot src = tm.getEdgeSource(e)
            Spot tgt = tm.getEdgeTarget(e)
            Double sf = src.getFeature(Spot.FRAME)
            Double tf = tgt.getFeature(Spot.FRAME)
            Double sx = src.getFeature(Spot.POSITION_X)
            Double sy = src.getFeature(Spot.POSITION_Y)
            Double tx = tgt.getFeature(Spot.POSITION_X)
            Double ty = tgt.getFeature(Spot.POSITION_Y)
            double wt = tm.getEdgeWeight(e)
            ew.println(tid + ',' + src.ID() + ',' + tgt.ID() + ',' + safeDouble(sf) + ',' + safeDouble(tf) + ',' + safeDouble(sx) + ',' + safeDouble(sy) + ',' + safeDouble(tx) + ',' + safeDouble(ty) + ',' + wt)
        }
    }
    ew.close()

    try {
        SelectionModel sm = new SelectionModel(model)
        DisplaySettings ds = DisplaySettings.defaultStyle()
        HyperStackDisplayer displayer = new HyperStackDisplayer(model, sm, imp2, ds)
        displayer.render()
        displayer.refresh()

        int midFrame = Math.max(1, (int)Math.round(nFrames / 2.0d))
        imp2.setT(midFrame)
        File qcPng = new File(figOutDir, qcPngName)
        ImagePlus qcFlat = imp2.flatten()
        IJ.saveAs(qcFlat, 'PNG', qcPng.getAbsolutePath())
        qcFlat.close()
    } catch (Exception overlayEx) {
        logMsg('WARNING: QC overlay export failed: ' + overlayEx.getMessage())
    }

    logMsg('SUCCESS: TrackMate full-movie dimspots batch completed.')
    logMsg('SUCCESS: total_detections=' + totalDetections + ', n_tracks=' + nTracks + ', n_frames=' + nFrames)
    logMsg('SUCCESS: Summary=' + new File(dataOutDir, summaryCsvName).getAbsolutePath())
    logMsg('SUCCESS: Spots=' + new File(dataOutDir, spotsCsvName).getAbsolutePath())
    logMsg('SUCCESS: Tracks=' + new File(dataOutDir, tracksCsvName).getAbsolutePath())
    logMsg('SUCCESS: Edges=' + new File(dataOutDir, edgesCsvName).getAbsolutePath())
    logMsg('SUCCESS: QC PNG=' + new File(figOutDir, qcPngName).getAbsolutePath())

} catch (Exception e) {
    logMsg('FAILURE: Unhandled exception: ' + e.getMessage())
    StringWriter sw = new StringWriter()
    e.printStackTrace(new PrintWriter(sw))
    logMsg(sw.toString())
} finally {
    try {
        if (imp2 != null) imp2.close()
    } catch (Exception ignored) {}
    try {
        if (impRaw != null) impRaw.close()
    } catch (Exception ignored) {}
    try {
        IJ.run('Close All')
    } catch (Exception ignored) {}
    IJ.runMacro('setBatchMode(false);')
    logMsg('INFO: Batch mode disabled. Script finished.')
    logWriter.close()
}
