package fr.ens.biologie.genomique.toullig.fast5tofastq;

import fr.ens.biologie.genomique.eoulsan.util.LocalReporter;
import ncsa.hdf.hdf5lib.exceptions.HDF5Exception;
import org.apache.commons.compress.compressors.bzip2.BZip2CompressorOutputStream;

import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.GZIPOutputStream;

/**
 * This class allow to process all the FAST5 files of a directory.
 */
public class DirectoryProcessor {

  private final boolean saveComplementSequence;
  private final boolean saveTemplateSequence;
  private final boolean saveConsensusSequence;
  private final boolean saveTranscriptSequence;

  private boolean metrichor = false;

  private final boolean saveCompressGZIP;
  private final boolean saveCompressBZIP2;
  private final File repertoryFastqOutput;

  private final List<File> listCorruptFast5Files = new ArrayList<>();

  private final Fast5.Version version;
  private final Fast5.Type type;
  private final Fast5.Status status;
  private final Fast5.ChemistryVersion chemistryVersion;
  private final Fast5.Basecaller basecaller;

  /**
   * Constructor.
   * @param repertoryFastqOutput FASTQ output directory
   * @param saveComplementSequence save complement sequence
   * @param saveTemplateSequence save template sequence
   * @param saveConsensusSequence save consensus sequence
   * @param saveTranscriptSequence save transcript sequence
   * @param saveCompressGZIP use GZIP compression for output file
   * @param saveCompressBZIP2 use BZIP2 compression for output file
   * @param status the status of the ONT run
   * @param basecaller the basecaller of the ONT run
   * @param version the version of the ONT run
   * @param type the type of the ONT run
   * @param chemistryVersion the chemistryVersion of the ONT run
   */
  DirectoryProcessor(final File repertoryFastqOutput,
      final boolean saveComplementSequence, final boolean saveTemplateSequence,
      final boolean saveConsensusSequence, final boolean saveTranscriptSequence,
      final boolean saveCompressGZIP, final boolean saveCompressBZIP2,
      Fast5.Status status, Fast5.Basecaller basecaller, Fast5.Version version,
      Fast5.Type type, Fast5.ChemistryVersion chemistryVersion,
      final boolean metrichor) {

    this.repertoryFastqOutput = repertoryFastqOutput;

    this.saveComplementSequence = saveComplementSequence;
    this.saveTemplateSequence = saveTemplateSequence;
    this.saveConsensusSequence = saveConsensusSequence;
    this.saveTranscriptSequence = saveTranscriptSequence;

    this.saveCompressGZIP = saveCompressGZIP;
    this.saveCompressBZIP2 = saveCompressBZIP2;

    this.status = status;
    this.basecaller = basecaller;
    this.version = version;
    this.type = type;
    this.chemistryVersion = chemistryVersion;

    this.metrichor = metrichor;
  }

  //
  // Inner class
  //

  /**
   * This class implement the compression of the fastq output.
   */
  private static class SynchronizedWriter extends OutputStreamWriter {

    /**
     * The constructor of the abstract class SynchronizedWriter.
     * @param file, the file to be compressed
     * @param compress, the type of compression
     * @throws IOException, test if the file can be compress
     */
    private SynchronizedWriter(File file, String compress) throws IOException {
      super(getOutputStream(file, compress));
    }

    /**
     * This method of the object of the class SynchronizedWriter compress a
     * file.
     * @param file, the file to be compressed
     * @param compress, the type of compression
     * @return the file compressed
     * @throws IOException, test if the file can be compress
     */
    private static OutputStream getOutputStream(File file, String compress)
        throws IOException {
      try {

        // compress sequence to gzip
        if (compress.equals("gzip")) {
          return new GZIPOutputStream(new FileOutputStream(file));
        }

        // compress sequence to bzip2
        if (compress.equals("bzip2")) {
          return new BZip2CompressorOutputStream(new FileOutputStream(file));
        } else {
          return new FileOutputStream(file);
        }

      } catch (IOException e) {
        throw new IOException("Could not create CompressorOutputStream", e);
      }
    }
  }

  //
  // Getter
  //

  /**
   * Get the list of corrupted files
   * @return a list with the corrupted files
   */
  public List<File> getListCorruptFast5Files() {
    return this.listCorruptFast5Files;
  }

  //
  //
  //

  /**
   * This method of the class Fast5ToFastq read fast5 files on a list and write
   * the fastq sequence.
   * @param listFast5Files is the list of fast5 file
   * @param complementWriter is the writer of the complement sequence
   * @param templateWriter is the writer of the template sequence
   * @param consensusWriter is the writer of the consensus sequence
   * @param transcriptWriter is the writer of the transcript sequence
   * @param status is the status of the fast5 file
   * @throws IOException, test the read of the file
   */
  private void readFast5WriteFastq(List<File> listFast5Files,
      Writer complementWriter, Writer templateWriter, Writer consensusWriter,
      Writer transcriptWriter, String status, LocalReporter localReporter)
      throws IOException {

    // read fast5 files
    for (File fast5File : listFast5Files) {

      // process the translation of a fast5 file to the fastq
      readFast5WriteFastq(fast5File, complementWriter, templateWriter,
          consensusWriter, transcriptWriter, status, localReporter);
    }
  }

  /**
   * This method of the class Fast5ToFastq process the type of sequence and
   * launch the read of fast5 and write of fastq sequence.
   * @param listFast5Files is a list of fast5 file
   * @param status is the status of fast5 file
   * @throws IOException, test the read of the file
   */
  public void processDirectory(List<File> listFast5Files, String status,
      LocalReporter localReporter) throws IOException {

    // test if the list of fast5 files is empty
    if (listFast5Files.isEmpty()) {
      return;
    }
    // Create writters

    Writer complementWriter = null;
    Writer templateWriter = null;
    Writer consensusWriter = null;
    Writer transcriptWriter = null;

    // test if the complement sequence is to process
    if (this.saveComplementSequence) {

      // create complement Writer
      complementWriter =
          createWriterFastq(listFast5Files.get(0), "complement", status);
    }

    // test if the template sequence is to process
    if (this.saveTemplateSequence) {

      // create template Writer
      templateWriter =
          createWriterFastq(listFast5Files.get(0), "template", status);
    }

    // test if the consensus sequence is to process
    if (this.saveConsensusSequence) {

      // create consensus Writer
      consensusWriter =
          createWriterFastq(listFast5Files.get(0), "consensus", status);
    }

    // test if the transcript sequence is to process
    if (this.saveTranscriptSequence) {

      // create transcript Writer
      transcriptWriter =
          createWriterFastq(listFast5Files.get(0), "transcript", status);
    }

    // Read all Fast5 files

    // get the time of the begin execution of the translation of a fast5
    // directory into a fastq
    long start1 = System.currentTimeMillis();

    // execution of the translation of a fast5 directory into a fastq
    readFast5WriteFastq(listFast5Files, complementWriter, templateWriter,
        consensusWriter, transcriptWriter, status, localReporter);

    // get the time of the end execution of the translation of a fast5 directory
    // into a fastq
    long end1 = System.currentTimeMillis();
    System.out.println("Time exe 1 thread:"
        + (end1 - start1) / 1000 + "s for a " + listFast5Files.size()
        + " number of fast5");

    // Close writters
    if (this.saveComplementSequence) {
      assert complementWriter != null;
      complementWriter.close();
    }
    if (this.saveTemplateSequence) {
      assert templateWriter != null;
      templateWriter.close();
    }
    if (this.saveConsensusSequence) {
      assert consensusWriter != null;
      consensusWriter.close();
    }

    if (this.saveTranscriptSequence) {
      assert transcriptWriter != null;
      transcriptWriter.close();
    }
  }

  //
  //
  // Create Writer Fastq Name
  //
  //

  /**
   * This method of the class Fast5ToFastq create the good name of the fastq
   * output file.
   * @param fast5File is the name of the first file of the list "listFast5Files"
   * @param typeSequence is the type of sequence (ex:complement)
   * @param status is the status of the fast5 file (ex:fail)
   * @return a writter with the correct output name for write a fastq sequence
   * @throws IOException, test if the compression or the writing is ok
   */
  private Writer createWriterFastq(File fast5File, String typeSequence,
      String status) throws IOException {

    String preNameFile;

    String[] part = fast5File.getName().split("_");

    // test if the filename contains the "mux_scan" field or the
    // "sequencing_run"
    if (part[4].equals("mux")
        && part[5].equals("scan") && !part[4].equals("sequencing")
        && !part[5].equals("run")) {

      // double substring on the "mux_scan" field and the 5 number id field.
      preNameFile = fast5File.getName().substring(0,
          fast5File.getName().indexOf("_mux_scan") + 1)
          + fast5File.getName().substring(fast5File.getName()
              .substring(0, fast5File.getName().indexOf("_scan_") + 6).length(),
              fast5File.getName().indexOf("_ch") - 5);
    } else {

      // double substring on the "sequencing_run" field and the 5 number id
      // field.
      preNameFile = fast5File.getName().substring(0,
          fast5File.getName().indexOf("_sequencing_run") + 1)
          + fast5File.getName().substring(fast5File.getName()
              .substring(0, fast5File.getName().indexOf("_run") + 5).length(),
              fast5File.getName().indexOf("_ch") - 5);
    }

    // create writer in bzip2 compression
    if (saveCompressBZIP2) {
      return new SynchronizedWriter(
          new File(this.repertoryFastqOutput
              + "/" + preNameFile + status + "_" + typeSequence + ".fastq.bz2"),
          "bzip2");
    }

    // create writer in gzip compression
    if (saveCompressGZIP) {
      return new SynchronizedWriter(
          new File(this.repertoryFastqOutput
              + "/" + preNameFile + status + "_" + typeSequence + ".fastq.gz"),
          "gzip");
    } else {
      return new SynchronizedWriter(
          new File(this.repertoryFastqOutput
              + "/" + preNameFile + status + "_" + typeSequence + ".fastq"),
          "");
    }
  }

  /**
   * Process a sequence.
   * @param sequence sequence to process
   * @param writer Writer to use to write the sequence
   * @param localReporter local reporter
   * @param counterName the counter name
   * @throws IOException if an error occurs while writing the sequence
   */
  private static void processSequence(final String sequence,
      final Writer writer, final LocalReporter localReporter,
      final String counterName) throws IOException {

    if (sequence == null) {
      return;
    }

    // get the part of the read sequence
    final int indexCR1 = sequence.indexOf('\n');
    final int indexCR2 = sequence.indexOf('\n', indexCR1 + 1);

    // test if the sequence is null
    if (sequence.substring(indexCR1 + 1, indexCR2).equals("")) {
      localReporter.incrCounter("numberSequenceWrite", counterName + "Null", 1);
    } else {
      writer.write(sequence);
      localReporter.incrCounter("numberSequenceWrite", counterName + "Write",
          1);
    }
  }

  /**
   * Fill the counters for a FAST5 file
   * @param f5 FAST5 file
   * @param status status value
   * @param localReporter local reporter
   */
  private static void fillCounters(final Fast5 f5, final String status,
      final LocalReporter localReporter) {

    //
    // Get the Workflows Informations in the SynchronizedCountWriteFastq
    // Object
    //

    // Barcode Workflow
    //
    String counterNameBarcode = status + "_barcodeWorkflow";
    String barcodingStatus = f5.getBarcodindFinalStatus();

    // test if the barcode counter contains this status
    if (localReporter.getCounterNames(counterNameBarcode)
        .contains(barcodingStatus)) {

      // incremente this status counter
      localReporter.incrCounter(counterNameBarcode, barcodingStatus, 1);
    } else {

      // add the status
      localReporter.setCounter(counterNameBarcode, barcodingStatus, 1);
    }

    // Basecall_1D Workflow
    //
    String counterNameBasecall1D = status + "_basecall1DWorkflow";
    String basecall1dStatus = f5.getBaseCall1DFinalStatus();

    // test if the barcode counter contains this status
    if (localReporter.getCounterNames(counterNameBasecall1D)
        .contains(basecall1dStatus)) {
      // incremente this status counter
      localReporter.incrCounter(counterNameBasecall1D, basecall1dStatus, 1);
    } else {

      // add the status
      localReporter.setCounter(counterNameBasecall1D, basecall1dStatus, 1);
    }

    // Basecall_2D Workflow
    //
    String counterNameBasecall2D = status + "_basecall2DWorkflow";
    String basecall2dStatus = f5.getBaseCall2DFinalStatus();

    // test if the barcode counter contains this status
    if (localReporter.getCounterNames(counterNameBasecall2D)
        .contains(basecall2dStatus)) {

      // incremente this status counter
      localReporter.incrCounter(counterNameBasecall2D, basecall2dStatus, 1);
    } else {

      // add the status
      localReporter.setCounter(counterNameBasecall2D, basecall2dStatus, 1);
    }

    // Calibration Strand Workflow
    //
    String counterNameCalibrationStrand = status + "_calibrationStrandWorkflow";
    String calibrationStrandStatus = f5.getCalibrationStrandFinalStatus();

    // test if the barcode counter contains this status
    if (localReporter.getCounterNames(counterNameCalibrationStrand)
        .contains(calibrationStrandStatus)) {

      // incremente this status counter
      localReporter.incrCounter(counterNameCalibrationStrand,
          calibrationStrandStatus, 1);

      // test if the status is "Calibration strand detected"
      if (calibrationStrandStatus.contains("Calibration strand detected")) {

        // add the status
        localReporter.incrCounter("numberFiles",
            "numberCalibrateStrandFast5Files", 1);
      }
    } else {
      localReporter.setCounter(counterNameCalibrationStrand,
          calibrationStrandStatus, 1);
    }

    // Event Detection Workflow
    //
    String counterNameEventDetection = status + "_eventDetectionWorkflow";
    String eventDetectionFinalStatus = f5.getEventDetectionFinalStatus();

    // test if the barcode counter contains this status
    if (localReporter.getCounterNames(counterNameEventDetection)
        .contains(eventDetectionFinalStatus)) {

      // incremente this status counter
      localReporter.incrCounter(counterNameEventDetection,
          eventDetectionFinalStatus, 1);
    } else {

      // add the status
      localReporter.setCounter(counterNameEventDetection,
          eventDetectionFinalStatus, 1);
    }

    // Hairpin Split Workflow
    //
    String counterNameHairpinSplit = status + "_hairpinSplitWorkflow";
    String hairpinSplitFinalStatus = f5.getHairpinSplitFinalStatus();

    // test if the barcode counter contains this status
    if (localReporter.getCounterNames(counterNameHairpinSplit)
        .contains(hairpinSplitFinalStatus)) {

      // incremente this status counter
      localReporter.incrCounter(counterNameHairpinSplit,
          hairpinSplitFinalStatus, 1);
    } else {

      // add the status
      localReporter.setCounter(counterNameHairpinSplit, hairpinSplitFinalStatus,
          1);
    }
  }

  /**
   * This method of the class Fast5ToFastq read the fast5 file and write in the
   * fastq file.
   * @param fast5File, the fast5 file to be read
   * @param complementWriter, a fastq output file
   * @param templateWriter, a fastq output file
   * @param consensusWriter, a fastq output file
   * @param transcriptWriter, a fastq output file
   * @param status, the name of the root classification of a minion run
   * @param localReporter, the object who stores log information
   * @throws IOException, test the read of the file
   */
  private void readFast5WriteFastq(File fast5File, Writer complementWriter,
      Writer templateWriter, Writer consensusWriter, Writer transcriptWriter,
      String status, LocalReporter localReporter) throws IOException {

    // test if the fast5 is corrupt or readable
    try (Fast5 f5 = new Fast5(fast5File, this.status, this.basecaller,
        this.version, this.type, this.chemistryVersion)) {

      // test if the complementWriter is not null and if the complement sequence
      // is not null
      if (complementWriter != null) {
        processSequence(f5.getComplementFastq(), complementWriter,
            localReporter, status + "_numberSequenceComplement");
      }

      // test if the templateWriter is not null and if the template sequence is
      // not null
      if (templateWriter != null) {
        processSequence(f5.getTemplateFastq(), templateWriter, localReporter,
            status + "_numberSequenceTemplate");
      }

      // test if the consensusWriter is not null and if the consensus sequence
      // is not null
      if (consensusWriter != null) {
        processSequence(f5.getConsensusFastq(), consensusWriter, localReporter,
            status + "_numberSequenceConsensus");
      }

      // test if the transcriptWriter is not null and if the transcript sequence
      // is not null
      if (transcriptWriter != null) {
        processSequence(f5.getTranscriptFastq(), transcriptWriter,
            localReporter, status + "_numberSequenceTranscript");
      }

      // test if the basecaller is Metrichor
      if (this.metrichor) {

        // Fill the counters
        fillCounters(f5, status, localReporter);
      }

    } catch (HDF5Exception e) {

      // incremente counter for corrupt files
      localReporter.incrCounter("numberFiles", "numberCorruptFast5Files", 1);
      this.listCorruptFast5Files.add(fast5File);
    }
  }

}
