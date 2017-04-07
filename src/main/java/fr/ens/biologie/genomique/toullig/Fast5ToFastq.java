package fr.ens.biologie.genomique.toullig;

import java.io.File;
import java.io.FileOutputStream;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.zip.GZIPOutputStream;

import fr.ens.biologie.genomique.eoulsan.util.LocalReporter;
import org.apache.commons.compress.compressors.bzip2.BZip2CompressorOutputStream;

import ncsa.hdf.hdf5lib.exceptions.HDF5Exception;

/**
 * This class read a minION run of Fast5 basecalled to extract fastq sequence.
 * @author Aurelien Birer
 */
public class Fast5ToFastq {

  private File fast5RunDirectory;
  private File repertoryFastqOutput;

  private final List<File> listCorruptFast5Files = new ArrayList<>();

  private boolean processMergeStatus;
  private boolean processFail;
  private boolean processPass;
  private boolean processUnclassified;

  private boolean saveComplementSequence;
  private boolean saveTemplateSequence;
  private boolean saveConsensusSequence;
  private boolean saveTranscriptSequence;

  private boolean saveCompressGZIP;
  private boolean saveCompressBZIP2;

  private final LocalReporter localReporter = new LocalReporter();

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
  //
  // The variable rootRepertoryFast5Run contains the folder downloads with the
  // basecalled reads and the folder uploaded with the prebasecalling reads.
  //
  //
  /**
   * The constructor of the class Fast5ToFastq.
   * @param fast5RunDirectory is the repertory of the ONT run
   * @param repertoryFastqOutput is the repertory to store fastq sequence
   * @throws IOException, to test if the fast5RunDirectory and
   *           repertoryFastqOutput File exist
   */
  public Fast5ToFastq(File fast5RunDirectory, File repertoryFastqOutput)
      throws IOException {

    // test if the run directory of fast5 exist
    if (!fast5RunDirectory.exists()) {
      throw new IOException(
          "The repertory " + fast5RunDirectory + " dont exist!");
    } else {
      this.fast5RunDirectory = fast5RunDirectory;
    }

    // test if the fastq output directory exist
    if (!repertoryFastqOutput.exists()) {
      throw new IOException(
          "The repertory " + repertoryFastqOutput + " dont exist!");
    } else {
      this.repertoryFastqOutput = repertoryFastqOutput;
    }
  }

  /**
   * Values of the object SequenceType that design the sequence type.
   * @author birer
   */
  public enum SequenceType {

    COMPLEMENT("complement"), TEMPLATE("template"), CONSENSUS("consensus"),
    TRANSCRIPT("transcript");

    private final String name;

    /**
     * This method of the object of the class Fast5ToFastq get the name of the
     * sequence type.
     * @return a string of the name of the sequence type
     */
    public String getName() {
      return this.name;
    }

    /**
     * The constructor of the object of the class Fast5ToFastq in the class
     * Fast5ToFastq.
     * @param name is the sequence type
     */
    SequenceType(String name) {
      this.name = name;
    }
  }

  /**
   * This method of the class Fast5ToFastq return the name of the directory of
   * the minION run.
   * @return a string of the name of root directory
   */
  public String getNameDirectoryRunFast5() {
    return this.fast5RunDirectory.getName();
  }

  //
  //
  // Record of Files
  //
  //
  /**
   * This method of the class Fast5ToFastq list the sub-directory of a
   * directory.
   * @param dir is a directory
   * @return a list of directory (must be contains fast5 files)
   */
  private List<File> listSubDir(File dir) {

    // create new List for the results
    List<File> result = new ArrayList<>();

    // get the directory in the result list
    for (File file : dir.listFiles()) {

      // test if the file is a directory
      if (file.isDirectory()) {
        result.add(file);
      }
    }

    return result;
  }

  /**
   * This method of the class Fast5ToFastq list the fast5 files of list
   * @param subdirname is a directory
   * @return a list of fast5 file
   */
  private List<File> listFast5(String subdirname) {
    return listFast5(new File(this.fast5RunDirectory, subdirname));
  }

  /**
   * This method of the class Fast5ToFastq list the fast5 files of a directory.
   * @param fast5Dir is a directory (must be contains fast5 files)
   * @return a list of fast5 file
   */
  private static List<File> listFast5(File fast5Dir) {

    return Arrays.asList(fast5Dir.listFiles(new FilenameFilter() {

      // test if the file is a fast5 file
      @Override
      public boolean accept(File dir, String name) {
        return name.endsWith(".fast5");
      }
    }));
  }

  /**
   * This method of the class Fast5ToFastq list all the fast5 files of a run.
   * @return a list of all fast5 file
   */
  private List<File> listAllFast5() {

    // create list of results
    List<File> listFast5Files = new ArrayList<>();

    // get the directory of list of files
    for (File directory : listSubDir(
        new File(this.fast5RunDirectory, "downloads"))) {

      // test if the directory is a directory
      if (directory.isDirectory()) {

        // list all sub directories
        List<File> subDirectories = listSubDir(directory);

        // list all fast5 file
        listFast5Files.addAll(listFast5(directory));

        // get the sub-directoryof a sub-directory
        for (File subDirectory : subDirectories) {

          // test if the sub-directory is a directory
          if (subDirectory.isDirectory()) {

            // list all fast5 file
            listFast5Files.addAll(listFast5(subDirectory));
          }
        }
      }
    }
    return listFast5Files;
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
    // the substring is done on the string "_ch" who correspond to the channel
    // number
    String preNameFile = fast5File.getName().substring(0,
        fast5File.getName().indexOf("_ch") + 1);

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

  //
  //
  // Getter Number Files
  //
  //

  /**
   * This method of the class Fast5ToFastq get the number of fast5 files.
   * @param localReporter, a localReporter object
   * @return a long of the number of fast5 files
   */
  public long getNumberFast5Files(LocalReporter localReporter) {

    // test if the number of fast5 files is not negatif
    if (localReporter.getCounterValue("numberFiles", "numberFast5Files") <= 0) {
      return 0;
    }
    return localReporter.getCounterValue("numberFiles", "numberFast5Files");
  }

  /**
   * This method of the class Fast5ToFastq get the number of pass files.
   * @param localReporter, a localReporter object
   * @return a long of the number of pass files
   */
  public long getNumberPassFast5Files(LocalReporter localReporter) {

    // test if the number of pass fast5 files is not negatif
    if (localReporter.getCounterValue("numberFiles",
        "numberPassFast5Files") <= 0) {
      return 0;
    }
    return localReporter.getCounterValue("numberFiles", "numberPassFast5Files");
  }

  /**
   * This method of the class Fast5ToFastq get the number of corrupt files.
   * @param localReporter, a localReporter object
   * @return a long of the number of corrupt files
   */
  public long getNumberCorruptFast5Files(LocalReporter localReporter) {

    // test if the number of corrupt fast5 files is not negatif
    if (localReporter.getCounterValue("numberFiles",
        "numberCorruptFast5Files") <= 0) {
      return 0;
    }
    return localReporter.getCounterValue("numberFiles",
        "numberCorruptFast5Files");
  }

  /**
   * This method of the class Fast5ToFastq get the number of fail files.
   * @param localReporter, a localReporter object
   * @return a long of the number of fail files
   */
  public long getNumberFailFast5Files(LocalReporter localReporter) {

    // test if the number of fail fast5 files is not negatif
    if (localReporter.getCounterValue("numberFiles",
        "numberFailFast5Files") <= 0) {
      return 0;
    }
    return localReporter.getCounterValue("numberFiles", "numberFailFast5Files");
  }

  /**
   * This method of the class Fast5ToFastq get the number of bad barcoded files.
   * @param localReporter, a localReporter object
   * @return a long of the number of bad barcoded files
   */
  public long getNumberUnclassifiedFast5Files(LocalReporter localReporter) {

    // test if the number of unclassified fast5 files is not negatif
    if (localReporter.getCounterValue("numberFiles",
        "numberUnclassifiedFast5Files") <= 0) {
      return 0;
    }
    return localReporter.getCounterValue("numberFiles",
        "numberUnclassifiedFast5Files");
  }

  /**
   * This method of the class Fast5ToFastq get the number of calibrate strand
   * files.
   * @param localReporter, a localReporter object
   * @return a long of the number of calibrate strand files
   */
  public long getNumberCalibrateStrandFast5Files(LocalReporter localReporter) {

    // test if the number of calibrate strand fast5 files is not negatif
    if (localReporter.getCounterValue("numberFiles",
        "numberCalibrateStrandFast5Files") <= 0) {
      return 0;
    }
    return localReporter.getCounterValue("numberFiles",
        "numberCalibrateStrandFast5Files");
  }

  /**
   * This method of the class Fast5ToFastq get a localReporter object.
   * @return a LocalReporter object
   */
  public LocalReporter getLocalReporter() {
    return this.localReporter;
  }

  //
  //
  // Getter List
  //
  //

  /**
   * This method of the class Fast5ToFastq get the list of corrupt Fast5.
   * @return a list of corrupt files
   */
  public List<File> getListCorruptFast5Files() {
    return this.listCorruptFast5Files;
  }

  //
  //
  // Setter Process
  //
  //

  /**
   * This method of the class Fast5ToFastq set the differenciation of the type
   * fast5 on process.
   * @param processMergeStatus, a boolean to process separately the types of
   *          fast5 file
   */
  public void setMergeAllStatusFast5(boolean processMergeStatus) {
    this.processMergeStatus = processMergeStatus;
  }

  /**
   * This method of the class Fast5ToFastq set the type of files Fail on
   * process.
   */
  public void setProcessFail() {
    this.processFail = true;
  }

  /**
   * This method of the class Fast5ToFastq set the type of files Pass on
   * process.
   */
  public void setProcessPass() {
    this.processPass = true;
  }

  /**
   * This method of the class Fast5ToFastq set the type of files Fail Barcoded
   * on process. files to process.
   */
  public void setProcessUnclassified() {
    this.processUnclassified = true;
  }

  /**
   * This method of the class Fast5ToFastq set the type of complement sequence
   * on process.
   */
  public void setSaveComplementSequence() {
    this.saveComplementSequence = true;
  }

  /**
   * This method of the class Fast5ToFastq set the type of template sequence on
   * process.
   */
  public void setSaveTemplateSequence() {
    this.saveTemplateSequence = true;
  }

  /**
   * This method of the class Fast5ToFastq set the type of consensus sequence on
   * process.
   */
  public void setSaveConsensusSequence() {
    this.saveConsensusSequence = true;
  }

  /**
   * This method of the class Fast5ToFastq set the type of transcript sequence
   * on process.
   */
  public void setSaveTranscriptSequence() {
    this.saveTranscriptSequence = true;
  }

  //
  // Compression format setters
  //

  /**
   * This method of the class Fast5ToFastq set the type of compression of fastq
   * output to gzip.
   */
  public void setGzipCompression() {
    this.saveCompressGZIP = true;
  }

  /**
   * This method of the class Fast5ToFastq set the type of compression of fastq
   * output to bzip2.
   */
  public void setBZip2Compression() {
    this.saveCompressBZIP2 = true;
  }

  //
  //
  // Important methods
  //
  //

  /**
   * This method of the class Fast5ToFastq launch processDirectory for a
   * directory of fast5.
   * @param fast5SubdirName is a directory of fast5 files
   * @param status is the status of fast5 file
   * @return an int who is the number of fast5 process
   * @throws IOException, test the read of the file
   */
  private int processDirectory(String fast5SubdirName, String status,
      LocalReporter localReporter) throws IOException {

    // list of fast5 files
    List<File> list = listFast5(fast5SubdirName);

    // process this list of fast5 files
    processDirectory(list, status, localReporter);
    return list.size();
  }

  /**
   * This method of the class Fast5ToFastq process the type of sequence and
   * launch the read of fast5 and write of fastq sequence.
   * @param listFast5Files is a list of fast5 file
   * @param status is the status of fast5 file
   * @throws IOException, test the read of the file
   */
  private void processDirectory(List<File> listFast5Files, String status,
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

  /**
   * This method of the class Fast5ToFastq execute processDirectory with a list
   * of barcode.
   * @param listBarcodeDir is the list of barcode of the run
   * @throws IOException, test the read of the file
   */
  private void processDirectories(List<File> listBarcodeDir,
      LocalReporter localReporter) throws IOException {

    // process each barcode directories
    for (File barcodeDirectory : listBarcodeDir) {

      // create a list of fast5 for a barcode
      List<File> listBarcodeFast5Files = listFast5(barcodeDirectory);

      // process this list
      processDirectory(listBarcodeFast5Files, barcodeDirectory.getName(),
          localReporter);

      // incremente numberPassFast5Files counter
      localReporter.incrCounter("numberFiles", "numberPassFast5Files",
          listBarcodeFast5Files.size());

      // incremente numberFast5Files counter
      localReporter.incrCounter("numberFiles", "numberFast5Files",
          listBarcodeFast5Files.size());
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
    try (Fast5 f5 = new Fast5(fast5File)) {

      // test if the complementWriter is not null and if the complement sequence
      // is not null
      if (complementWriter != null && f5.getComplementFastq() != null) {

        // get the complement sequence
        String sequence = f5.getComplementFastq();

        // get the part of the read sequence
        String[] part = sequence.split("\n");

        // test if the sequence is null
        if (part[1].equals("")) {
          String counterName = status + "_numberSequenceComplementNull";
          localReporter.incrCounter("numberSequenceWrite", counterName, 1);
        } else {
          complementWriter.write(sequence);
          String counterName = status + "_numberSequenceComplementWrite";
          localReporter.incrCounter("numberSequenceWrite", counterName, 1);
        }

      }

      // test if the templateWriter is not null and if the template sequence is
      // not null
      if (templateWriter != null && f5.getTemplateFastq() != null) {

        // get the complement sequence
        String sequence = f5.getTemplateFastq();

        // get the part of the read sequence
        String[] part = sequence.split("\n");

        // test if the sequence is null
        if (part[1].equals("")) {
          String counterName = status + "_numberSequenceTemplateNull";
          localReporter.incrCounter("numberSequenceWrite", counterName, 1);
        } else {
          templateWriter.write(sequence);
          String counterName = status + "_numberSequenceTemplateWrite";
          localReporter.incrCounter("numberSequenceWrite", counterName, 1);
        }

      }

      // test if the consensusWriter is not null and if the consensus sequence
      // is not null
      if (consensusWriter != null && f5.getConsensusFastq() != null) {

        // get the complement sequence
        String sequence = f5.getConsensusFastq();

        // get the part of the read sequence
        String[] part = sequence.split("\n");

        // test if the sequence is null
        if (part[1].equals("")) {
          String counterName = status + "_numberSequenceConsensusNull";
          localReporter.incrCounter("numberSequenceWrite", counterName, 1);
        } else {
          consensusWriter.write(sequence);
          String counterName = status + "_numberSequenceConsensusWrite";
          localReporter.incrCounter("numberSequenceWrite", counterName, 1);
        }

      }

      // test if the transcriptWriter is not null and if the transcript sequence
      // is not null
      if (transcriptWriter != null && f5.getTranscriptFastq() != null) {

        // get the complement sequence
        String sequence = f5.getTranscriptFastq();

        // get the part of the read sequence
        String[] part = sequence.split("\n");

        // test if the sequence is null
        if (part[1].equals("")) {
          String counterName = status + "_numberSequenceTranscriptNull";
          localReporter.incrCounter("numberSequenceWrite", counterName, 1);
        } else {
          transcriptWriter.write(sequence);
          String counterName = status + "_numberSequenceTranscriptWrite";
          localReporter.incrCounter("numberSequenceWrite", counterName, 1);
        }

      }

      //
      // Get the Workflows Informations in the SynchronizedCountWriteFastq
      // Object
      //

      // Barcode Workflow
      //
      String counterNameBarcode = status + "_barcodeWorkflow";

      // test if the barcode counter contains this status
      if (localReporter.getCounterNames(counterNameBarcode)
          .contains(f5.getBarcodindFinalStatus())) {

        // incremente this status counter
        localReporter.incrCounter(counterNameBarcode,
            f5.getBarcodindFinalStatus(), 1);
      } else {

        // add the status
        localReporter.setCounter(counterNameBarcode,
            f5.getBarcodindFinalStatus(), 1);
      }

      // Basecall_1D Workflow
      //
      String counterNameBasecall1D = status + "_basecall1DWorkflow";

      // test if the barcode counter contains this status
      if (localReporter.getCounterNames(counterNameBasecall1D)
          .contains(f5.getBaseCall1DFinalStatus())) {
        // incremente this status counter
        localReporter.incrCounter(counterNameBasecall1D,
            f5.getBaseCall1DFinalStatus(), 1);
      } else {

        // add the status
        localReporter.setCounter(counterNameBasecall1D,
            f5.getBaseCall1DFinalStatus(), 1);
      }

      // Basecall_2D Workflow
      //
      String counterNameBasecall2D = status + "_basecall2DWorkflow";

      // test if the barcode counter contains this status
      if (localReporter.getCounterNames(counterNameBasecall2D)
          .contains(f5.getBaseCall2DFinalStatus())) {

        // incremente this status counter
        localReporter.incrCounter(counterNameBasecall2D,
            f5.getBaseCall2DFinalStatus(), 1);
      } else {

        // add the status
        localReporter.setCounter(counterNameBasecall2D,
            f5.getBaseCall2DFinalStatus(), 1);
      }

      // Calibration Strand Workflow
      //
      String counterNameCalibrationStrand =
          status + "_calibrationStrandWorkflow";

      // test if the barcode counter contains this status
      if (localReporter.getCounterNames(counterNameCalibrationStrand)
          .contains(f5.getCalibrationStrandFinalStatus())) {

        // incremente this status counter
        localReporter.incrCounter(counterNameCalibrationStrand,
            f5.getCalibrationStrandFinalStatus(), 1);

        // test if the status is "Calibration strand detected"
        if (f5.getCalibrationStrandFinalStatus()
            .contains("Calibration strand detected")) {

          // add the status
          localReporter.incrCounter("numberFiles",
              "numberCalibrateStrandFast5Files", 1);
        }
      } else {
        localReporter.setCounter(counterNameCalibrationStrand,
            f5.getCalibrationStrandFinalStatus(), 1);
      }

      // Event Detection Workflow
      //
      String counterNameEventDetection = status + "_eventDetectionWorkflow";

      // test if the barcode counter contains this status
      if (localReporter.getCounterNames(counterNameEventDetection)
          .contains(f5.getEventDetectionFinalStatus())) {

        // incremente this status counter
        localReporter.incrCounter(counterNameEventDetection,
            f5.getEventDetectionFinalStatus(), 1);
      } else {

        // add the status
        localReporter.setCounter(counterNameEventDetection,
            f5.getEventDetectionFinalStatus(), 1);
      }

      // Hairpin Split Workflow
      //
      String counterNameHairpinSplit = status + "_hairpinSplitWorkflow";

      // test if the barcode counter contains this status
      if (localReporter.getCounterNames(counterNameHairpinSplit)
          .contains(f5.getHairpinSplitFinalStatus())) {

        // incremente this status counter
        localReporter.incrCounter(counterNameHairpinSplit,
            f5.getHairpinSplitFinalStatus(), 1);
      } else {

        // add the status
        localReporter.setCounter(counterNameHairpinSplit,
            f5.getHairpinSplitFinalStatus(), 1);
      }

      //
      //
      //

    } catch (HDF5Exception e) {

      // incremente counter for corrupt files
      localReporter.incrCounter("numberFiles", "numberCorruptFast5Files", 1);
      this.listCorruptFast5Files.add(fast5File);
    }
  }

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

  //
  //
  // Main function
  //
  //

  /**
   * This method of the class Fast5ToFastq execute the process to retrieve the
   * fastq sequence on fast5 file.
   * @throws IOException, test the read of the file
   */
  public void execute() throws IOException {

    // test if the merge of fastq is enable
    if (this.processMergeStatus) {

      // process all fast5
      processDirectory(listAllFast5(), "merge_status", this.localReporter);
      return;
    }

    // test if the fail fast5 is to process
    if (this.processFail) {

      // get the list of fail fast5 files
      int numberFailFast5Files =
          processDirectory("downloads/fail", "fail", this.localReporter);

      // incremente fail fast5 files counter
      this.localReporter.incrCounter("numberFiles", "numberFailFast5Files",
          numberFailFast5Files);

      // incremente fast5 files counter
      this.localReporter.incrCounter("numberFiles", "numberFast5Files",
          numberFailFast5Files);
    }

    // test if the pass fast5 is to process
    if (this.processPass) {

      // get the list of pass fast5 files
      int numberPassFast5Files =
          processDirectory("downloads/pass", "pass", this.localReporter);

      // test if the number of pass fast5 files is not null
      if (numberPassFast5Files != 0) {

        // incremente pass fast5 files counter
        this.localReporter.incrCounter("numberFiles", "numberPassFast5Files",
            numberPassFast5Files);

        // incremente fast5 files counter
        this.localReporter.incrCounter("numberFiles", "numberFast5Files",
            numberPassFast5Files);
      } else {

        // get the list of barcode pass fast5 files
        List<File> listBarcodeFast5Dir =
            listSubDir(new File(this.fast5RunDirectory, "downloads/pass"));

        // process the barcode directories
        processDirectories(listBarcodeFast5Dir, this.localReporter);
      }

    }

    // test if the unclassified fast5 is to process
    if (this.processUnclassified) {

      // get the list of unclassified fast5 files
      int numberUnclassifiedFast5Files = processDirectory(
          "downloads/fail/unclassified", "unclassified", this.localReporter);

      // incremente unclassified fast5 files counter
      this.localReporter.incrCounter("numberFiles",
          "numberUnclassifiedFast5Files", numberUnclassifiedFast5Files);

      // incremente fast5 files counter
      this.localReporter.incrCounter("numberFiles", "numberFast5Files",
          numberUnclassifiedFast5Files);
    }
  }
}