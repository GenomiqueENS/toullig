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
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.zip.GZIPOutputStream;

import org.apache.commons.compress.compressors.bzip2.BZip2CompressorOutputStream;

import ncsa.hdf.hdf5lib.exceptions.HDF5Exception;

/**
 * This class read a minION run of Fast5 basecalled to extract fastq sequence.
 * @author Aurelien Birer
 */
public class Fast5ToFastq {

  private File fast5RunDirectory;
  private File repertoryFastqOutput;

  private List<File> listCorruptFast5Files = new ArrayList<>();

  private boolean processMergeStatus;
  private boolean processFail;
  private boolean processPass;
  private boolean processUnclassified;
  private boolean processPassBarcode;

  private boolean saveComplementSequence;
  private boolean saveTemplateSequence;
  private boolean saveConsensusSequence;
  private boolean saveTranscriptSequence;

  private boolean saveCompressGZIP;
  private boolean saveCompressBZIP2;

  private LocalReporter localReporter = new LocalReporter();

  /**
   * This class implement the compression of the fastq output.
   */
  private static class SynchronizedWriter extends OutputStreamWriter {

    /**
     *  The constructor of the abstract class SynchronizedWriter.
     * @param file, the file to be compressed
     * @param compress, the type of compression
     * @throws IOException, test if the file can be compress
     */
    private SynchronizedWriter(File file, String compress) throws IOException {
      super(getOutputStream(file, compress));
    }

    /**
     * This method of the object of the class SynchronizedWriter compress a file.
     * @param file, the file to be compressed
     * @param compress, the type of compression
     * @return the file compressed
     * @throws IOException, test if the file can be compress
     */
    private static OutputStream getOutputStream(File file, String compress)
        throws IOException {
      try {
        if (compress.equals("gzip")) {
          return new GZIPOutputStream(new FileOutputStream(file));
        }
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
   * @throws IOException, to test if the fast5RunDirectory and repertoryFastqOutput File exist
   */
  public Fast5ToFastq(File fast5RunDirectory, File repertoryFastqOutput)
      throws IOException {

    if (!fast5RunDirectory.exists()) {
      throw new IOException(
          "The repertory " + fast5RunDirectory + " dont exist!");
    } else {
      this.fast5RunDirectory = fast5RunDirectory;
    }
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

    COMPLEMENT("complement"), TEMPLATE("template"), CONSENSUS("consensus"), TRANSCRIPT("transcript");

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
    List<File> result = new ArrayList<>();
    for (File file : dir.listFiles()) {
      if (file.isDirectory()) {
        result.add(file);
      }
    }
    return result;
  }

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
    List<File> listFast5Files = new ArrayList<File>();
    for (File directory : listSubDir(
        new File(this.fast5RunDirectory, "downloads"))) {
      if (directory.isDirectory()) {
        List<File> subDirectories = listSubDir(directory);
        listFast5Files.addAll(listFast5(directory));
        for (File subDirectory : subDirectories) {
          if (subDirectory.isDirectory()) {
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

    if (saveCompressBZIP2) {
      return new SynchronizedWriter(
          new File(this.repertoryFastqOutput
              + "/" + preNameFile + status + "_" + typeSequence + ".bzip2"),
          "bzip2");
    }
    if (saveCompressGZIP) {
      return new SynchronizedWriter(
          new File(this.repertoryFastqOutput
              + "/" + preNameFile + status + "_" + typeSequence + ".gzip"),
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
    if(localReporter.getCounterValue("numberFiles",
            "numberFast5Files")<=0){
      return       0;
    }
    return       localReporter.getCounterValue("numberFiles",
            "numberFast5Files");
  }

  /**
   * This method of the class Fast5ToFastq get the number of pass files.
   * @param localReporter, a localReporter object
   * @return a long of the number of pass files
   */
  public long getNumberPassFast5Files(LocalReporter localReporter) {
    if(localReporter.getCounterValue("numberFiles",
            "numberPassFast5Files")<=0){
      return       0;
    }
    return localReporter.getCounterValue("numberFiles",
            "numberPassFast5Files");
  }

  /**
   * This method of the class Fast5ToFastq get the number of corrupt files.
   * @param localReporter, a localReporter object
   * @return a long of the number of corrupt files
   */
  public long getNumberCorruptFast5Files(LocalReporter localReporter) {
    if(localReporter.getCounterValue("numberFiles",
            "numberCorruptFast5Files")<=0){
      return       0;
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
    if(localReporter.getCounterValue("numberFiles",
            "numberFailFast5Files")<=0){
      return       0;
    }
    return localReporter.getCounterValue("numberFiles",
            "numberFailFast5Files");
  }

  /**
   * This method of the class Fast5ToFastq get the number of bad barcoded files.
   * @param localReporter, a localReporter object
   * @return a long of the number of bad barcoded files
   */
  public long getNumberUnclassifiedFast5Files(LocalReporter localReporter) {
    if(localReporter.getCounterValue("numberFiles",
            "numberUnclassifiedFast5Files")<=0){
      return       0;
    }
    return localReporter.getCounterValue("numberFiles",
            "numberUnclassifiedFast5Files");
  }

  /**
   * This method of the class Fast5ToFastq get the number of barcoded files.
   * @param localReporter, a localReporter object
   * @return a long of the number of barcoded files
   */
  public long getNumberBarcodeFast5Files(LocalReporter localReporter) {
    if(localReporter.getCounterValue("numberFiles",
            "numberBarcodeFast5Files")<=0){
      return       0;
    }
    return localReporter.getCounterValue("numberFiles",
            "numberBarcodeFast5Files");
  }

  /**
   * This method of the class Fast5ToFastq get the number of calibrate strand files.
   * @param localReporter, a localReporter object
   * @return a long of the number of calibrate strand files
   */
  public long getNumberCalibrateStrandFast5Files(LocalReporter localReporter) {
    if(localReporter.getCounterValue("numberFiles",
            "numberCalibrateStrandFast5Files")<=0){
      return       0;
    }
    return localReporter.getCounterValue("numberFiles",
            "numberCalibrateStrandFast5Files");
  }

  /**
   * This method of the class Fast5ToFastq get a localReporter object.
   * @return a LocalReporter object
   */
  public LocalReporter getLocalReporter(){
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
  public List<File> getListCorruptFast5Files(){
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
   * @param processMergeStatus, a boolean to process separately the types of fast5
   *          file
   */
  public void setMergeAllStatusFast5(boolean processMergeStatus) {
      this.processMergeStatus = processMergeStatus;
  }

  /**
   * This method of the class Fast5ToFastq set the type of files Fail on
   * process.
   * @param processFail, a boolean to process the type of files Fail
   */
  public void setProcessFail(boolean processFail) {
      this.processFail = processFail;
  }

  /**
   * This method of the class Fast5ToFastq set the type of files Pass on
   * process.
   * @param processPass, a boolean to process the type of files Pass
   */
  public void setProcessPass(boolean processPass) {
      this.processPass = processPass;
  }

  /**
   * This method of the class Fast5ToFastq set the type of files Fail Barcoded
   * on process. files to process.
   * @param processUnclassified, a boolean to process the type of files Fail
   *          Barcoded
   */
  public void setProcessUnclassified(boolean processUnclassified) {
      this.processUnclassified = processUnclassified;
  }

  /**
   * This method of the class Fast5ToFastq set the type of files Pass Barcoded
   * on process.
   * @param processPassBarcode, a boolean to process the type of files Pass
   *          Barcoded
   */
  public void setProcessPassBarcode(boolean processPassBarcode) {
      this.processPassBarcode = processPassBarcode;
  }

  /**
   * This method of the class Fast5ToFastq set the type of complement sequence
   * on process.
   * @param saveComplementSequence, a boolean to process the type of complement
   *          sequences
   */
  public void setSaveComplementSequence(boolean saveComplementSequence) {
      this.saveComplementSequence = saveComplementSequence;
  }

  /**
   * This method of the class Fast5ToFastq set the type of template sequence on
   * process.
   * @param saveTemplateSequence, a boolean to process the type of template
   *          sequences
   */
  public void setSaveTemplateSequence(boolean saveTemplateSequence) {
      this.saveTemplateSequence = saveTemplateSequence;
  }

  /**
   * This method of the class Fast5ToFastq set the type of consensus sequence on
   * process.
   * @param saveConsensusSequence, a boolean to process the type of barcode
   *          sequences
   */
  public void setSaveConsensusSequence(boolean saveConsensusSequence) {
      this.saveConsensusSequence = saveConsensusSequence;
  }

  /**
   * This method of the class Fast5ToFastq set the type of transcript sequence on
   * process.
   * @param saveTranscriptSequence, a boolean to process the type of barcode
   *          sequences
   */
  public void setSaveTranscriptSequence(boolean saveTranscriptSequence) {
    this.saveTranscriptSequence = saveTranscriptSequence;
  }


  //
  // Compression format setters
  //

  /**
   * This method of the class Fast5ToFastq set the type of compression of fastq output to gzip.
   * @param saveCompressGZIP, boolean to process the type of compression
   */
  public void setGzipCompression(boolean saveCompressGZIP) {
      this.saveCompressGZIP = saveCompressGZIP;
  }

  /**
   * This method of the class Fast5ToFastq set the type of compression of fastq output to bzip2.
   * @param saveCompressBZIP2, boolean to process the type of compression
   */
  public void setBZip2Compression(boolean saveCompressBZIP2) {
      this.saveCompressBZIP2 = saveCompressBZIP2;
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

    List<File> list = listFast5(fast5SubdirName);
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

    if (listFast5Files.isEmpty()) {
      return;
    }
    // Create writters
    Writer complementWriter = null;
    Writer templateWriter = null;
    Writer consensusWriter = null;
    Writer transcriptWriter = null;

    if (this.saveComplementSequence) {
      complementWriter =
          createWriterFastq(listFast5Files.get(0), "complement", status);
    }
    if (this.saveTemplateSequence) {
      templateWriter =
          createWriterFastq(listFast5Files.get(0), "template", status);
    }
    if (this.saveConsensusSequence) {
      consensusWriter =
          createWriterFastq(listFast5Files.get(0), "consensus", status);
    }
    if (this.saveTranscriptSequence) {
      transcriptWriter =
              createWriterFastq(listFast5Files.get(0), "transcript", status);
    }


    // Read all Fast5 files

    long start1 = System.currentTimeMillis();

    readFast5WriteFastq(listFast5Files, complementWriter, templateWriter,
        consensusWriter, transcriptWriter, status, localReporter);

    long end1 = System.currentTimeMillis();
    System.out.println("Time exe 1 thread:"
        + (end1 - start1) / 1000 + "s for a " + listFast5Files.size()
        + " number of fast5");


//    for (int i = 0; i <= listFast5Files.size(); i += 50000) {
//
//     long start2 = System.currentTimeMillis();
//     multiThreadReadFast5WriteFastq(listFast5Files.subList(0, i),
//     complementWriter, templateWriter, consensusWriter, transcriptWriter, status, localReporter);
//
//     long end2 = System.currentTimeMillis();
//
//     System.out.println("Time exe multi thread :"
//     + (end2 - start2) / 1000 + "s for a " + i + " number of fast5");
//     }
//     long start3 = System.currentTimeMillis();
//     multiThreadReadFast5WriteFastq(listFast5Files, complementWriter,
//     templateWriter, consensusWriter, transcriptWriter, status, localReporter);
//
//     long end3 = System.currentTimeMillis();
//
//     System.out
//     .println("Time exe multi thread :" + (end3 - start3) / 1000 + "s");
//
//    long start4 = System.currentTimeMillis();
//
//    readFast5WriteFastq(listFast5Files, complementWriter, templateWriter,
//            consensusWriter, transcriptWriter, status, localReporter);
//
//    long end4 = System.currentTimeMillis();
//    System.out.println("Time exe 1 thread:"
//            + (end4 - start4) / 1000 + "s for a " + listFast5Files.size()
//            + " number of fast5");

    // Close writters
    if (this.saveComplementSequence) {
      complementWriter.close();
    }
    if (this.saveTemplateSequence) {
      templateWriter.close();
    }
    if (this.saveConsensusSequence) {
      consensusWriter.close();
    }

    if (this.saveTranscriptSequence) {
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
    for (File barcodeDirectory : listBarcodeDir) {
      List<File> listBarcodeFast5Files = listFast5(barcodeDirectory);
      processDirectory(listBarcodeFast5Files, barcodeDirectory.getName(),
          localReporter);

      localReporter.incrCounter("numberFiles",
              "numberBarcodeFast5Files", listBarcodeFast5Files.size());
      localReporter.incrCounter("numberFiles",
              "numberFast5Files", listBarcodeFast5Files.size());
    }
  }

  /**
   * This method of the class Fast5ToFastq read the fast5 file and write in the fastq file.
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
      Writer templateWriter, Writer consensusWriter, Writer transcriptWriter, String status,
      LocalReporter localReporter) throws IOException {

    try (Fast5 f5 = new Fast5(fast5File)) {
      if (complementWriter != null && f5.getComplementFastq() != null) {
        complementWriter.write(f5.getComplementFastq());
        String counterName = status+"_numberSequenceComplementWrite";
        localReporter.incrCounter("numberSequenceWrite",
                counterName, 1);
      }
      if (templateWriter != null && f5.getTemplateFastq() != null) {
        templateWriter.write(f5.getTemplateFastq());
        String counterName = status + "_numberSequenceTemplateWrite";
        localReporter.incrCounter("numberSequenceWrite",
                counterName, 1);
      }
      if (consensusWriter != null && f5.getConsensusFastq() != null) {
        consensusWriter.write(f5.getConsensusFastq());
        String counterName = status+"_numberSequenceConsensusWrite";
        localReporter.incrCounter("numberSequenceWrite",
                counterName, 1);
      }
      if (transcriptWriter != null && f5.getTranscriptFastq() != null) {
        transcriptWriter.write(f5.getTranscriptFastq());
        String counterName = status+"_numberSequenceTranscriptWrite";
        localReporter.incrCounter("numberSequenceWrite",
                counterName, 1);
      }

      //
      // Get the Workflows Informations in the SynchronizedCountWriteFastq
      // Object
      //

      // Barcode Workflow
      //
      String counterNameBarcode = status+"_barcodeWorkflow";
      if (localReporter.getCounterNames(counterNameBarcode)
          .contains(f5.getBarcodindFinalStatus())) {
        localReporter.incrCounter(counterNameBarcode, f5.getBarcodindFinalStatus(),
            1);
      } else {
        localReporter.setCounter(counterNameBarcode, f5.getBarcodindFinalStatus(),
            1);
      }

      // Basecall_1D Workflow
      //
      String counterNameBasecall1D = status+"_basecall1DWorkflow";
      if (localReporter.getCounterNames(counterNameBasecall1D)
          .contains(f5.getBaseCall1DFinalStatus())) {
        localReporter.incrCounter(counterNameBasecall1D,
            f5.getBaseCall1DFinalStatus(), 1);
      } else {
        localReporter.setCounter(counterNameBasecall1D, f5.getBaseCall1DFinalStatus(),
            1);
      }

      // Basecall_2D Workflow
      //
      String counterNameBasecall2D = status+"_basecall2DWorkflow";
      if (localReporter.getCounterNames(counterNameBasecall2D)
          .contains(f5.getBaseCall2DFinalStatus())) {
        localReporter.incrCounter(counterNameBasecall2D,
            f5.getBaseCall2DFinalStatus(), 1);
      } else {
        localReporter.setCounter(counterNameBasecall2D, f5.getBaseCall2DFinalStatus(),
            1);
      }

      // Calibration Strand Workflow
      //
      String counterNameCalibrationStrand = status+"_calibrationStrandWorkflow";
      if (localReporter.getCounterNames(counterNameCalibrationStrand)
          .contains(f5.getCalibrationStrandFinalStatus())) {
        localReporter.incrCounter(counterNameCalibrationStrand,
            f5.getCalibrationStrandFinalStatus(), 1);
        if(f5.getCalibrationStrandFinalStatus().contains("Calibration strand detected")){
          localReporter.incrCounter("numberFiles",
                  "numberCalibrateStrandFast5Files", 1);
        }
      } else {
        localReporter.setCounter(counterNameCalibrationStrand,
            f5.getCalibrationStrandFinalStatus(), 1);
      }

      // Event Detection Workflow
      //
      String counterNameEventDetection = status+"_eventDetectionWorkflow";
      if (localReporter.getCounterNames(counterNameEventDetection)
          .contains(f5.getEventDetectionFinalStatus())) {
        localReporter.incrCounter(counterNameEventDetection,
            f5.getEventDetectionFinalStatus(), 1);
      } else {
        localReporter.setCounter(counterNameEventDetection,
            f5.getEventDetectionFinalStatus(), 1);
      }

      // Hairpin Split Workflow
      //
      String counterNameHairpinSplit = status+"_hairpinSplitWorkflow";
      if (localReporter.getCounterNames(counterNameHairpinSplit)
          .contains(f5.getHairpinSplitFinalStatus())) {
        localReporter.incrCounter(counterNameHairpinSplit,
            f5.getHairpinSplitFinalStatus(), 1);
      } else {
        localReporter.setCounter(counterNameHairpinSplit,
            f5.getHairpinSplitFinalStatus(), 1);
      }

      //
      //
      //

    } catch (HDF5Exception e) {
      localReporter.incrCounter("numberFiles",
              "numberCorruptFast5Files", 1);
      this.listCorruptFast5Files.add(fast5File);
    }
  }

  /**
   * This method of the class Fast5ToFastq read the fast5 file and write in the fastq file in multi-threading.
   * @param listFast5Files, the fast5 file to be read
   * @param complementWriter, a fastq output file
   * @param templateWriter, a fastq output file
   * @param consensusWriter, a fastq output file
   * @param transcriptWriter, a fastq output file
   * @param status, the name of the root classification of a minion run
   * @param localReporter, the object who stores log information
   * @throws IOException, test the read of the file
   */
  private void multiThreadReadFast5WriteFastq(List<File> listFast5Files,
      final Writer complementWriter, final Writer templateWriter,
      final Writer consensusWriter, final Writer transcriptWriter, final String status,
      final LocalReporter localReporter) throws IOException {

    ExecutorService executor = Executors
        .newFixedThreadPool(Runtime.getRuntime().availableProcessors());

    for (final File fast5File : listFast5Files) {
      executor.submit(new Runnable() {

        @Override
        public void run() {
          try {
            readFast5WriteFastq(fast5File, complementWriter, templateWriter,
                consensusWriter, transcriptWriter, status, localReporter);
          } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
          }
        }
      });
    }

    executor.shutdown();
    while (!executor.isTerminated()) {
      try {
        Thread.sleep(2000);
      } catch (InterruptedException e) {
        // TODO Auto-generated catch block
        e.printStackTrace();
      }
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
      Writer complementWriter, Writer templateWriter, Writer consensusWriter, Writer transcriptWriter,
      String status, LocalReporter localReporter) throws IOException {

    for (File fast5File : listFast5Files) {
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

    if (this.processMergeStatus) {
      processDirectory(listAllFast5(), "merge_status", this.localReporter);
      return;
    }

    if (this.processFail) {
      int numberFailFast5Files =
          processDirectory("downloads/fail", "fail", this.localReporter);
      this.localReporter.incrCounter("numberFiles",
              "numberFailFast5Files", numberFailFast5Files);
      this.localReporter.incrCounter("numberFiles",
              "numberFast5Files", numberFailFast5Files);
    }
    if (this.processPass) {
      int numberPassFast5Files =
          processDirectory("downloads/pass", "pass", this.localReporter);
      this.localReporter.incrCounter("numberFiles",
              "numberPassFast5Files", numberPassFast5Files);
      this.localReporter.incrCounter("numberFiles",
              "numberFast5Files", numberPassFast5Files);
    }
    if (this.processUnclassified) {
      int numberUnclassifiedFast5Files = processDirectory(
          "downloads/fail/unclassified", "unclassified", this.localReporter);
      this.localReporter.incrCounter("numberFiles",
              "numberUnclassifiedFast5Files", numberUnclassifiedFast5Files);
      this.localReporter.incrCounter("numberFiles",
              "numberFast5Files", numberUnclassifiedFast5Files);
    }
    if (this.processPassBarcode) {
      List<File> listBarcodeFast5Dir =
          listSubDir(new File(this.fast5RunDirectory, "downloads/pass"));
      processDirectories(listBarcodeFast5Dir, this.localReporter);
    }
  }
}
