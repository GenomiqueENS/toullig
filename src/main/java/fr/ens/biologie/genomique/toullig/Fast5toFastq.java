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
import java.util.Set;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.zip.GZIPOutputStream;

import org.apache.commons.compress.compressors.bzip2.BZip2CompressorOutputStream;

import ncsa.hdf.hdf5lib.exceptions.HDF5Exception;

/**
 * This class read a minION run of Fast5 basecalled to extract fastq sequence.
 * @author Aurelien Birer
 */
public class Fast5toFastq {

  private File fast5RunDirectory;
  private File repertoryFastqOutput;

  private List<File> listFast5Files = new ArrayList<File>();
  private List<File> listCorruptFast5Files = new ArrayList<File>();

  private int numberFast5Files;
  private int numberFailFast5Files;
  private int numberPassFast5Files;
  private int numberFailBarcodeFast5Files;
  private int numberBarcodeFast5Files;
  private int numberCorruptFast5Files;

  private boolean processMergeStatus = false;
  private boolean processFail = false;
  private boolean processPass = false;
  private boolean processFailBarcode = false;
  private boolean processPassBarcode = false;

  private boolean saveComplementSequence = false;
  private boolean saveTemplateSequence = false;
  private boolean saveBarcodeSequence = false;

  private boolean saveCompressGZIP = false;
  private boolean saveCompressBZIP2 = false;

  private List<String> listWriteSequenceLog = new ArrayList<String>();
  private List<String> listWorkflowStatusLog = new ArrayList<String>();

  private static class SynchronizedWriter extends OutputStreamWriter {

    private SynchronizedWriter(File file, String compress) throws Exception {
      super(getOutputStream(file, compress));
    }

    private static OutputStream getOutputStream(File file, String compress)
        throws Exception {
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
        throw new Exception("Could not create CompressorOutputStream", e);
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
   * The constructor of the class Fast5toFastq.
   * @param fast5RunDirectory is the repertory of the ONT run
   * @param repertoryFastqOutput is the repertory to store fastq sequence
   * @throws IOException
   */
  public Fast5toFastq(File fast5RunDirectory, File repertoryFastqOutput)
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

    COMPLEMENT("complement"), TEMPLATE("template"), BARCODE("barcode");

    private final String name;

    /**
     * This method of the object of the class Fast5toFastq get the name of the
     * sequence type.
     * @return a string of the name of the sequence type
     */
    public String getName() {
      return this.name;
    }

    /**
     * The constructor of the object of the class Fast5toFastq in the class
     * Fast5toFastq.
     * @param name is the sequence type
     */
    SequenceType(String name) {
      this.name = name;
    }
  }

  /**
   * This method of the class Fast5toFastq return the name of the directory of
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
   * This method of the class Fast5toFastq list the sub-directory of a
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
   * This method of the class Fast5toFastq list the fast5 files of a directory.
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
   * This method of the class Fast5toFastq list all the fast5 files of a run.
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
   * This method of the class Fast5toFastq create the good name of the fastq
   * output file.
   * @param fast5File is the name of the first file of the list "listFast5Files"
   * @param typeSequence is the type of sequence (ex:complement)
   * @param status is the status of the fast5 file (ex:fail)
   * @return a writter with the correct output name for write a fastq sequence
   * @throws Exception
   */
  private Writer createWriterFastq(File fast5File, String typeSequence,
      String status) throws Exception {
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
   * This method of the class Fast5toFastq get the number of fast5 files.
   * @return an int of the number of fast5 files
   */
  public int getNumberFast5Files() {
    return this.numberFast5Files;
  }

  /**
   * This method of the class Fast5toFastq get the number of pass files.
   * @return an int of the number of pass files
   */
  public int getNumberPassFast5Files() {
    return this.numberPassFast5Files;
  }

  /**
   * This method of the class Fast5toFastq get the number of corrupt files.
   * @return an int of the number of corrupt files
   */
  public int getNumberCorruptFast5Files() {
    return this.numberCorruptFast5Files;
  }

  /**
   * This method of the class Fast5toFastq get the number of fail files.
   * @return an int of the number of fail files
   */
  public int getNumberFailFast5Files() {
    return this.numberFailFast5Files;
  }

  /**
   * This method of the class Fast5toFastq get the number of bad barcoded files.
   * @return an int of the number of bad barcoded files
   */
  public int getNumberBadBarcodeFast5Files() {
    return this.numberFailBarcodeFast5Files;
  }

  /**
   * This method of the class Fast5toFastq get the number of barcoded files.
   * @return an int of the number of barcoded files
   */
  public int getNumberBarcodeFast5Files() {
    return this.numberBarcodeFast5Files;
  }

  //
  //
  // Getter Process
  //
  //

  /**
   * This method of the class Fast5toFastq get boolean value of the fail
   * process.
   * @return a boolean
   */
  public boolean getProcessFail() {
    return this.processFail;
  }

  /**
   * This method of the class Fast5toFastq get boolean value of the pass
   * process.
   * @return a boolean
   */
  public boolean getProcessPass() {
    return this.processPass;
  }

  /**
   * This method of the class Fast5toFastq get boolean value of the fail barcode
   * process.
   * @return a boolean
   */
  public boolean getProcessFailBarcode() {
    return this.processFailBarcode;
  }

  /**
   * This method of the class Fast5toFastq get boolean value of the pass barcode
   * process.
   * @return a boolean
   */
  public boolean getProcessPassBarcode() {
    return this.processPassBarcode;
  }

  //
  //
  // Setter Process
  //
  //

  /**
   * This method of the class Fast5toFastq set the differenciation of the type
   * fast5 on process.
   * @param processStatus, a boolean to process separately the types of fast5
   *          file
   */
  public void setMergeAllStatusFast5(boolean processMergeStatus) {
    if (processMergeStatus) {
      this.processMergeStatus = true;
    }
  }

  /**
   * This method of the class Fast5toFastq set the type of files Fail on
   * process.
   * @param processFail, a boolean to process the type of files Fail
   */
  public void setProcessFail(boolean processFail) {
    if (processFail) {
      this.processFail = true;
    }
  }

  /**
   * This method of the class Fast5toFastq set the type of files Pass on
   * process.
   * @param processPass, a boolean to process the type of files Pass
   */
  public void setProcessPass(boolean processPass) {
    if (processPass) {
      this.processPass = true;
    }
  }

  /**
   * This method of the class Fast5toFastq set the type of files Fail Barcoded
   * on process. files to process.
   * @param processFailBarcode, a boolean to process the type of files Fail
   *          Barcoded
   */
  public void setProcessFailBarcode(boolean processFailBarcode) {
    if (processFailBarcode) {
      this.processFailBarcode = true;
    }
  }

  /**
   * This method of the class Fast5toFastq set the type of files Pass Barcoded
   * on process.
   * @param processPassBarcode, a boolean to process the type of files Pass
   *          Barcoded
   */
  public void setProcessPassBarcode(boolean processPassBarcode) {
    if (processPassBarcode) {
      this.processPassBarcode = true;
    }
  }

  /**
   * This method of the class Fast5toFastq set the type of complement sequence
   * on process.
   * @param saveComplementSequence, a boolean to process the type of complement
   *          sequences
   */
  public void setSaveComplementSequence(boolean saveComplementSequence) {
    if (saveComplementSequence) {
      this.saveComplementSequence = true;
    }
  }

  /**
   * This method of the class Fast5toFastq set the type of template sequence on
   * process.
   * @param saveTemplateSequence, a boolean to process the type of template
   *          sequences
   */
  public void setSaveTemplateSequence(boolean saveTemplateSequence) {
    if (saveTemplateSequence) {
      this.saveTemplateSequence = true;
    }
  }

  /**
   * This method of the class Fast5toFastq set the type of barcode sequence on
   * process.
   * @param saveBarcodeSequence, a boolean to process the type of barcode
   *          sequences
   */
  public void setSaveBarcodeSequence(boolean saveBarcodeSequence) {
    if (saveBarcodeSequence) {
      this.saveBarcodeSequence = true;
    }
  }

  //
  // Compression format setters
  //

  public void setCompressGZIP(boolean saveCompressGZIP) {
    if (saveCompressGZIP) {
      this.saveCompressGZIP = true;
    }
  }

  public void setCompressBZIP2(boolean saveCompressBZIP2) {
    if (saveCompressBZIP2) {
      this.saveCompressBZIP2 = true;
    }
  }

  //
  //
  // Create Log
  //
  //

  /**
   * This method of the class Fast5toFastq get a list of log on the number of
   * files.
   * @return a list of string
   */
  public List<String> getListLog() {
    List<String> listLog = new ArrayList<String>();
    listLog.add("Number of fast5 files read " + getNumberFast5Files());
    listLog.add("Number of corrupt files " + getNumberCorruptFast5Files());
    listLog.add("Number of fail files read " + getNumberFailFast5Files());
    listLog.add("Number of pass files read " + getNumberPassFast5Files());
    listLog.add("Number of fail attribution barcode files read "
        + getNumberBadBarcodeFast5Files());
    listLog.add(
        "Number of pass barcode file " + getNumberBarcodeFast5Files() + "\n");
    for (String element : this.listWriteSequenceLog) {
      listLog.add(element);
    }
    return listLog;
  }

  /**
   * This method of the class Fast5toFastq get a list of log on the status of
   * the workflows.
   * @return a list of string
   */
  public List<String> getListLogStatusWorkflow() {
    List<String> listLog = new ArrayList<String>();
    for (String element : this.listWorkflowStatusLog) {
      listLog.add(element);
    }
    return listLog;
  }

  /**
   * This method of the class Fast5toFastq get a list of corrupt fast5 files
   * @return a list of string
   */
  public List<String> getListCorruptFileLog() {
    List<String> listCorruptFile = new ArrayList<String>();
    for (File file : this.listCorruptFast5Files) {
      listCorruptFile.add(file.toString());
    }
    return listCorruptFile;
  }

  /**
   * This method of the class Fast5toFastq get All Log for the conversion of
   * fast5 to fastq file.
   * @param status a string of the folder fast5
   * @param logGroup a LocalReporter
   */
  private void getAllListForLog(String status, LocalReporter logGroup) {
    //
    // Count of write sequence in fastq file
    //
    this.listWriteSequenceLog.add("In the file "
        + status + " template the number of total sequence write "
        + logGroup.getCounterValue("numberSequenceWrite",
            "numberSequenceTemplateWrite"));
    this.listWriteSequenceLog.add("In the file "
        + status + " complement the number of total sequence write "
        + logGroup.getCounterValue("numberSequenceWrite",
            "numberSequenceComplementWrite"));
    this.listWriteSequenceLog.add("In the file "
        + status + " barcode the number of total sequence write "
        + logGroup.getCounterValue("numberSequenceWrite",
            "numberSequenceBarcodeWrite"));

    //
    // Workflow status getters
    //
    this.listWorkflowStatusLog.add("\n###################################\n"
        + status + " :\n###################################");

    // Event Detection Workflow
    //
    this.listWorkflowStatusLog.add("\nEvent Detection workflow :\n");
    stackHashWorkflow(logGroup, "eventDetectionWorkflow", status,
        "Event Detection");

    // Hairpin Split Workflow
    //
    this.listWorkflowStatusLog.add("\nHairpin Split workflow :\n");
    stackHashWorkflow(logGroup, "hairpinSplitWorkflow", status,
        "Hairpin Split");

    // Basecall_1D Workflow
    //
    this.listWorkflowStatusLog.add("\nBasecall_1D workflow :\n");
    stackHashWorkflow(logGroup, "basecall1DWorkflow", status, "Basecall1D");

    // Basecall_2D Workflow
    //
    this.listWorkflowStatusLog.add("\nBasecall_2D workflow :\n");
    stackHashWorkflow(logGroup, "basecall2DWorkflow", status, "Basecall2D");

    // Calibration Strand Workflow
    //
    this.listWorkflowStatusLog.add("\nCalibration Strand :\n");
    stackHashWorkflow(logGroup, "calibrationStrandWorkflow", status,
        "Calibration Strand");

    // Barcode Workflow
    //
    this.listWorkflowStatusLog.add("\nBarcode workflow :\n");
    stackHashWorkflow(logGroup, "barcodeWorkflow", status, "Barcode");

  }

  /**
   * This method of the class Fast5toFastq stack the list of workflow status
   * with information of workflow status hash.
   * @param statusHash, a hash of final status
   * @param status, the name of folder fast5 status
   * @param nameWorkflow, the name of the workflow
   */
  private void stackHashWorkflow(LocalReporter logGroup, String group,
      String status, String nameWorkflow) {
    Set<String> keys = logGroup.getCounterNames(group);

    if (keys.size() >= 1) {
      for (String key : keys) {
        this.listWorkflowStatusLog.add("The status \""
            + key + "\" is present " + logGroup.getCounterValue(group, key)
            + " times in the folder " + status);
      }
    } else {
      this.listWorkflowStatusLog.add("No status for the "
          + nameWorkflow + " Workflow in the folder " + status);
    }
  }

  //
  //
  // Important methods
  //
  //

  /**
   * This method of the class Fast5toFastq launch processDirectory for a
   * directory of fast5.
   * @param fast5SubdirName is a directory of fast5 files
   * @param status is the status of fast5 file
   * @return an int who is the number of fast5 process
   * @throws Exception
   */
  private int processDirectory(String fast5SubdirName, String status,
      LocalReporter logGroup) throws Exception {

    List<File> list = listFast5(fast5SubdirName);
    processDirectory(list, status, logGroup);
    return list.size();
  }

  /**
   * This method of the class Fast5toFastq process the type of sequence and
   * launch the read of fast5 and write of fastq sequence.
   * @param listFast5Files is a list of fast5 file
   * @param status is the status of fast5 file
   * @throws Exception
   */
  private void processDirectory(List<File> listFast5Files, String status,
      LocalReporter logGroup) throws Exception {

    if (listFast5Files.isEmpty()) {
      return;
    }
    // Create writters
    Writer complementWriter = null;
    Writer templateWriter = null;
    Writer barcodeWriter = null;

    if (this.saveComplementSequence) {
      complementWriter =
          createWriterFastq(listFast5Files.get(0), "complement", status);
    }
    if (this.saveTemplateSequence) {
      templateWriter =
          createWriterFastq(listFast5Files.get(0), "template", status);
    }
    if (this.saveBarcodeSequence) {
      barcodeWriter =
          createWriterFastq(listFast5Files.get(0), "barcode", status);
    }

    // Read all Fast5 files

    long start1 = System.currentTimeMillis();

    readFast5WriteFastq(listFast5Files, complementWriter, templateWriter,
        barcodeWriter, status, logGroup);

    long end1 = System.currentTimeMillis();
    System.out.println("Time exe 1 thread:"
        + (end1 - start1) / 1000 + "s for a " + listFast5Files.size()
        + " number of fast5");

    // for (int i = 0; i <= listFast5Files.size(); i += 10000) {
    //
    // long start2 = System.currentTimeMillis();
    // multiThreadReadFast5WriteFastq(listFast5Files.subList(0, i),
    // complementWriter, templateWriter, barcodeWriter, status, logGroup);
    //
    // long end2 = System.currentTimeMillis();
    //
    // System.out.println("Time exe multi thread :"
    // + (end2 - start2) / 1000 + "s for a " + i + " number of fast5");
    // }
    // long start3 = System.currentTimeMillis();
    // multiThreadReadFast5WriteFastq(listFast5Files, complementWriter,
    // templateWriter, barcodeWriter, status, logGroup);
    //
    // long end3 = System.currentTimeMillis();
    //
    // System.out
    // .println("Time exe multi thread :" + (end3 - start3) / 1000 + "s");

    // Close writters
    if (this.saveComplementSequence) {
      complementWriter.close();
    }
    if (this.saveTemplateSequence) {
      templateWriter.close();
    }
    if (this.saveBarcodeSequence) {
      barcodeWriter.close();
    }
  }

  /**
   * This method of the class Fast5toFastq execute processDirectory with a list
   * of barcode.
   * @param listBarcodeDir is the list of barcode of the run
   * @throws Exception
   */
  private void processDirectories(List<File> listBarcodeDir,
      LocalReporter logGroup) throws Exception {
    for (File barcodeDirectory : listBarcodeDir) {
      List<File> listBarcodeFast5Files = listFast5(barcodeDirectory);
      processDirectory(listBarcodeFast5Files, barcodeDirectory.getName(),
          logGroup);
      this.numberBarcodeFast5Files += listBarcodeFast5Files.size();
      this.numberFast5Files += listBarcodeFast5Files.size();
    }
  }

  private void readFast5WriteFastq(File fast5File, Writer complementWriter,
      Writer templateWriter, Writer barcodeWriter, String status,
      LocalReporter logGroup) throws IOException {

    try (Fast5 f5 = new Fast5(fast5File)) {
      if (complementWriter != null && f5.getComplementFastq() != null) {
        complementWriter.write(f5.getComplementFastq());
        logGroup.incrCounter("numberSequenceWrite",
            "numberSequenceComplementWrite", 1);
      }
      if (templateWriter != null && f5.getTemplateFastq() != null) {
        templateWriter.write(f5.getTemplateFastq());
        logGroup.incrCounter("numberSequenceWrite",
            "numberSequenceTemplateWrite", 1);
      }
      if (barcodeWriter != null && f5.getLongBarcodingFastq() != null) {
        barcodeWriter.write(f5.getLongBarcodingFastq());
        logGroup.incrCounter("numberSequenceWrite",
            "numberSequenceBarcodeWrite", 1);
      }

      //
      // Get the Workflows Informations in the SynchronizedCountWriteFastq
      // Object
      //

      // Barcode Workflow
      //
      if (logGroup.getCounterNames("barcodeWorkflow")
          .contains(f5.getBarcodindFinalStatus())) {
        logGroup.incrCounter("barcodeWorkflow", f5.getBarcodindFinalStatus(),
            1);
      } else {
        logGroup.setCounter("barcodeWorkflow", f5.getBarcodindFinalStatus(),
            1);
      }

      // Basecall_1D Workflow
      //
      if (logGroup.getCounterNames("basecall1DWorkflow")
          .contains(f5.getBaseCall1DFinalStatus())) {
        logGroup.incrCounter("basecall1DWorkflow",
            f5.getBaseCall1DFinalStatus(), 1);
      } else {
        logGroup.setCounter("basecall1DWorkflow", f5.getBaseCall1DFinalStatus(),
            1);
      }

      // Basecall_2D Workflow
      //
      if (logGroup.getCounterNames("basecall2DWorkflow")
          .contains(f5.getBaseCall2DFinalStatus())) {
        logGroup.incrCounter("basecall2DWorkflow",
            f5.getBaseCall2DFinalStatus(), 1);
      } else {
        logGroup.setCounter("basecall2DWorkflow", f5.getBaseCall2DFinalStatus(),
            1);
      }

      // Calibration Strand Workflow
      //
      if (logGroup.getCounterNames("calibrationStrandWorkflow")
          .contains(f5.getCalibrationStrandFinalStatus())) {
        logGroup.incrCounter("calibrationStrandWorkflow",
            f5.getCalibrationStrandFinalStatus(), 1);
      } else {
        logGroup.setCounter("calibrationStrandWorkflow",
            f5.getCalibrationStrandFinalStatus(), 1);
      }

      // Event Detection Workflow
      //
      if (logGroup.getCounterNames("eventDetectionWorkflow")
          .contains(f5.getEventDetectionFinalStatus())) {
        logGroup.incrCounter("eventDetectionWorkflow",
            f5.getEventDetectionFinalStatus(), 1);
      } else {
        logGroup.setCounter("eventDetectionWorkflow",
            f5.getEventDetectionFinalStatus(), 1);
      }

      // Hairpin Split Workflow
      //
      if (logGroup.getCounterNames("hairpinSplitWorkflow")
          .contains(f5.getHairpinSplitFinalStatus())) {
        logGroup.incrCounter("hairpinSplitWorkflow",
            f5.getHairpinSplitFinalStatus(), 1);
      } else {
        logGroup.setCounter("hairpinSplitWorkflow",
            f5.getHairpinSplitFinalStatus(), 1);
      }

      //
      //
      //

    } catch (HDF5Exception e) {
      this.numberCorruptFast5Files++;
      this.listCorruptFast5Files.add(fast5File);
    }
  }

  private void multiThreadReadFast5WriteFastq(List<File> listFast5Files,
      final Writer complementWriter, final Writer templateWriter,
      final Writer barcodeWriter, final String status,
      final LocalReporter logGroup) throws IOException {

    ExecutorService executor = Executors
        .newFixedThreadPool(Runtime.getRuntime().availableProcessors());

    for (final File fast5File : listFast5Files) {
      executor.submit(new Runnable() {

        @Override
        public void run() {
          try {
            readFast5WriteFastq(fast5File, complementWriter, templateWriter,
                barcodeWriter, status, logGroup);
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
      ;
    }
    getAllListForLog(status, logGroup);

  }

  /**
   * This method of the class Fast5toFastq read fast5 files on a list and write
   * the fastq sequence.
   * @param listFast5Files is the list of fast5 file
   * @param complementWriter is the writer of the complement sequence
   * @param templateWriter is the writer of the template sequence
   * @param barcodeWriter is the writer of the barcode sequence
   * @param status is the status of the fast5 file
   * @throws IOException
   */
  private void readFast5WriteFastq(List<File> listFast5Files,
      Writer complementWriter, Writer templateWriter, Writer barcodeWriter,
      String status, LocalReporter logGroup) throws IOException {

    for (File fast5File : listFast5Files) {
      readFast5WriteFastq(fast5File, complementWriter, templateWriter,
          barcodeWriter, status, logGroup);
    }

    getAllListForLog(status, logGroup);

  }

  //
  //
  // Main function
  //
  //

  /**
   * This method of the class Fast5toFastq execute the process to retrieve the
   * fastq sequence on fast5 file.
   * @throws Exception
   */
  public void execute() throws Exception {

    final LocalReporter logGroup = new LocalReporter();

    if (this.processMergeStatus) {
      this.listFast5Files = listAllFast5();
      processDirectory(this.listFast5Files, "merge_status", logGroup);
      return;
    }

    if (this.processFail) {
      this.numberFailFast5Files =
          processDirectory("downloads/fail", "fail", logGroup);
      this.numberFast5Files += this.numberFailFast5Files;
      logGroup.clear();
    }
    if (this.processPass) {
      this.numberPassFast5Files =
          processDirectory("downloads/pass", "pass", logGroup);
      this.numberFast5Files += this.numberPassFast5Files;
      logGroup.clear();
    }
    if (this.processFailBarcode) {
      this.numberFailBarcodeFast5Files = processDirectory(
          "downloads/fail/unclassified", "unclassified", logGroup);
      this.numberFast5Files += this.numberFailBarcodeFast5Files;
      logGroup.clear();
    }
    if (this.processPassBarcode) {
      List<File> listBarcodeFast5Dir =
          listSubDir(new File(this.fast5RunDirectory, "downloads/pass"));
      processDirectories(listBarcodeFast5Dir, logGroup);
      logGroup.clear();
    }
  }
}
