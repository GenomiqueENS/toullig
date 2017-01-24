package fr.ens.biologie.genomique.nanopore;

import java.io.File;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.Writer;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import ncsa.hdf.hdf5lib.exceptions.HDF5Exception;

/**
 * This class read a minION run of Fast5 basecalled to extract fastq sequence.
 * @author Aurelien Birer
 */
public class Fast5toFastq {

  private File fast5RunDirectory;
  private File repertoryFastqOutput;

  private List<File> listFast5Files = new ArrayList<File>();

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
  private boolean saveHairpinSequence = false;
  private boolean saveBarcodeSequence = false;

  private List<String> listWriteSequenceLog = new ArrayList<String>();

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

    COMPLEMENT("complement"), TEMPLATE("template"), HAIRPIN("hairpin"),
    BARCODE("barcode");

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
   * @throws IOException
   */
  private Writer createWriterFastq(File fast5File, String typeSequence,
      String status) throws IOException {
    // the substring is done on the string "_ch" who correspond to the channel
    // number
    String preNameFile = fast5File.getName().substring(0,
        fast5File.getName().indexOf("_ch") + 1);
    Writer fastqFile = new FileWriter(new File(this.repertoryFastqOutput
        + "/" + preNameFile + status + "_" + typeSequence + ".fastq"));
    return fastqFile;
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
   * This method of the class Fast5toFastq set the type of hairpin sequence on
   * process.
   * @param saveHairpinSequence, a boolean to process the type of hairpin
   *          sequences
   */
  public void setSaveHairpinSequence(boolean saveHairpinSequence) {
    if (saveHairpinSequence) {
      this.saveHairpinSequence = true;
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
  //
  // Create Log in a Hash
  //
  //

  /**
   * This method of the class Fast5toFastq get a list of log on the number of
   * files.
   * @return a hash
   */
  public List<String> getListLog() {
    List<String> listLog = new ArrayList<String>();
    listLog.add("Number of fast5 files read " + getNumberFast5Files());
    listLog.add("Number of fast5 files read " + getNumberFast5Files());
    listLog.add("Number of corrupt files " + getNumberCorruptFast5Files());
    listLog.add("Number of fail files read " + getNumberFailFast5Files());
    listLog.add("Number of pass files read " + getNumberPassFast5Files());
    listLog.add("Number of fail attribution barcode files read "
        + getNumberBadBarcodeFast5Files());
    listLog
        .add("Number of Barcode file " + getNumberBarcodeFast5Files() + "\n");
    for (String element : this.listWriteSequenceLog) {
      listLog.add(element);
    }
    return listLog;
  }

  //
  //
  // Important methods
  //
  //

  /**
   * This method of the class Fast5toFastq retrun the name of the directory of the minION run.
   * @return a string of the name of root directory
   */
  public String getNameDirectoryRunFast5() {
    return this.fast5RunDirectory.getName();
  }

  /**
   * This method of the class Fast5toFastq launch processDirectory for a directory of fast5.
   * @param fast5SubdirName is a directory of fast5 files
   * @param status is the status of fast5 file
   * @return an int who is the number of fast5 process
   * @throws IOException
   */
  private int processDirectory(String fast5SubdirName, String status)
      throws IOException {

    List<File> list = listFast5(fast5SubdirName);
    processDirectory(list, status);
    return list.size();
  }

  /**
   * This method of the class Fast5toFastq process the type of sequence and
   * launch the read of fast5 and write of fastq sequence.
   * @param listFast5Files is a list of fast5 file
   * @param status is the status of fast5 file
   * @throws IOException
   */
  private void processDirectory(List<File> listFast5Files, String status)
      throws IOException {

    if (listFast5Files.isEmpty()) {
      return;
    }
    // Create writters
    Writer complementWriter = null;
    Writer templateWriter = null;
    Writer hairpinWriter = null;
    Writer barcodeWriter = null;

    if (this.saveComplementSequence) {
      complementWriter =
          createWriterFastq(listFast5Files.get(0), "complement", status);
    }
    if (this.saveTemplateSequence) {
      templateWriter =
          createWriterFastq(listFast5Files.get(0), "template", status);
    }
    if (this.saveHairpinSequence) {
      hairpinWriter =
          createWriterFastq(listFast5Files.get(0), "hairpin", status);
    }
    if (this.saveBarcodeSequence) {
      barcodeWriter =
          createWriterFastq(listFast5Files.get(0), "barcode", status);
    }
    // Read all Fast5 files
    readFast5WriteFastq(listFast5Files, complementWriter, templateWriter,
        hairpinWriter, barcodeWriter, status);

    // Close writters
    if (this.saveComplementSequence) {
      complementWriter.close();
    }
    if (this.saveTemplateSequence) {
      templateWriter.close();
    }
    if (this.saveHairpinSequence) {
      hairpinWriter.close();
    }
    if (this.saveBarcodeSequence) {
      barcodeWriter.close();
    }
  }

  /**
   * This method of the class Fast5toFastq execute processDirectory with a list
   * of barcode.
   * @param listBarcodeDir is the list of barcode of the run
   * @throws IOException
   */
  private void processDirectories(List<File> listBarcodeDir)
      throws IOException {
    for (File barcodeDirectory : listBarcodeDir) {
      List<File> listBarcodeFast5Files = listFast5(barcodeDirectory);
      processDirectory(listBarcodeFast5Files, barcodeDirectory.getName());
      this.numberBarcodeFast5Files += listBarcodeFast5Files.size();
    }
  }

  /**
   * This method of the class Fast5toFastq read fast5 files on a list and write
   * the fastq sequence.
   * @param listFast5Files is the list of fast5 file
   * @param complementWriter is the writer of the complement sequence
   * @param templateWriter is the writer of the template sequence
   * @param hairpinWriter is the writer of the hairpin sequence
   * @param barcodeWriter is the writer of the barcode sequence
   * @param status is the status of the fast5 file
   * @throws IOException
   */
  private void readFast5WriteFastq(List<File> listFast5Files,
      Writer complementWriter, Writer templateWriter, Writer hairpinWriter,
      Writer barcodeWriter, String status) throws IOException {
    int numberSequenceComplementWrite = 0;
    int numberSequenceTemplateWrite = 0;
    int numberSequenceHairpinWrite = 0;
    int numberSequenceBarcodeWrite = 0;
    for (File fast5File : listFast5Files) {
      try (Fast5 f5 = new Fast5(fast5File)) {
        if (complementWriter != null && f5.getComplementFastq() != null) {
          complementWriter.write(f5.getComplementFastq());
          numberSequenceComplementWrite++;
        }
        if (templateWriter != null && f5.getTemplateFastq() != null) {
          templateWriter.write(f5.getTemplateFastq());
          numberSequenceTemplateWrite++;
        }
        if (hairpinWriter != null && f5.getHairpinFastq() != null) {
          hairpinWriter.write(f5.getHairpinFastq());
          numberSequenceHairpinWrite++;
        }
        if (barcodeWriter != null && f5.getBarcodingFastq() != null) {
          barcodeWriter.write(f5.getBarcodingFastq());
          numberSequenceBarcodeWrite++;
        }
      } catch (HDF5Exception e) {
        this.numberCorruptFast5Files++;
      }
    }
    this.listWriteSequenceLog.add("In the file "
        + status + " complement the number of total sequence write "
        + numberSequenceComplementWrite);
    this.listWriteSequenceLog.add("In the file "
        + status + " template the number of total sequence write "
        + numberSequenceTemplateWrite);
    this.listWriteSequenceLog.add("In the file "
        + status + " hairpin the number of total sequence write "
        + numberSequenceHairpinWrite);
    this.listWriteSequenceLog.add("In the file "
        + status + " barcode the number of total sequence write "
        + numberSequenceBarcodeWrite);
  }

  //
  //
  // Main function
  //
  //

  /**
   * This method of the class Fast5toFastq execute the process to retrieve the
   * fastq sequence on fast5 file.
   * @throws IOException
   */
  public void execution() throws IOException {

    if (this.processMergeStatus) {

      this.listFast5Files = listAllFast5();
      processDirectory(this.listFast5Files, "merge_status");
      this.numberFast5Files = this.listFast5Files.size();
      return;
    }

    if (this.processFail) {
      this.numberFailFast5Files = processDirectory("downloads/fail", "fail");
    }
    if (this.processPass) {
      this.numberPassFast5Files = processDirectory("downloads/pass", "pass");
    }
    if (this.processFailBarcode) {
      this.numberFailBarcodeFast5Files =
          processDirectory("downloads/fail/unclassified", "unclassified");
    }
    if (this.processPassBarcode) {
      List<File> listBarcodeFast5Dir =
          listSubDir(new File(this.fast5RunDirectory, "downloads/pass"));
      processDirectories(listBarcodeFast5Dir);
    }
  }
}
