package fr.ens.biologie.genomique.nanopore;

import java.io.File;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.Writer;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;

import ncsa.hdf.hdf5lib.exceptions.HDF5Exception;

/**
 * @author Aurelien Birer
 */
public class Fast5toFastq {

  private static final class OutputFile {

    private final File file;
    private final TypeFastq type;

    @Override
    public int hashCode() {
      return Objects.hash(file, type);
    }

    @Override
    public boolean equals(Object obj) {
      if (this == obj) {
        return true;
      }
      if (!(obj instanceof OutputFile)) {
        return false;
      }
      OutputFile that = (OutputFile) obj;
      return Objects.equals(this.file, that.file) && this.type == that.type;
    }

    /**
     * The constructor of the class outputFile of the class Fast5toFastq.
     * @param file is the name of the fastq writer
     * @param type is the sequence type
     */
    public OutputFile(File file, TypeFastq type) {
      this.file = file;
      this.type = type;
    }

    /**
     * This method of the class outputFile of the class Fast5toFastq get the
     * file to write fastq sequences.
     * @return a file of the fastq writer
     */
    public File getFile() {
      return this.file;
    }

    /**
     * This method of the class outputFile of the class Fast5toFastq get the
     * type of the sequence
     * @return a typefastq of the sequence
     */
    public TypeFastq getType() {
      return this.type;
    }
  }

  private File rootRepertoryFast5Run;
  private File repertoryFastqOutput;

  private List<File> listFast5Files = new ArrayList<File>();
  private List<File> listFailFast5Files = new ArrayList<File>();
  private List<File> listPassFast5Files = new ArrayList<File>();
  private List<File> listCorruptFast5Files = new ArrayList<File>();
  private List<File> listBarcodeFast5Dir = new ArrayList<File>();
  private List<File> listBarcodeFast5Files = new ArrayList<File>();
  private List<File> listBadBarcodeFast5Files = new ArrayList<File>();
  private List<String> listBarcode = new ArrayList<String>();

  private boolean processStatus = false;
  private boolean processFail = false;
  private boolean processPass = false;
  private boolean processFailBarcode = false;
  private boolean processPassBarcode = false;

  private boolean saveComplementSequence = false;
  private boolean saveTemplateSequence = false;
  private boolean saveHairpinSequence = false;
  private boolean saveBarcodeSequence = false;

  private Map<OutputFile, List<File>> hashProcess =
      new HashMap<OutputFile, List<File>>();

  private Map<String, Integer> hashWriteSequenceLog =
      new HashMap<String, Integer>();
  

  //
  //
  // The variable rootRepertoryFast5Run contains the folder downloads with the
  // basecalled reads and the folder uploaded with the prebasecalling reads.
  //
  //
  /**
   * The constructor of the class Fast5toFastq.
   * @param rootRepertoryFast5Run is the repertory of the ONT run
   * @param repertoryFastqOutput is the repertory to store fastq sequence
   * @throws IOException
   */
  public Fast5toFastq(File rootRepertoryFast5Run, File repertoryFastqOutput)
      throws IOException {

    if (!rootRepertoryFast5Run.exists()) {
      throw new IOException(
          "The repertory " + rootRepertoryFast5Run + " dont exist!");
    } else {
      this.rootRepertoryFast5Run = rootRepertoryFast5Run;
    }
    if (!repertoryFastqOutput.exists()) {
      throw new IOException(
          "The repertory " + repertoryFastqOutput + " dont exist!");
    } else {
      this.repertoryFastqOutput = repertoryFastqOutput;
    }
  }

  /**
   * Values of the object TypeFastq that design the sequence type.
   * @author birer
   */
  public enum TypeFastq {

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
    TypeFastq(String name) {
      this.name = name;
    }
  }

  //
  //
  // Create Hash Fast5 Files
  //
  //

  /**
   * This method of the class Fast5toFastq input value in the hashProcess in
   * function of the sequence type.
   */
  public void createHashFast5Files() {
    if (this.saveComplementSequence) {
      TypeFastq type = TypeFastq.COMPLEMENT;
      createHashNotBarcodedFast5Files(type);
      if (this.processPassBarcode) {
        createHashBarcodedFast5Files(type);
      }
    }
    if (this.saveTemplateSequence) {
      TypeFastq type = TypeFastq.TEMPLATE;
      createHashNotBarcodedFast5Files(type);
      if (this.processPassBarcode) {
        createHashBarcodedFast5Files(type);
      }
    }
    if (this.saveHairpinSequence) {
      TypeFastq type = TypeFastq.HAIRPIN;
      createHashNotBarcodedFast5Files(type);
      if (this.processPassBarcode) {
        createHashBarcodedFast5Files(type);
      }
    }
    if (this.saveBarcodeSequence) {
      TypeFastq type = TypeFastq.BARCODE;
      createHashNotBarcodedFast5Files(type);
      if (this.processPassBarcode) {
        createHashBarcodedFast5Files(type);
      }
    }
  }

  /**
   * This method of the class Fast5toFastq input value in the hashProcess of
   * barcoded fast5 files.
   * @param type is the type of sequence
   */
  public void createHashBarcodedFast5Files(TypeFastq type) {
    for (File directory : this.listBarcodeFast5Dir) {
      List<File> listBarcodeFast5File = listFast5(directory);
      this.listBarcodeFast5Files.addAll(listBarcodeFast5File);
      this.listFast5Files.addAll(listBarcodeFast5File);
      String name = directory.getName() + "_" + type.getName();
      OutputFile FastqFile = new OutputFile(
          createFastqFile(listBarcodeFast5File.get(0), name), type);
      hashProcess.put(FastqFile, listBarcodeFast5File);
    }
  }

  /**
   * This method of the class Fast5toFastq input value in the hashProcess of not
   * barcoded fast5 files.
   * @param type is the type of sequence
   */
  public void createHashNotBarcodedFast5Files(TypeFastq type) {

    String name = "";
    if (this.processStatus) {
      if (this.processFail) {
        name = this.listFailFast5Files.get(0).getParentFile().getName()
            + "_" + type.getName();
        OutputFile FastqFile = new OutputFile(
            createFastqFile(this.listFailFast5Files.get(0), name), type);
        hashProcess.put(FastqFile, this.listFailFast5Files);
      }
      if (this.processPass) {
        name = this.listPassFast5Files.get(0).getParentFile().getName()
            + "_" + type.getName();
        OutputFile FastqFile = new OutputFile(
            createFastqFile(this.listPassFast5Files.get(0), name), type);
        hashProcess.put(FastqFile, this.listPassFast5Files);
      }
      if (this.processFailBarcode) {
        name = this.listBadBarcodeFast5Files.get(0).getParentFile().getName()
            + "_" + type.getName();
        OutputFile FastqFile = new OutputFile(
            createFastqFile(this.listBadBarcodeFast5Files.get(0), name), type);
        hashProcess.put(FastqFile, this.listBadBarcodeFast5Files);
      }
    } else {
      name = type.getName();
      OutputFile FastqFile = new OutputFile(
          createFastqFile(this.listFast5Files.get(0), name), type);
      hashProcess.put(FastqFile, this.listFast5Files);
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

  /**
   * This method of the class Fast5toFastq list the fast5 files of a directory.
   * @param fast5Dir is a directory (must be contains fast5 files)
   * @return a list of fast5 file
   */
  private List<File> listFast5(File fast5Dir) {
    return Arrays.asList(fast5Dir.listFiles(new FilenameFilter() {

      @Override
      public boolean accept(File dir, String name) {
        return name.endsWith(".fast5");
      }
    }));
  }

  //
  //
  // Writer Fastq
  //
  //

  /**
   * This method of the class Fast5toFastq create the fastq output file.
   * @param fastqFile is the path to store the fastq files
   * @param hashProcess is the hash of all information for processing the write
   *          of fastq (name of output file + type of sequence + list of fast5
   *          file)
   * @throws IOException if the fast5 file is corrupt
   */
  private void createFastq(File fastqFile,
      Map<OutputFile, List<File>> hashProcess) throws IOException {
    for (OutputFile outputFile : hashProcess.keySet()) {
      try {
        Writer writer = new FileWriter(
            new File(fastqFile, outputFile.getFile().toString()));
        String seq = null;
        int numberSequenceWrite = 0;
        for (File file : hashProcess.get(outputFile)) {
          try (Fast5 f5 = new Fast5(file)) {
            switch (outputFile.getType()) {
            case COMPLEMENT:
              seq = f5.getComplementFastq();
              break;
            case TEMPLATE:
              seq = f5.getTemplateFastq();
              break;
            case HAIRPIN:
              seq = f5.getHairpinFastq();
              break;
            case BARCODE:
              seq = f5.getBarcodingFastq();
              break;
            default:
              break;
            }
            f5.close();
          } catch (HDF5Exception e) {
            // e.printStackTrace();
            this.listCorruptFast5Files.add(file);
          }
          if (seq != null) {
            numberSequenceWrite++;
            writer.write(seq);
          }
        }
        this.hashWriteSequenceLog.put("In the file "
            + outputFile.getFile().toString()
            + " the number of total sequence write", numberSequenceWrite);
        writer.close();
      } catch (IOException e) {
        e.printStackTrace();
      }
    }
  }

  //
  //
  // Create Writer Fastq Name
  //
  //

  /**
   * This method of the class Fast5toFastq create the good name of the fastq
   * output file.
   * @param file is the name of the first file of the list "listFast5Files"
   * @param typeSequence is the type of sequence (ex:complement)
   * @return a file with the correct output name for write a fastq sequence
   */
  public File createFastqFile(File file, String typeSequence) {
    String preNameFile =
        file.getName().substring(0, file.getName().indexOf("_ch") + 1);// the
    // substring
    // is
    // done
    // on
    // the
    // string
    // "_ch"
    // who
    // correspond
    // to
    // the
    // channel
    // number
    File fastqFile = new File(preNameFile + typeSequence + ".fastq");
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
    return this.listFast5Files.size();
  }

  /**
   * This method of the class Fast5toFastq get the number of pass files.
   * @return an int of the number of pass files
   */
  public int getNumberPassFast5Files() {
    return this.listPassFast5Files.size();
  }

  /**
   * This method of the class Fast5toFastq get the number of corrupt files.
   * @return an int of the number of corrupt files
   */
  public int getNumberCorruptFast5Files() {
    return this.listCorruptFast5Files.size();
  }

  /**
   * This method of the class Fast5toFastq get the number of fail files.
   * @return an int of the number of fail files
   */
  public int getNumberFailFast5Files() {
    return this.listFailFast5Files.size();
  }

  /**
   * This method of the class Fast5toFastq get the number of bad barcoded files.
   * @return an int of the number of bad barcoded files
   */
  public int getNumberBadBarcodeFast5Files() {
    return this.listBadBarcodeFast5Files.size();
  }

  /**
   * This method of the class Fast5toFastq get the number of barcoded files.
   * @return an int of the number of barcoded files
   */
  public int getNumberBarcodeFast5Files() {
    return this.listBarcodeFast5Files.size();
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
  public void setStatusDifferenciate(boolean processStatus) {
    if (processStatus) {
      this.processStatus = true;
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
  // Getter List Files
  //
  //

  /**
   * This method of the class Fast5toFastq obtain the list of fast5 files.
   * @return a list of fast5 files
   */
  public List<File> getListFast5Files() {
    return this.listFast5Files;
  }

  /**
   * This method of the class Fast5toFastq obtain the list of fail fast5 files.
   * @return a list of fail fast5 files
   */
  public List<File> getListFailFast5Files() {
    return listFailFast5Files;
  }

  /**
   * This method of the class Fast5toFastq obtain the list of pass fast5 files.
   * @return a list of pass fast5 files
   */
  public List<File> getListPassFast5Files() {
    return listPassFast5Files;
  }

  /**
   * This method of the class Fast5toFastq obtain the list of corrupt fast5
   * files.
   * @return a list of corrupt fast5 files
   */
  public List<File> getListCorruptFast5Files() {
    return this.listCorruptFast5Files;
  }

  /**
   * This method of the class Fast5toFastq obtain the list of barcoded fast5
   * files.
   * @return a list of barcoded fast5 files
   */
  public List<File> getListBarcodeFast5Files() {
    return this.listBarcodeFast5Dir;
  }

  /**
   * This method of the class Fast5toFastq obtain the list of bad barcoded fast5
   * files.
   * @return a list of bad barcoded fast5 files
   */
  public List<File> getListBadBarcodedFast5Files() {
    return this.listBadBarcodeFast5Files;
  }

  /**
   * This method of the class Fast5toFastq obtain the list of barcode.
   * @return a list of bardcode
   */
  public List<String> getListBarcode() {
    return this.listBarcode;
  }

  //
  //
  // Create Log in a Hash
  //
  //

  /**
   * This method of the class Fast5toFastq get a hash of log on the number of files.
   * @return a hash
   */
  public Map<String, Integer> getHashLog() {
    Map<String, Integer> hashLog = new HashMap<>();
    hashLog.put("Number of fast5 files read", getNumberFast5Files());
    hashLog.put("Number of corrupt files", getNumberCorruptFast5Files());
    hashLog.put("Number of fail files read", getNumberFailFast5Files());
    hashLog.put("Number of pass files read", getNumberPassFast5Files());
    hashLog.put("Number of fail attribution barcode files read",
        getNumberBadBarcodeFast5Files());
    hashLog.put("Number of Barcode file", getNumberBarcodeFast5Files());
    for(String key : this.hashWriteSequenceLog.keySet()){
      hashLog.put(key, this.hashWriteSequenceLog.get(key));
    }
    return hashLog;
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

    if (this.processFail) {
      this.listFailFast5Files.addAll(
          listFast5(new File(this.rootRepertoryFast5Run, "downloads/fail")));
      if (listFailFast5Files.isEmpty()) {
        throw new IOException("The repertory "
            + rootRepertoryFast5Run + "downloads/fail is empty!");
      }
      this.listFast5Files.addAll(this.listFailFast5Files);
    }
    if (this.processPass) {
      this.listPassFast5Files.addAll(
          listFast5(new File(this.rootRepertoryFast5Run, "downloads/pass")));
      if (listPassFast5Files.isEmpty()) {
        throw new IOException("The repertory "
            + rootRepertoryFast5Run + "downloads/pass is empty!");
      }
      this.listFast5Files.addAll(this.listPassFast5Files);
    }
    if (this.processFailBarcode) {
      this.listBadBarcodeFast5Files.addAll(listFast5(
          new File(this.rootRepertoryFast5Run, "downloads/fail/unclassified")));
      if (listBadBarcodeFast5Files.isEmpty()) {
        throw new IOException("The repertory "
            + rootRepertoryFast5Run + "downloads/fail/unclassified is empty!");
      }
      this.listFast5Files.addAll(this.listBadBarcodeFast5Files);
    }
    if (this.processPassBarcode) {
      this.listBarcodeFast5Dir.addAll(
          listSubDir(new File(this.rootRepertoryFast5Run, "downloads/pass")));
      if (listBarcodeFast5Dir.isEmpty()) {
        throw new IOException("The sub-repertories "
            + rootRepertoryFast5Run + "downloads/pass are empty!");
      }
    }

    if (this.listFast5Files.isEmpty()) {
      throw new IOException("No files to process!");
    }

    createHashFast5Files(); // Create the hash to process

    createFastq(this.repertoryFastqOutput, this.hashProcess);// Exploit
    // all
    // fast5
    // files
    // of the
    // hash.

  }
}
