package fr.ens.biologie.genomique.toullig.fast5tofastq;

import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.nio.file.*;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import fr.ens.biologie.genomique.eoulsan.util.LocalReporter;

/**
 * This class read a minION run of Fast5 basecalled to extract fastq sequence.
 * @author Aurelien Birer
 */
public class Fast5ToFastq {

  private File fast5RunDirectory;
  private File repertoryFastqOutput;

  private DirectoryProcessor processor;

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
  private List<File> listSubDir(File dir, String extension) {

    // create new List for the results
    List<File> result = new ArrayList<>();

    // get the directory in the result list
    for (File file : dir.listFiles()) {

      // test if the file is a directory
      if (file.isDirectory() && file.toString().contains(extension)) {
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
        new File(this.fast5RunDirectory, "downloads"), "")) {

      // test if the directory is a directory
      if (directory.isDirectory()) {

        // list all sub directories
        List<File> subDirectories = listSubDir(directory, "");

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

    if (this.processor == null) {
      return Collections.emptyList();
    }

    return this.processor.getListCorruptFast5Files();
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
  private int processDirectory(File fast5SubdirName, String status,
      LocalReporter localReporter) throws IOException {

    boolean subDirectoryNumber4000 = false;

    if (fast5SubdirName.listFiles() != null) {
      for (File f1 : fast5SubdirName.listFiles()) {
        subDirectoryNumber4000 = f1.isDirectory();
        break;
      }
    }

    // for search in sub Directory
    if (subDirectoryNumber4000) {

      List<File> list = new ArrayList<>();

      // list subdir
      for (File f1 : fast5SubdirName.listFiles()) {

        list.addAll(listFast5(f1));
        if (list.isEmpty()) {
          return 0;
        }
      }
      // process this list of fast5 files
      this.processor.processDirectory(list, status, localReporter);
      return list.size();

    } else {

      List<File> list = listFast5(fast5SubdirName);

      // process this list of fast5 files
      this.processor.processDirectory(list, status, localReporter);
      return list.size();
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

    boolean subDirectoryNumber4000 = false;

    // process each barcode directories
    for (File barcodeDirectory : listBarcodeDir) {

      for (File f1 : barcodeDirectory.listFiles()) {
        subDirectoryNumber4000 = f1.isDirectory();
        break;
      }

      // for search in sub Directory
      if (subDirectoryNumber4000) {

        List<File> listBarcodeFast5Files = new ArrayList<>();

        // list subdir
        for (File f1 : barcodeDirectory.listFiles()) {

          // create a list of fast5 for a barcode
          listBarcodeFast5Files.addAll(listFast5(f1));

        }
        // process this list
        this.processor.processDirectory(listBarcodeFast5Files,
            barcodeDirectory.getName(), localReporter);

        // incremente numberPassFast5Files counter
        localReporter.incrCounter("numberFiles", "numberPassFast5Files",
            listBarcodeFast5Files.size());

        // incremente numberFast5Files counter
        localReporter.incrCounter("numberFiles", "numberFast5Files",
            listBarcodeFast5Files.size());

      } else {
        // create a list of fast5 for a barcode
        List<File> listBarcodeFast5Files = listFast5(barcodeDirectory);

        // process this list
        this.processor.processDirectory(listBarcodeFast5Files,
            barcodeDirectory.getName(), localReporter);

        // incremente numberPassFast5Files counter
        localReporter.incrCounter("numberFiles", "numberPassFast5Files",
            listBarcodeFast5Files.size());

        // incremente numberFast5Files counter
        localReporter.incrCounter("numberFiles", "numberFast5Files",
            listBarcodeFast5Files.size());
      }
    }
  }

  /**
   * This method of the class Fast5ToFastq get a fast5 file.
   * @return a fast5 file
   */
  private File getAFast5File() throws IOException {

    File dir = this.fast5RunDirectory;

    // create new List for the results
    List<File> result1 = new ArrayList<>();
    List<File> result2 = new ArrayList<>();
    List<File> result3 = new ArrayList<>();

    // get the directory in the result list
    for (File file : dir.listFiles()) {

      // test if the file is a directory
      if (file.isDirectory()) {
        result1.add(file);
      }

      // test if is it's a file and have an extension ".fast5"
      if (file.isFile() && file.toString().contains(".fast5")) {

        try (DirectoryStream<Path> stream =
            Files.newDirectoryStream(dir.toPath(), "*.{fast5}")) {

          for (Path entry : stream) {

            try (Fast5 f5 = new Fast5(entry.toFile())) {

              // test if the fast5 file is basecalled
              if (f5.isBasecalled()) {

                return entry.toFile();

              }

            } catch (Exception e) {
              e.printStackTrace();
            }
          }
        } catch (DirectoryIteratorException ex) {
          // I/O error encounted during the iteration, the cause is an
          // IOException
          throw ex.getCause();
        }
      }
    }

    //
    //
    //

    // get the directory in the result list
    for (File resultDirectory : result1) {

      boolean notBasecalled = false;

      for (File file : resultDirectory.listFiles()) {

        // test if the file is a directory
        if (file.isDirectory()) {
          result2.add(file);

        }

        // test if is it's a file and have an extension ".fast5"
        if (file.isFile() && file.toString().contains(".fast5")) {

          try (DirectoryStream<Path> stream =
              Files.newDirectoryStream(resultDirectory.toPath(), "*.{fast5}")) {

            if (notBasecalled == true) {
              break;
            }

            for (Path entry : stream) {

              try (Fast5 f5 = new Fast5(entry.toFile())) {

                // test if the fast5 file is basecalled
                if (f5.isBasecalled()) {

                  return entry.toFile();

                } else {
                  notBasecalled = true;
                  break;
                }

              } catch (Exception e) {
                e.printStackTrace();
              }
            }
          } catch (DirectoryIteratorException ex) {
            // I/O error encounted during the iteration, the cause is an
            // IOException
            throw ex.getCause();
          }
        }
      }
    }

    //
    //
    //

    // get the directory in the result list
    for (File resultDirectory : result2) {

      for (File file : resultDirectory.listFiles()) {

        // test if the file is a directory
        if (file.isDirectory()) {
          result3.add(file);
        }

        // test if is it's a file and have an extension ".fast5"
        if (file.isFile() && file.toString().contains(".fast5")) {

          try (DirectoryStream<Path> stream =
              Files.newDirectoryStream(resultDirectory.toPath(), "*.{fast5}")) {

            for (Path entry : stream) {

              try (Fast5 f5 = new Fast5(entry.toFile())) {

                // test if the fast5 file is basecalled
                if (f5.isBasecalled()) {

                  return entry.toFile();

                }

              } catch (Exception e) {
                e.printStackTrace();
              }
            }
          } catch (DirectoryIteratorException ex) {
            // I/O error encounted during the iteration, the cause is an
            // IOException
            throw ex.getCause();
          }
        }
      }
    }

    //
    //
    //

    // get the directory in the result list
    for (File resultDirectory : result3) {

      for (File file : resultDirectory.listFiles()) {

        // test if is it's a file and have an extension ".fast5"
        if (file.isFile() && file.toString().contains(".fast5")) {

          try (DirectoryStream<Path> stream =
              Files.newDirectoryStream(resultDirectory.toPath(), "*.{fast5}")) {

            for (Path entry : stream) {

              try (Fast5 f5 = new Fast5(entry.toFile())) {

                // test if the fast5 file is basecalled
                if (f5.isBasecalled()) {

                  return entry.toFile();

                }

              } catch (Exception e) {
                e.printStackTrace();
              }
            }
          } catch (DirectoryIteratorException ex) {
            // I/O error encounted during the iteration, the cause is an
            // IOException
            throw ex.getCause();
          }
        }
      }
    }

    return null;
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

    System.out.println("Sampling of a Fast5 file !");

    // get a fast5 file
    File sampleFast5File = getAFast5File();

    System.out.println("End of the sampling of a Fast5 file");

    try (Fast5 f5 = new Fast5(sampleFast5File)) {

      Fast5.Status status = f5.getStatus();
      Fast5.Basecaller basecaller = f5.getBasecaller();
      Fast5.Version version = f5.getVersion();
      Fast5.Type type = f5.getType();
      Fast5.ChemistryVersion chemistryVersion = f5.getChemistryVersion();

      this.processor = new DirectoryProcessor(repertoryFastqOutput,
          saveComplementSequence, saveTemplateSequence, saveConsensusSequence,
          saveTranscriptSequence, saveCompressGZIP, saveCompressBZIP2, status,
          basecaller, version, type, chemistryVersion,
          basecaller == Fast5.Basecaller.METRICHOR);

      if (basecaller == Fast5.Basecaller.METRICHOR || basecaller == null) {

        // execution for the basecaller Metrichor classification
        executeBasecallerMetrichor();
      }

      if (basecaller == Fast5.Basecaller.ALBACORE) {

        // execution for the basecaller Metrichor classification
        executeBasecallerAlbacore();
      }

    } catch (Exception e) {
      e.printStackTrace();
    }

    System.out.println("Conversion finish !");
  }

  /**
   * This method of the class Fast5ToFastq execute the process to retrieve the
   * fastq sequence on fast5 file for the basecaller Metrichor.
   * @throws IOException, test the read of the file
   */
  public void executeBasecallerMetrichor() throws IOException {

    // test if the merge of fastq is enable
    if (this.processMergeStatus) {

      // process all fast5
      this.processor.processDirectory(listAllFast5(), "merge_status",
          this.localReporter);
      return;
    }

    // test if the fail fast5 is to process
    if (this.processFail) {

      // get the list of fail fast5 files
      int numberFailFast5Files = processDirectory(
          new File(
              this.fast5RunDirectory.toPath().toString() + "/downloads/fail"),
          "fail", this.localReporter);

      // incremente fail fast5 files counter
      this.localReporter.incrCounter("numberFiles", "numberFailFast5Files",
          numberFailFast5Files);

      // incremente fast5 files counter
      this.localReporter.incrCounter("numberFiles", "numberFast5Files",
          numberFailFast5Files);
    }

    // test if the pass fast5 is to process
    if (this.processPass) {

      // get the list of barcode pass fast5 files
      List<File> listBarcodeFast5Dir =
          listSubDir(new File(this.fast5RunDirectory, "/downloads/pass"), "");

      // test if the number of pass fast5 files is not null
      if (listBarcodeFast5Dir.size() == 0) {

        // get the list of pass fast5 files
        int numberPassFast5Files = processDirectory(
            new File(
                this.fast5RunDirectory.toPath().toString() + "/downloads/pass"),
            "pass", this.localReporter);

        // incremente pass fast5 files counter
        this.localReporter.incrCounter("numberFiles", "numberPassFast5Files",
            numberPassFast5Files);

        // incremente fast5 files counter
        this.localReporter.incrCounter("numberFiles", "numberFast5Files",
            numberPassFast5Files);
      } else {

        // process the barcode directories
        processDirectories(listBarcodeFast5Dir, this.localReporter);
      }

    }

    // test if the unclassified fast5 is to process
    if (this.processUnclassified) {

      // get the list of unclassified fast5 files
      int numberUnclassifiedFast5Files = processDirectory(
          new File(this.fast5RunDirectory.toPath().toString()
              + "/downloads/fail/unclassified"),
          "unclassified", this.localReporter);

      // incremente unclassified fast5 files counter
      this.localReporter.incrCounter("numberFiles",
          "numberUnclassifiedFast5Files", numberUnclassifiedFast5Files);

      // incremente fast5 files counter
      this.localReporter.incrCounter("numberFiles", "numberFast5Files",
          numberUnclassifiedFast5Files);
    }
  }

  /**
   * This method of the class Fast5ToFastq execute the process to retrieve the
   * fastq sequence on fast5 file for the basecaller Albacore.
   * @throws IOException, test the read of the file
   */
  public void executeBasecallerAlbacore() throws IOException {

    // test if the merge of fastq is enable
    if (this.processMergeStatus) {

      // process all fast5
      this.processor.processDirectory(listAllFast5(), "merge_status",
          this.localReporter);
      return;
    }

    // test if the fail fast5 is to process
    if (this.processFail) {

      // get the list of fail fast5 files
      int numberFailFast5Files =
          processDirectory(this.fast5RunDirectory, "fail", this.localReporter);

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
          processDirectory(this.fast5RunDirectory, "pass", this.localReporter);

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
            listSubDir(this.fast5RunDirectory, "barcode");

        // process the barcode directories
        processDirectories(listBarcodeFast5Dir, this.localReporter);
      }

    }

    // test if the unclassified fast5 is to process
    if (this.processUnclassified) {

      // get the list of unclassified fast5 files
      int numberUnclassifiedFast5Files = processDirectory(
          new File(
              this.fast5RunDirectory.toPath().toString() + "/unclassified"),
          "unclassified", this.localReporter);

      // incremente unclassified fast5 files counter
      this.localReporter.incrCounter("numberFiles",
          "numberUnclassifiedFast5Files", numberUnclassifiedFast5Files);

      // incremente fast5 files counter
      this.localReporter.incrCounter("numberFiles", "numberFast5Files",
          numberUnclassifiedFast5Files);
    }

  }

}