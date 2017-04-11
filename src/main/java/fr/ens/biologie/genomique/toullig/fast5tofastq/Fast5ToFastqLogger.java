package fr.ens.biologie.genomique.toullig.fast5tofastq;

import fr.ens.biologie.genomique.eoulsan.util.LocalReporter;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.*;

/**
 * This class make Log files.
 * @author Aurelien Birer
 */
public class Fast5ToFastqLogger {

  private final List<String> listWriteSequenceLog = new ArrayList<>();
  private final List<String> listWorkflowStatusLog = new ArrayList<>();

  private LocalReporter localReporter = new LocalReporter();

  private final Fast5ToFastq f5;
  private Writer logWriter;
  private Writer logCorruptWriter;
  private Writer logWorkflowWriter;

  /**
   * Constructor of the Fast5ToFastqLogger class.
   * @param runConversionFastq, the run directory of fast5
   * @param fastqDir, the output directory of fastq
   * @throws IOException if an IO error occur
   */
  public Fast5ToFastqLogger(Fast5ToFastq runConversionFastq, File fastqDir)
      throws IOException {

    // get the localReporter of the run of fast5tofastq
    this.localReporter = runConversionFastq.getLocalReporter();

    // get all the info in the list of localReporter
    getAllListForLog(localReporter);

    // get the object of the run of fast5tofastq
    this.f5 = runConversionFastq;
    this.logWriter =
        new FileWriter(new File(fastqDir + "/logConversionFastq.txt"));
    this.logCorruptWriter =
        new FileWriter(new File(fastqDir + "/logCorruptFast5Files.txt"));
    this.logWorkflowWriter =
        new FileWriter(new File(fastqDir + "/logWorkflow.txt"));
  }

  /**
   * Method to create a log about the number fastq sequence reads and write in
   * all differents status of fast5 files.
   * @throws IOException if an IO error occur
   */
  public void createLogConversionFastq(Date beginDate, Date endDate,
      List<String> arguments) throws IOException {
    try {

      // write the begin date of the fast5tofastq execution
      this.logWriter.write(beginDate + "\n\n");
      this.logWriter
          .write("Command line: " + arguments.toString().replace(",", ""));
      this.logWriter.write("Log of the fast5 run "
          + this.f5.getNameDirectoryRunFast5()
          + " to extract fastq sequences\n\n");

      // list the element of the conversion fast5 to fastq
      for (String element : getListLog()) {
        this.logWriter.write(element + "\n");
      }

      // write the end date of the fast5tofastq execution
      this.logWriter.write("\n" + endDate);
      this.logWriter.close();
    } catch (Exception e) {
      throw new IOException(e);
    }
  }

  /**
   * Method to create a log about the final message status of each common
   * workflows made by the basecalling.
   * @throws IOException if an IO error occur
   */
  public void createLogWorkflow() throws IOException {
    try {
      this.logWorkflowWriter.write(
          "\n The normal way of workflows use by basecaller follow this order : \n"
              + " - eventDetectionWorkflow \n" + " - hairpinSplitWorkflow \n"
              + " - basecall1DWorkflow \n" + " - basecall2DWorkflow \n"
              + " - calibrationStrandWorkflow \n" + " - barcodeWorkflow \n");

      // list the element of the workflow
      for (String element : getListLogStatusWorkflow()) {
        this.logWorkflowWriter.write(element + "\n");
      }
      this.logWorkflowWriter.close();
    } catch (Exception e) {
      throw new IOException(e);
    }
  }

  /**
   * Method to create a log about the complète list of fast5 files corromp find
   * in the analysis.
   * @throws IOException if an IO error occur
   */
  public void createLogCorruptFile() throws IOException {
    try {

      // write the path of corrupt fast5 file
      for (File file : this.f5.getListCorruptFast5Files()) {
        this.logCorruptWriter.write(file.toString() + "\n");
      }
      this.logCorruptWriter.close();
    } catch (Exception e) {
      throw new IOException(e);
    }
  }

  //
  //
  // Create Log
  //
  //

  /**
   * This method of the class Fast5ToFastq get a list of log on the number of
   * files.
   * @return a list of string
   */
  private List<String> getListLog() {
    List<String> listLog = new ArrayList<>();
    listLog.add("Input fast5 files read: "
        + this.f5.getNumberFast5Files(this.localReporter));
    listLog.add("Input corrupt files: "
        + this.f5.getNumberCorruptFast5Files(this.localReporter));
    listLog.add("Input calibrate strand files read: "
        + this.f5.getNumberCalibrateStrandFast5Files(this.localReporter));
    listLog.add("Input unclassified files read: "
        + this.f5.getNumberUnclassifiedFast5Files(this.localReporter));
    listLog.add("Input fail files read: "
        + this.f5.getNumberFailFast5Files(this.localReporter));
    listLog.add("Input pass files read: "
        + this.f5.getNumberPassFast5Files(this.localReporter));

    // add to the list log the write file per type
    for (String element : this.listWriteSequenceLog) {
      listLog.add(element);
    }
    return listLog;
  }

  /**
   * This method of the class Fast5ToFastq get a list of log on the status of
   * the workflows.
   * @return a list of string
   */
  private List<String> getListLogStatusWorkflow() {
    return this.listWorkflowStatusLog;
  }

  /**
   * This method of the class Fast5ToFastq get All Log for the conversion of
   * fast5 to fastq file.
   * @param localReporter a LocalReporter
   */
  private void getAllListForLog(LocalReporter localReporter) {

    //
    // Count of write sequence in fastq file
    //
    List<String> groupNameWrite = new ArrayList<>();

    // get the number of sequence write per type
    for (String element : localReporter
        .getCounterNames("numberSequenceWrite")) {
      groupNameWrite.add(element);
    }

    // sort by name
    Collections.sort(groupNameWrite);
    String oldStatusWrite = "";

    // get the element on the group name type
    for (String element : groupNameWrite) {

      long numberSequence =
          localReporter.getCounterValue("numberSequenceWrite", element);

      // test if the number of sequence equal -1 (eoulsan local reporter return)
      if (numberSequence == -1) {
        numberSequence = 0;
      }
      String[] typeSequence = element.split("Sequence")[1].split("Write");
      String[] partSplit = element.split("_");

      // test if the status is already write
      if (!oldStatusWrite.equals(partSplit[0])) {
        this.listWriteSequenceLog.add("");
      }

      // test if the status the sequence is null
      if (typeSequence[0].contains("Null")) {
        String[] partSplitNull = element.split("Nu");
        this.listWriteSequenceLog.add("Barcode "
            + partSplit[0] + ", " + partSplitNull[0]
            + " empty sequence(s) found: " + numberSequence);
      } else {
        this.listWriteSequenceLog.add("Barcode "
            + partSplit[0] + ", " + typeSequence[0] + " sequence(s) writted: "
            + numberSequence);
      }
      oldStatusWrite = partSplit[0];
    }

    // create a list of type
    List<String> groupName = new ArrayList<>();

    // get the list of type
    for (String element : localReporter.getCounterGroups()) {
      groupName.add(element);
    }

    // sort by name
    Collections.sort(groupName);
    String oldStatus = "";

    // get the element on the group name type
    for (String element : groupName) {

      String[] part = element.split("_");

      // test if the element contains a separator '_'
      if (element.contains("_")) {
        String status = part[0];

        // test if the status is not already write
        if (!oldStatus.equals(status)) {
          this.listWorkflowStatusLog
              .add("\n###################################\n"
                  + status + " :\n###################################");
        }

        this.listWorkflowStatusLog.add("\n" + part[1] + " :\n");
        stackHashWorkflow(localReporter, element, part[0]);

        oldStatus = status;
      }

    }

  }

  /**
   * This method of the class Fast5ToFastq stack the list of workflow status.
   * @param localReporter, name of object LocalReporter
   * @param group, name of group
   * @param status, the name of folder fast5 status
   */
  private void stackHashWorkflow(LocalReporter localReporter, String group,
      String status) {
    Set<String> keys = localReporter.getCounterNames(group);

    // test if the number of keys is not null
    if (keys.size() >= 1) {

      // get all key
      for (String key : keys) {

        // test if the status is fail and the key contains"Calibration strand
        // detected."
        if (status.equals("fail")
            && key.contains("Calibration strand detected.")) {
          localReporter.incrCounter("numberFiles",
              "numberCalibrateStrandFast5Files", 1);
        }
        this.listWorkflowStatusLog.add("The status \""
            + key + "\" is present " + localReporter.getCounterValue(group, key)
            + " times in the folder " + status);
      }
    } else {
      this.listWorkflowStatusLog.add(
          "No status for the " + group + " Workflow in the folder " + status);
    }
  }

}