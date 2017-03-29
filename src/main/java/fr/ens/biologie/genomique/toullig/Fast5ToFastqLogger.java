package fr.ens.biologie.genomique.toullig;

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

  private List<String> listWriteSequenceLog = new ArrayList<String>();
  private List<String> listWorkflowStatusLog = new ArrayList<String>();

  private LocalReporter localReporter = new LocalReporter();

  private Fast5ToFastq f5;
  private Writer logWriter;
  private Writer logCorruptWriter;
  private Writer logWorkflowWriter;

  /**
   * Constructor of the Fast5ToFastqLogger class.
   * @param runConversionFastq, the run directory of fast5
   * @param fastqDir, the output directory of fastq
   * @throws IOException
   */
  public Fast5ToFastqLogger(Fast5ToFastq runConversionFastq, File fastqDir)
      throws IOException {

    this.localReporter = runConversionFastq.getLocalReporter();
    getAllListForLog(localReporter);
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
   * @throws IOException
   */
  public void createLogConversionFastq(Date dateDeb, Date dateEnd)
      throws IOException {
    try {
      this.logWriter.write(dateDeb + "\n\n");
      this.logWriter.write("Log of the fast5 run "
          + this.f5.getNameDirectoryRunFast5()
          + " to extract fastq sequences\n\n");
      for (String element : getListLog()) {
        this.logWriter.write(element + "\n");
      }
      this.logWriter.write("\n" + dateEnd);
      this.logWriter.close();
    } catch (Exception e) {
      throw new IOException(e);
    }
  }

  /**
   * Method to create a log about the final message status of each common
   * workflows made by the basecalling.
   * @throws IOException
   */
  public void createLogWorkflow() throws IOException {
    try {
      this.logWorkflowWriter.write(
          "\n The normal way of workflows use by basecaller follow this order : \n"
              + " - eventDetectionWorkflow \n" + " - hairpinSplitWorkflow \n"
              + " - basecall1DWorkflow \n" + " - basecall2DWorkflow \n"
              + " - calibrationStrandWorkflow \n" + " - barcodeWorkflow \n");
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
   * @throws IOException
   */
  public void createLogCorruptFile() throws IOException {
    try {
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
  public List<String> getListLog() {
    List<String> listLog = new ArrayList<String>();
    listLog.add("Number of fast5 files read "
        + this.f5.getNumberFast5Files(this.localReporter));
    listLog.add("Number of corrupt files "
        + this.f5.getNumberCorruptFast5Files(this.localReporter));
    listLog.add("Number of calibrate strand files read "
        + this.f5.getNumberCalibrateStrandFast5Files(this.localReporter));
    listLog.add("Number of unclassified files read "
        + this.f5.getNumberUnclassifiedFast5Files(this.localReporter));
    listLog.add("Number of fail files read "
        + this.f5.getNumberFailFast5Files(this.localReporter));
    listLog.add("Number of pass files read "
        + this.f5.getNumberPassFast5Files(this.localReporter));

    listLog.add("Number of pass barcode file "
        + this.f5.getNumberBarcodeFast5Files(this.localReporter) + "\n");
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
  public List<String> getListLogStatusWorkflow() {
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
    for (String element : localReporter
        .getCounterNames("numberSequenceWrite")) {
      groupNameWrite.add(element);
    }

    Collections.sort(groupNameWrite);
    String oldStatusWrite = "";
    for (String element : groupNameWrite) {

      long numberSequence =
          localReporter.getCounterValue("numberSequenceWrite", element);
      if (numberSequence == -1) {
        numberSequence = 0;
      }
      String[] typeSequence = element.split("Sequence")[1].split("Write");
      String[] partSplit = element.split("_");
      if (!oldStatusWrite.equals(partSplit[0])) {
        this.listWriteSequenceLog.add("\n");
      }
      this.listWriteSequenceLog.add("In the file "
          + partSplit[0] + " " + typeSequence[0]
          + " the number of total sequence write " + numberSequence);
      oldStatusWrite = partSplit[0];
    }
    List<String> groupName = new ArrayList<>();
    for (String element : localReporter.getCounterGroups()) {
      groupName.add(element);
    }

    Collections.sort(groupName);
    String oldStatus = "";
    for (String element : groupName) {

      String[] part = element.split("_");
      if (element.contains("_")) {
        String status = part[0];
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

    if (keys.size() >= 1) {
      for (String key : keys) {

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
