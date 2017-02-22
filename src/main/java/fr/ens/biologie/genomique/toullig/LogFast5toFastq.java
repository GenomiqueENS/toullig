package fr.ens.biologie.genomique.toullig;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.Date;

/**
 * This class make Log files.
 * @author Aurelien Birer
 */
public class LogFast5toFastq {

  private LocalReportertoLog logF5;
  private Fast5toFastq f5;
  private Writer logWriter;
  private Writer logCorruptWriter;
  private Writer logWorkflowWriter;

  /**
   * Constructor of the LogFast5toFastq class.
   * @param runConversionFastq, the run directory of fast5
   * @param fastqDir, the output directory of fastq
   * @throws IOException
   */
  public LogFast5toFastq(Fast5toFastq runConversionFastq, File fastqDir)
      throws IOException {

    this.f5= runConversionFastq;
    this.logF5 = new LocalReportertoLog(this.f5);
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
  public void createLogConversionFastq(Date dateDeb, Date dateEnd) throws IOException {
    try {
      this.logWriter.write(dateDeb+"\n\n");
      this.logWriter.write("Log of the fast5 run "
          + this.f5.getNameDirectoryRunFast5()
          + " to extract fastq sequences\n\n");
      for (String element : this.logF5.getListLog()) {
        this.logWriter.write(element + "\n");
      }
      this.logWriter.write("\n"+dateEnd);
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
      for (String element : this.logF5.getListLogStatusWorkflow()) {
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

}
