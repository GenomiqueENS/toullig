package fr.ens.biologie.genomique.nanopore;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;

public class LogFast5toFastq {

  private final Fast5toFastq f5;
  private final Writer logWriter;

  public LogFast5toFastq(Fast5toFastq runConversionFastq, File fastqDir)
      throws IOException {
    this.f5 = runConversionFastq;
    this.logWriter =
        new FileWriter(new File(fastqDir + "/logConversionFastq.txt"));
  }

  public void createLogConversionFastq() throws IOException {
    try {
      this.logWriter.write("Log of the fast5 run "+this.f5.getNameDirectoryRunFast5()+" to extract fastq sequences\n\n");
      for (String element : f5.getListLog()){
        this.logWriter.write(element+"\n");
      }
      this.logWriter.close();
    } catch (Exception e) {
      throw new IOException(e);
    }
  }
}
