package fr.ens.biologie.genomique.toullig.trimming;

import org.usadellab.trimmomatic.fastq.FastqRecord;
import org.usadellab.trimmomatic.trim.IlluminaClippingTrimmer;
import org.usadellab.trimmomatic.trim.Trimmer;
import org.usadellab.trimmomatic.util.Logger;

import java.io.File;
import java.io.IOException;

/**
 * Class to execute the Trimmomatic trimmer. Created by birer on 27/03/17.
 */
class TrimmomaticTrimmer {

  private Trimmer trimer;

  TrimmomaticTrimmer(File adaptorFile, int seedMismatchesTrimmomatic,
      int palindromeClipThresholdTrimmomatic, int simpleClipThreshold) {

    System.out.println(adaptorFile.getPath()
        + ":" + seedMismatchesTrimmomatic + ":"
        + palindromeClipThresholdTrimmomatic + ":" + simpleClipThreshold);

    try {
      Logger logger = new Logger(true, true, true);
      this.trimer = IlluminaClippingTrimmer.makeIlluminaClippingTrimmer(logger,
          adaptorFile.getPath()
              + ":" + seedMismatchesTrimmomatic + ":"
              + palindromeClipThresholdTrimmomatic + ":" + simpleClipThreshold);

    } catch (IOException e) {
      e.printStackTrace();
    }

  }

  //
  // Trimmomatic
  //

  /**
   * Method of the class TrimmomaticTrimmer to trim a sequence with trimmomatic.
   * @param sequence, a sequence of DNA
   * @param quality, a q-quality of fastq
   * @return , a string
   */
  String TrimmomaticTrim(String sequence, String quality) {

    FastqRecord record = new FastqRecord("name", sequence, "", quality, 33);
    FastqRecord[] result =
        this.trimer.processRecords(new FastqRecord[] {record});
    return result[0].getSequence();
  }
}
