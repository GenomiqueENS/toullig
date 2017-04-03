package fr.ens.biologie.genomique.toullig.trimming;

import fr.ens.biologie.genomique.eoulsan.bio.ReadSequence;
import fr.ens.biologie.genomique.eoulsan.bio.io.FastqWriter;
import org.usadellab.trimmomatic.fastq.FastqRecord;
import org.usadellab.trimmomatic.trim.IlluminaClippingTrimmer;
import org.usadellab.trimmomatic.util.Logger;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import static fr.ens.biologie.genomique.toullig.Utils.reverse;
import static fr.ens.biologie.genomique.toullig.trimming.UtilsTrimming.*;

/**
 * Class to execute the Trimmomatic trimmer. Created by birer on 27/03/17.
 */
public class TrimmomaticTrimmer implements Trimmer {

  private org.usadellab.trimmomatic.trim.Trimmer trimer;
  private final File nameOutputFastq;

  public TrimmomaticTrimmer(File adaptorFile, File nameOutputFastq,
      int seedMismatchesTrimmomatic, int palindromeClipThresholdTrimmomatic,
      int simpleClipThreshold) {

    this.nameOutputFastq = nameOutputFastq;

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
   * Method of the class TrimmomaticTrimmer to trimming with the Trimmer
   * interface.
   * @param leftLengthOutlier , the length of the left outlier
   * @param rightLengthOutlier , the length of the left outlier
   * @param sequence, the sequence of the read
   * @param id , the id of the read
   * @param quality , the quality of the read
   */
  public void preProcessTrimming(int leftLengthOutlier, int rightLengthOutlier,
      String sequence, String id, String quality) {

    try (BufferedWriter fastqOutputFile =
        new BufferedWriter(new FileWriter(this.nameOutputFastq))) {

      // Process for trimming
      String leftOutlierSequence =
          getOutlierLeftSequence(rightLengthOutlier, sequence);
      String rightOutlierSequence =
          getOutlierRightSequence(leftLengthOutlier, sequence);
      String leftOutlierScore =
          getOutlierLeftQuality(rightLengthOutlier, quality);
      String rightOutlierScore =
          getOutlierRightQuality(leftLengthOutlier, quality);

      String mainSequence = sequence.substring(rightLengthOutlier,
          sequence.length() - leftLengthOutlier);
      String mainScore = quality.substring(rightLengthOutlier,
          quality.length() - leftLengthOutlier);

      // Trimming with trimmomatic
      String rigthTrimSequence =
          trimmingTrimmomatic(rightOutlierSequence, rightOutlierScore);
      String leftTrimSequence = trimmingTrimmomatic(
          reverse(leftOutlierSequence), reverse(leftOutlierScore));

      // Process for writting trimming result
      String rigthTrimScore = quality.substring(mainScore.length(),
          mainScore.length() + rightOutlierSequence.length());
      String leftTrimScore = quality.substring(
          leftOutlierScore.length() - leftTrimSequence.length(),
          leftOutlierScore.length());

      String sequenceTranscript =
          leftTrimSequence + mainSequence + rigthTrimSequence;
      String scoreTranscript = leftTrimScore + mainScore + rigthTrimScore;

      ReadSequence fastq = new ReadSequence();
      fastq.setName(id);
      fastq.setSequence(sequenceTranscript);
      fastq.setQuality(scoreTranscript);
      FastqWriter fastqWriter = new FastqWriter(fastqOutputFile);
      fastqWriter.write(fastq);
    } catch (IOException e) {
      e.printStackTrace();
    }
  }

  /**
   * Method of the class TrimmomaticTrimmer to trim a sequence with trimmomatic.
   * @return , a string
   */
  String trimmingTrimmomatic(String sequence, String quality) {

    FastqRecord record = new FastqRecord("name", sequence, "", quality, 33);
    FastqRecord[] result =
        this.trimer.processRecords(new FastqRecord[] {record});
    return result[0].getSequence();
  }

  /**
   * Method of the class TrimmomaticTrimmer to trimming with the Trimmer
   * interface
   */
  public void trimming() {
  }

}
